function [finalModel, Afinal] = rFastcormics4cobra(model, discretized, rownames, dico, consensusProportion, epsilon, optionalSettings, biomassReactionName, adaptiveScalingFlag)
% The rFASTCORMICS is a context-specific building algorithm for
% reconstructing a tissue, a cell-specific, or any context-specific model from RNAseq data
% and a generic reconstruction (Pacheco et al. 2019)

% USAGE:
%
%    A = rFastcormics4cobra(model, discretized, rownames, dico, consensusProportion, epsilon, optionalSettings, biomassReactionName)
% REQUIREMENTS:    Matlab
%                         * Statistics and Machine Learning Toolbox
%                         * Curve fitting toolbox
% INPUTS:
% model:                 object - the following fields are required - (others can be supplied)
%                         * S  - `m x 1` Stoichiometric matrix
%                         * lb - `n x 1` Lower bounds
%                         * ub - `n x 1` Upper bounds
%                         * rxns   - `n x 1` cell array of reaction abbreviations
%                         * metFormulas m*1 metabolite Formulas
% discretized:           discretized values for the samples, size(discretized,1) =
%                        number of genes, size(discretized,2)= number of
%                        samples
% rownames:              cell array with the gene IDs
% dico:                  table which contains corresponding gene identifier information. Needed
%                        to map the rownames to the genes in the model.

% OPTIONAL INPUTS:
% consensusProportion:   the rate of samples that have to express or not to express a gene for the
%                        gene to be considered expressed or not in the
%                        context of interest (default 0.9)
% epsilon:               smallest flux that is considered nonzero (default getCobraSolverParams('LP', 'feasTol')*100 = 1e-4)
% optionalSettings:      object
%                        * func - cell array of reaction abbreviations that should carry a flux
%                        * medium - cell array of metabolites abbreviation that defines metabolites
%                          in the growth medium of cells to constrain the model, see example medium_example.mat
%                        * notMediumConstrained - react ????
% biomassReactionName:   cell array or all other functions that might require some exchange reactions to ??? (default '')
% adaptiveScalingFlag:   0 for inactive(default), 1 for active

% OUTPUT:
% A:                     indices of the active reactions

% .. Authors:
%       - Maria Pires Pacheco, Thomas Sauter, 2016, University of Luxembourg
%       - Maria Pires Pacheco, Thomas Sauter, 2023, adaptation of the code to the Cobra toolbox

if nargin < 9
    adaptiveScalingFlag = 0;
end
if nargin < 8
    biomassReactionName = [];
end
if nargin < 7
    optionalSettings = [];
end
if nargin < 6 || isempty(epsilon)
    epsilon = getCobraSolverParams('LP', 'feasTol')*100;
end
if nargin<5
    consensusProportion = 0.9;
end

origModel = model;

% checking model.subSystems format 
if ischar(model.subSystems{1})
else
    model.subSystems = vertcat(model.subSystems{:});
end

% fix reversibilities
model = fixIrrRFASTCORMICS(model);

[A, ~, ~] = fastcc(model, epsilon, 0, 0, 'original');

%fix rules
if ~isfield(model, 'rules')
    disp('creating model.rules')
    model = generateRules(model);
end

% create a consistent model
consistentModel = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A)));

% Data discretization
% The discretized values will be discretized into expressed, not expressed, and unknown
% expression status.
%disp(size(discretized))
mapping = mapExpressionToModel(consistentModel, discretized, dico, rownames, 1);
%disp(mapping);
mapping = sparse(mapping);
disp(size(mapping))

if sum(mapping) == 0
    disp('no genes were mapped, check again') %no gene expressed according to the RNAseq data
    return
end

if numel(rownames) ~= size(discretized, 1)
    disp('data and iDs do not correspond')
    return
end
if size(mapping,1)~= numel(consistentModel.rxns)
    disp('when the option already_mapped is used the size of the data has to correpond to the number of reaction of the model')
    return
end

requiredToolboxes = {'curve_fitting_toolbox', 'optimization_toolbox'};
if sum(isstruct(prepareTest('requiredToolboxes', requiredToolboxes)))
%% Optionally the model can be constrained in function of the medium
%composition
if ~isempty(optionalSettings) && isfield(optionalSettings, 'medium')
    if isfield(optionalSettings, 'notMediumConstrained')
        notMediumConstrained = optionalSettings.notMediumConstrained;
    else
        notMediumConstrained = [];
    end
    mediumMets = optionalSettings.medium;
    if ~isfield(optionalSettings,'func')
        optionalSettings.func = '';
    end
    
    mediumConstrainedmodel =  constrainModelRFASTCORMICS(consistentModel, mediumMets, notMediumConstrained, biomassReactionName, optionalSettings.func);
    
    % Is the medium sufficient?
    % conditions for a medium to be sufficient : contains obj
    % function/model is consistent/FBA result is OK
    
elseif isempty(optionalSettings)
    mediumConstrainedmodel = consistentModel; 
    warning('No optional settings detected')
else
    mediumConstrainedmodel = consistentModel;
    warning('No medium set')
end

%% Identify reactions that are under the control of expressed genes
numberOfArrayPerModel = size(discretized,2);
C =  find(sum(mapping,2) >= (consensusProportion * numberOfArrayPerModel)); % find(sum(mapping == 1,2)?
%% Removing transporters from the core
%finding transporters in the Core
ModelTransRxns = findTransRxns(mediumConstrainedmodel);
[~, TransIDs] = ismember(ModelTransRxns, mediumConstrainedmodel.rxns);
C = setdiff(C, TransIDs); % we remove transporters from the core
%% Additions of the reactions needed for a given function to carry a flux
% to the core set

if isfield(optionalSettings, 'func')
    functionKeep = optionalSettings.func;
else
    functionKeep = '';
end

if ~isempty(functionKeep)
   
    B = find(ismember(mediumConstrainedmodel.rxns,functionKeep)); %reactions to keep
    if isempty(B)
        warning('no functions set to be kept')
    elseif numel(B) ~= numel(functionKeep)
        warning('Not all functions set to be kept were found in the model')
    end
    
    %disp("first fastcore");
    BiomassRelatedRxns = fastcore(mediumConstrainedmodel, B, epsilon, 0, adaptiveScalingFlag, C); % on ne pénalise pas le core mais on pénalise le reste
    %disp("end of the first fastcore");
    BiomassRelatedRxnsID = find(ismember(consistentModel.rxns,BiomassRelatedRxns.rxns)); % why consistentModel and not mediumConstrainedmodel?
    C = union(C,BiomassRelatedRxnsID); % add reactions that support the ATP and biomass reactions to the core set
    
else
    BiomassRelatedRxnsID = []; 
end

%% Identification of the inactive reactions set (medium composition and
% reactions under control of unexpressed genes the or  and removing of
% inactive branches 

% Pourquoi on fait pas cette étape dès le départ ? --> faut que ça soit
% là au départ pour éventuellement compléter le core ?
notexpressed = find(sum(mapping,2) <= (- consensusProportion * numberOfArrayPerModel)); % find(sum(mapping == -1,2)>= (consensusProportion * numberOfArrayPerModel)
% on garde les reactions qui sont not expressed si elles sont associées à la
% biomass
notexpressed = setdiff(notexpressed, BiomassRelatedRxnsID); % est-ce qu'on ne devrait pas laisser ouvertes les reactions qui sont dans le core ? même si pas expressed ?
mediumConstrainedmodel.lb(notexpressed) = 0;
mediumConstrainedmodel.ub(notexpressed) = 0;

% Building of a consistent model
[A, ~, ~] = fastcc(mediumConstrainedmodel, epsilon, 0, 0, 'original');
ConsistentmediumConstrainedmodel = removeRxns(mediumConstrainedmodel,mediumConstrainedmodel.rxns(setdiff(1:numel(mediumConstrainedmodel.rxns), A)));

%% Establishment of the reactions in the core set

[~, ~, IB] = intersect(consistentModel.rxns(C), ConsistentmediumConstrainedmodel.rxns); % on garde seulement les réactions du core qui sont dans le medium constraint model
C = IB; % indices des rxns

%% building of the context-specific model
[output_model_fastcore, A2] = fastcore(ConsistentmediumConstrainedmodel, C, 1e-4, 0, adaptiveScalingFlag);
%disp(A2);
%disp(size(A2));
indices = 1:numel(ConsistentmediumConstrainedmodel.rxns);
%disp(indices);
keep_rxn = indices(A2 == 1); % pourquoi on prend pas direct Consistentmediumconstrainedmodel.rxns(A2) ?
%disp(keep_rxn);

%
%% building of the context-specific model
finalModel = removeRxns(ConsistentmediumConstrainedmodel, ...
    ConsistentmediumConstrainedmodel.rxns(setdiff(1:numel(ConsistentmediumConstrainedmodel.rxns), keep_rxn))); % rename context-specific model ?
Afinal = find(ismember(origModel.rxns, ConsistentmediumConstrainedmodel.rxns(keep_rxn)));

end
