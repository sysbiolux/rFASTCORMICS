function [contextSpecificModel, Afinal] = rFastcormics4cobra_v2(model, discretized, rownames, dico, consensusProportion, epsilon, optionalSettings, biomassReactionName, adaptiveScalingFlag)
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

%% Initializing the arguments
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
if nargin < 5
    consensusProportion = 0.9;
end

if isfield(optionalSettings, 'func')
    functionKeep = optionalSettings.func;
else
    functionKeep = '';
end
%% Saving the original model
origModel = model;
%% Checking required solver and tooboxes
% check if it works if one solver/toolbox is missing
requiredSolvers = {'ibm_cplex', 'gurobi'};
requiredToolboxes = {'curve_fitting_toolbox', 'optimization_toolbox'};
solversPkgs = prepareTest('requireOneSolverOf', requiredSolvers, 'requiredToolboxes', requiredToolboxes);
%% Checking the format of the discretized values
if numel(rownames) ~= size(discretized, 1)
    disp('The number of gene IDs between the discretized data and the given IDs ("rownames") do not correspond')
    return
end
%% Correcting the format of the model

% checking model.subSystems format 
if ~ischar(model.subSystems{1})
    model.subSystems = vertcat(model.subSystems{:});
end

% fixing rules
if ~isfield(model, 'rules')
    disp('Creating model.rules as it was missing')
    model = generateRules(model);
end

% fixing reversibilities
model = fixIrrRFASTCORMICS(model); %need to check that

% creating a consistent model
[A, ~, ~] = fastcc(model, epsilon, 0, 0, 'original'); %indexes of the consistent reactions
consistentModel = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns), A))); %removing the non consistent reactions
%% Optionally the model can be constrained in function of the medium
if ~isempty(optionalSettings) && isfield(optionalSettings, 'medium')
    %check if the medium constrain works
    if isfield(optionalSettings, 'notMediumConstrained')
        notMediumConstrained = optionalSettings.notMediumConstrained;
    else
        notMediumConstrained = [];
    end
    mediumMets = optionalSettings.medium;

    mediumConstrainedModel =  constrainModelRFASTCORMICS(consistentModel, mediumMets, notMediumConstrained, biomassReactionName, functionKeep);
    %watch out, the model here could be unconsistent
    
    % fillingMedium shoud be added here
    
elseif isempty(optionalSettings)
    mediumConstrainedModel = consistentModel; 
    warning('No optional settings detected')
else
    mediumConstrainedModel = consistentModel;
    warning('No given medium')
end
%% Mapping the reactions to the model
mapping = mapExpressionToModel(mediumConstrainedModel, discretized, dico, rownames, 1);
mapping = sparse(mapping);

if sum(mapping) == 0
    disp('No expressed genes were mapped, check again') %no gene expressed according to the RNAseq data
    return
end

%% Discretization
% The discretized values will be discretized into expressed, not expressed, and unknown
% expression status.
%Core is defined as the reactions that are under the control of expressed genes
numberOfSamples = size(discretized,2);
initialCore = find(sum(mapping == 1, 2) >= (consensusProportion * numberOfSamples)); 

% Finding transporters in the Core
ModelTransRxns = findTransRxns(mediumConstrainedModel);
[~, TransIDs] = ismember(ModelTransRxns, mediumConstrainedModel.rxns);
% Removing transporters from the core
CoreWithoutTrans = setdiff(initialCore, TransIDs); % we remove transporters from the core

notExpressed = find(sum(mapping == -1, 2) <= (consensusProportion * numberOfSamples)); % find(sum(mapping == -1,2)>= (consensusProportion * numberOfArrayPerModel)
%notExpressed = setdiff(notExpressed, BiomassRelatedRxnsID); % si ça sert pour le core, on laisse même si c'est pas exprimé
mediumConstrainedModel.lb(notExpressed) = 0;
mediumConstrainedModel.ub(notExpressed) = 0;

% Building of a consistent model
[A, ~, ~] = fastcc(mediumConstrainedModel, epsilon, 0, 0, 'original');
consistentMediumConstrainedModel = removeRxns(mediumConstrainedModel, mediumConstrainedmodel.rxns(setdiff(1:numel(mediumConstrainedmodel.rxns), A)));

%%
if ~isempty(functionKeep)
    B = find(ismember(consistentMediumConstrainedModel.rxns, functionKeep)); %reactions to keep
    if isempty(B)
        warning('No functions set to be kept')
    elseif numel(B) ~= numel(functionKeep)
        warning('Not all functions set to be kept were found in the model')
    end
    completedCore = union(CoreWithoutTrans, B);
else
    completedCore = CoreWithoutTrans;
end
%% building of the context-specific model
[contextSpecificModel, A2] = fastcore(consistentMediumConstrainedModel, completedCore, 1e-4, 0, adaptiveScalingFlag);
indices = 1:numel(consistentMediumConstrainedModel.rxns);
keepRxns = indices(A2 == 1);

Afinal = find(ismember(origModel.rxns, consistentMediumConstrainedModel.rxns(keepRxns)));
