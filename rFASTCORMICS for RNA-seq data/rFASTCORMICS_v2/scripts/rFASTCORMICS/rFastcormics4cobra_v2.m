function [contextSpecificModel, retainedRxns] = rFastcormics4cobra_v2(model, discretized, rownames, dico, consensusProportion, epsilon, optionalSettings, biomassReactionName, fillingMediumFlag, adaptiveScalingFlag)
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
% biomassReactionName:   cell array or all other functions that might require some exchange reactions to ??? (default '')/obj function ??
% fillingMediumFlag      filling the medium in case it's not sufficient for the model to fulfill the objective function. 1 for active (default), 0 for inactive
% adaptiveScalingFlag:   0 for inactive (default), 1 for active

% OUTPUT:
% contextSpecificModel        context-specific model, reduced to the retained reactions 
% Afinal:                     indices of the active reactions

% .. Authors:
%       - Maria Pires Pacheco, Thomas Sauter, 2016, University of Luxembourg
%       - Maria Pires Pacheco, Thomas Sauter, 2023, adaptation of the code to the Cobra toolbox
%       - Vanille Lejal, 2025, removing the transporters from the core,
%       removing the nonPen argument, switching to one fastcore, cleaning
%       of the code

%% Initializing the arguments
if nargin < 10
    adaptiveScalingFlag = 0;
end
if nargin < 9
    fillingMediumFlag = 1;
end
if nargin < 8
    biomassReactionName = 'biomass_reaction';
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
[consistentRxns, ~, ~] = fastcc(model, epsilon, 0, 0, 'original'); %indexes of the consistent reactions
consistentModel = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns), consistentRxns))); %removing the non consistent reactions

%% Optionally the model can be constrained in function of the medium
if ~isempty(optionalSettings) && isfield(optionalSettings, 'medium')
    %check if the medium constrain works
    if isfield(optionalSettings, 'notMediumConstrained')
        notMediumConstrained = optionalSettings.notMediumConstrained;
    else
        notMediumConstrained = [];
    end
    mediumMets = optionalSettings.medium;
    % finding the exchange reactions associated with the medium
    [excRxnsBool, ~] = findExcRxns(consistentModel);
    excRxns = consistentModel.rxns(excRxnsBool);
    [mediumRxns, ~] = findRxnsFromMets(consistentModel, mediumMets);
    excMediumRxns = intersect(excRxns, mediumRxns);

    mediumConstrainedModel =  constrainModelRFASTCORMICS(consistentModel, mediumMets, notMediumConstrained, biomassReactionName, functionKeep);
    % watch out, the model here could be unconsistent
    
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

%% Defining the core
% The discretized values will be discretized into expressed, not expressed, and unknown
% expression status.
%Core is defined as the reactions that are under the control of expressed genes
numberOfSamples = size(discretized,2);
initialCore = find(sum(mapping == 1, 2) >= (consensusProportion * numberOfSamples)); 

% Finding transporters in the Core
modelTransRxns = findTransRxns(mediumConstrainedModel);
[~, TransIDs] = ismember(modelTransRxns, mediumConstrainedModel.rxns);
% Removing transporters from the core
CoreWithoutTrans = setdiff(initialCore, TransIDs); % we remove transporters from the core

if ~isempty(functionKeep)
    B = find(ismember(mediumConstrainedModel.rxns, functionKeep)); %reactions to keep
    %problem: model has lost function keep
    if isempty(B)
        warning('No reactions to retain were found in the (medium) constrained model')
    elseif numel(B) ~= numel(functionKeep)
        warning('Not all functions set to be kept were found in the (medium) constrained model')
    end
    completedCore = union(CoreWithoutTrans, B);
else
    completedCore = CoreWithoutTrans;
end

% Get the names of the core reactions
rxnNamesCompletedCore = mediumConstrainedModel.rxns(completedCore);

%% Removing the reactions under the control of unexpressed genes
notExpressed = find(sum(mapping == -1, 2) >= (consensusProportion * numberOfSamples));
%notExpressed = setdiff(notExpressed, BiomassRelatedRxnsID); % si ça sert pour le core, on laisse même si c'est pas exprimé
mediumConstrainedModel.lb(notExpressed) = 0;
mediumConstrainedModel.ub(notExpressed) = 0;

% Building of a consistent model
[consistentRxnsAfterMedium, ~, ~] = fastcc(mediumConstrainedModel, epsilon, 0, 0, 'original');
consistentMediumConstrainedModel = removeRxns(mediumConstrainedModel, mediumConstrainedModel.rxns(setdiff(1:numel(mediumConstrainedModel.rxns), consistentRxnsAfterMedium)));

%% Checking medium sufficiency + add this loop after fastcore also
if any(strcmp(consistentMediumConstrainedModel.rxns, biomassReactionName))
    consistentMediumConstrainedModel = changeObjective(consistentMediumConstrainedModel, biomassReactionName);
    fbaResults = optimizeCbModel(consistentMediumConstrainedModel, 'max', 'zero');

    if ~isempty(fbaResults.f) && ~isnan(fbaResults.f)
        disp(['The model still contains ' biomassReactionName ' after the application of medium constraints, and FBA result is not null.' newline 'Medium is sufficient. Continuing with fastcore.']);
        needMediumFilling = false;
    else
        disp(['The model still contains ' biomassReactionName ' after the application of medium constraints, but FBA result is ' num2str(fbaResults.f) '. Medium is not sufficient.']);
        needMediumFilling = true;
    end
else
    disp(['The model lost ' biomassReactionName ' after the application of medium constraints.']);
    needMediumFilling = true;
end

%% In case medium if not sufficient
%Checking if all the exc rxns associated with the medium are in. If not, force them to be in the core
if needMediumFilling
    [excRxnsAfterMediumBool, ~] = findExcRxns(consistentMediumConstrainedModel);
    excRxnsAfterMedium = consistentMediumConstrainedModel.rxns(excRxnsAfterMediumBool);
    notPresentExcMediumRxns = setdiff(excMediumRxns, excRxnsAfterMedium);
    if ~isempty(notPresentExcMediumRxns)
        disp(['The following exchange reactions: ' newline strjoin(notPresentExcMediumRxns, newline) newline 'were initially not included in the consistent medium constraint model, even though ther were associated with the medium metabolites.' newline 'In the case of medium filling, they will be forced to be included in the model.']);
    end
end

%% Get the indices of the core reactions in the consistent constrained model
correctIndicesCompletedCore = find(ismember(consistentMediumConstrainedModel.rxns, rxnNamesCompletedCore));

%% building of the context-specific model
[contextSpecificModel, retainedRxnsBool] = fastcore(consistentMediumConstrainedModel, correctIndicesCompletedCore, 1e-4, 0, adaptiveScalingFlag);
indices = 1:numel(consistentMediumConstrainedModel.rxns);
keepRxns = indices(retainedRxnsBool == 1);
retainedRxns = find(ismember(origModel.rxns, consistentMediumConstrainedModel.rxns(keepRxns)));

%% Proceed to medium filling if needed and required
if fillingMediumFlag == 1 && needMediumFilling
    disp('Proceeding to medium filling.');
    
    % copying the bounds of rxns in the context specific model in the consistent model to preserve the medium constraints
    [~, idxConsistentModel, idxContextSpecificModel] = intersect(consistentModel.rxns, contextSpecificModel.rxns);
    consistentModel.lb(idxConsistentModel) = contextSpecificModel.lb(idxContextSpecificModel);
    consistentModel.ub(idxConsistentModel) = contextSpecificModel.ub(idxContextSpecificModel);
    
    % finding the exchange reactions associated with the medium
    [excRxnsCtxtSpeModelBool, ~] = findExcRxns(contextSpecificModel);
    excRxnsCtxtSpeModel = contextSpecificModel.rxns(excRxnsCtxtSpeModelBool);
    notPresentExcMediumRxnsSpe = setdiff(excMediumRxns, excRxnsCtxtSpeModel);
    disp(['The following exchange reactions: ' newline strjoin(notPresentExcMediumRxnsSpe, newline) newline 'will first be used to fill the medium the model as they are initially associated with the medium metabolites.']);
    
    % creating a new core with all the rxns of the contextSpecificModel,
    % the biomass, the .func rxns, and the missing medium exchange rxns
    mediumFillingCore = unique([contextSpecificModel.rxns; biomassReactionName; functionKeep; notPresentExcMediumRxnsSpe]);
    indicesMediumFillingCore = find(ismember(consistentModel.rxns, mediumFillingCore));
    
    % completing the medium and the model 
    [contextSpecificModel, retainedRxnsMediumFilledBool] = fastcore(consistentModel, indicesMediumFillingCore, 1e-4, 0, adaptiveScalingFlag);
    indices = 1:numel(consistentModel.rxns);
    keepRxns = indices(retainedRxnsMediumFilledBool == 1);
    retainedRxnsMediumFilling = find(ismember(origModel.rxns, consistentModel.rxns(keepRxns)));
    supplementaryRxns = setdiff(retainedRxnsMediumFilling, retainedRxns);
    disp(['Here are the final reactions that have been added to the model in order to fill the medium: ' newline strjoin(origModel.rxns(supplementaryRxns), newline) newline]);
    retainedRxns = retainedRxnsMediumFilling;
end

%%
disp('rFastcormics is done.');    
end
