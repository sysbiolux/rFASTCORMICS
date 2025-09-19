function [contextSpecificModel, retainedRxns, indicesCompletedCoreOrig] = rFastcormics4cobra_v2(model, discretized, rownames, dico, consensusProportion, epsilon, optionalSettings, biomassReactionName, fillingMediumFlag, adaptiveScalingFlag)
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
% retainedRxns:               indices of the retained reactions in the input model
% rxnNamesCompletedCore       names of the 

% .. Authors:
%       - Maria Pires Pacheco, Thomas Sauter, 2016, University of Luxembourg
%       - Maria Pires Pacheco, Thomas Sauter, 2023, adaptation of the code to the Cobra toolbox
%       - Vanille Lejal, 2025, University of Luxembourg, reorganization of the code, removing the transporters from the core, switching to one fastcore,
%       integration of an option to fill a missing medium, cleaning
%       of the code, adaptation of the code for a use with gurobi and
%       Matlab 2024

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
    fprintf("Number of exchange reactions associated with the medium metabolites: %d.\n", numel(excMediumRxns));

    mediumConstrainedModel =  constrainModelRFASTCORMICS(consistentModel, mediumMets, notMediumConstrained, biomassReactionName, functionKeep);
    % watch out, the model here could be unconsistent
    
elseif isempty(optionalSettings)
    mediumConstrainedModel = consistentModel; 
    warning('No optional settings detected.')
else
    mediumConstrainedModel = consistentModel;
    warning('No given medium.')
end

%% Mapping the reactions to the model
mapping = mapExpressionToModel(mediumConstrainedModel, discretized, dico, rownames, 1);
mapping = sparse(mapping);

if sum(mapping) == 0
    disp('No expressed genes were mapped, check again.') %no gene expressed according to the RNAseq data
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
coreWithoutTrans = setdiff(initialCore, TransIDs); % we remove transporters from the core
fprintf("Number of core reactions (after mapping, transporters removed, without .func reactions): %d\n", numel(coreWithoutTrans));

if ~isempty(functionKeep)
    foundFunctionKeep = find(ismember(mediumConstrainedModel.rxns, functionKeep)); %reactions to keep
    if isempty(foundFunctionKeep)
        warning('No reactions from the .func set were found in the model.')
    elseif numel(foundFunctionKeep) ~= numel(functionKeep)
        warning('Part of the reactions from the .func set were not found in the model.')
    end
    completedCore = union(coreWithoutTrans, foundFunctionKeep);
else
    completedCore = coreWithoutTrans;
end

% Get the names of the core reactions
rxnNamesCompletedCore = mediumConstrainedModel.rxns(completedCore);
% Get the indices of the core reactions in the original model
indicesCompletedCoreOrig = find(ismember(origModel.rxns, rxnNamesCompletedCore));

%% Removing the reactions under the control of unexpressed genes
notExpressed = find(sum(mapping == -1, 2) >= (consensusProportion * numberOfSamples));
fprintf("Number of non expressed reactions (after mapping): %d.\n", numel(notExpressed));
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
    disp(['The model lost ' biomassReactionName ' after the application of medium constraints.']); % pb of display, check
    needMediumFilling = true;
end

%% In case medium if not sufficient
%Checking if all the exc rxns associated with the medium are in. If not, force them to be in the core
if needMediumFilling
    [excRxnsAfterMediumBool, ~] = findExcRxns(consistentMediumConstrainedModel);
    excRxnsAfterMedium = consistentMediumConstrainedModel.rxns(excRxnsAfterMediumBool);
    notPresentExcMediumRxns = setdiff(excMediumRxns, excRxnsAfterMedium);
    disp(notPresentExcMediumRxns);
    if ~isempty(notPresentExcMediumRxns)
        disp(['The following exchange reactions: ' newline strjoin(notPresentExcMediumRxns, newline) newline 'were initially not included in the consistent medium constraint model, even though ther were associated with the medium metabolites.' newline 'In the case of medium filling, their inclusion in the model will not be penalized.']);
    end
end

%% Get the indices of the core reactions in the consistent constrained model
correctIndicesCompletedCore = find(ismember(consistentMediumConstrainedModel.rxns, rxnNamesCompletedCore));

%% Building of the context-specific model
[contextSpecificModel, retainedRxnsBool] = fastcore(consistentMediumConstrainedModel, correctIndicesCompletedCore, 1e-4, 0, adaptiveScalingFlag);
indices = 1:numel(consistentMediumConstrainedModel.rxns);
keepRxns = indices(retainedRxnsBool == 1);
retainedRxns = find(ismember(origModel.rxns, consistentMediumConstrainedModel.rxns(keepRxns)));

%% Check if all the .func reactions are included in the model
missingFunctionKeep = functionKeep(~ismember(functionKeep, contextSpecificModel.rxns));
if ~isempty(missingFunctionKeep)
    disp(['The following .func reactions are missing:' newline strjoin(missingFunctionKeep, newline) newline ' and will be added to the context-specific model through a second fastcore run.']);
    functionKeepFlag = true;
else
    functionKeepFlag = false;
end

%% Proceed to medium filling if needed and required
if fillingMediumFlag == 1 && needMediumFilling || functionKeepFlag
    
    % we will use the consistent model as input as it still contains all the
    % reactions of the original model
    % copying the bounds of rxns in the context specific model in the consistent model to preserve the medium constraints
    [~, idxConsistentModel, idxContextSpecificModel] = intersect(consistentModel.rxns, contextSpecificModel.rxns);
    consistentModel.lb(idxConsistentModel) = contextSpecificModel.lb(idxContextSpecificModel);
    consistentModel.ub(idxConsistentModel) = contextSpecificModel.ub(idxContextSpecificModel);
    
    % Initializing a new core for a second fastcore
    fillingCoreRxns = unique([contextSpecificModel.rxns; biomassReactionName]);
    
    if fillingMediumFlag == 1 && needMediumFilling
        disp('Proceeding to medium filling.');
        % finding the exchange reactions associated with the medium
        [excRxnsCtxtSpeModelBool, ~] = findExcRxns(contextSpecificModel);
        excRxnsCtxtSpeModel = contextSpecificModel.rxns(excRxnsCtxtSpeModelBool);
        notPresentExcMediumRxnsSpe = setdiff(excMediumRxns, excRxnsCtxtSpeModel);
        if ~isempty(notPresentExcMediumRxnsSpe)
            disp(['The inclusion of the following exchange reactions: ' newline strjoin(notPresentExcMediumRxnsSpe, newline) newline 'will not be penalized to fill the model as they are initially associated with the medium metabolites.']);
        end
    else
       notPresentExcMediumRxnsSpe = ''; 
    end 
    
    if functionKeepFlag
        disp('Proceeding to .func filling.');
        fillingCoreRxns = unique([fillingCoreRxns; missingFunctionKeep]);
    end

    indicesFillingCore = find(ismember(consistentModel.rxns, fillingCoreRxns));
    indicesNonPenMediumExcRxns = find(ismember(consistentModel.rxns, notPresentExcMediumRxnsSpe));
    
    % completing the model 
    [contextSpecificModel, retainedRxnsFilledModelBool] = fastcore(consistentModel, indicesFillingCore, 1e-4, 0, adaptiveScalingFlag, indicesNonPenMediumExcRxns);
    indices = 1:numel(consistentModel.rxns);
    keepRxns = indices(retainedRxnsFilledModelBool == 1);
    retainedRxnsFilledModel = find(ismember(origModel.rxns, consistentModel.rxns(keepRxns)));
    supplementaryRxns = setdiff(retainedRxnsFilledModel, retainedRxns);
    disp(['Here are the reactions that have been added to the model during the second fastcore run in order to fill it: ' newline strjoin(origModel.rxns(supplementaryRxns), newline) newline]);
    retainedRxns = retainedRxnsFilledModel;
end

%% Uptakes that do not come from the medium
if isfield(optionalSettings, 'medium')
    [~, uptRxnsBool] = findExcRxns(contextSpecificModel);
    uptRxns = contextSpecificModel.rxns(uptRxnsBool);
    additionalMedium = setdiff(uptRxns,excMediumRxns);
    disp(['Model additional uptakes (not provided by the medium): ' newline strjoin(additionalMedium, newline) newline]);
end
%%
disp('rFastcormics is done.');    
end
