%%
clear all % This clears all your workspace variables. It helps you to have a fresh start in your script
feature astheightlimit 2000;
%changeCobraSolver('ibm_cplex');
changeCobraSolver('gurobi');

%% loading workspace
load workspaceExample.mat

%%  Checking consistency of the original model
consistentRxnsBool = fastcc(model, 1e-4, 1);

consistentModel = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns), consistentRxnsBool))); % create a consistent model based on the vector A
fastcc(consistentModel, 1e-4, 1);

consistentModel = changeObjective(consistentModel, "biomass_reaction");

%% Discretization
fpkm = table2array(fpkm); % transform table to array
discretized = discretizeFPKM(fpkm, colnames);

%% Settings
load optionalSettings.mat
biomassRxn = 'biomass_reaction';
consensusProportion = 0.9;
epsilon = 1e-4;

%% Running rFASTCORMICS_v2
[BCRAmodel, retainedRxns] = rFastcormics(consistentModel, discretized, rownames, dico, biomassRxn, ...
    consensusProportion, epsilon, optionalSettings);

%% FBA
fbaControl = optimizeCbModel(consistentModel, 'max', 'zero');
fbaBCRA = optimizeCbModel(BCRAmodel, 'max', 'one');

%% In _silico_ gene deletions

% On the control model
[grRatioControl, grRateKOControl, ...
    grRateWTControl, ~, ~, geneListControl] = singleGeneDeletion(consistentModel, 'FBA', [], 0, 1);

% For the context-specific model:
[grRatioBCRA, grRateKOBCRA, ...
    grRateWTBCRA, ~, ~, geneListBCRA] = singleGeneDeletion(BCRAmodel, 'FBA', [], 0, 1);

%% Example with a non-sufficient medium
optionalSettingsBis = optionalSettings;
toRemove = {'glc_D[e]', 'o2[e]', 'ala_L[e]', 'arg_L[e]', 'asn_L[e]', 'asp_L[e]', 'Lcystin[e]', 'glu_L[e]', 'gly[e]', 'his_L[e]', 'ile_L[e]', 'leu_L[e]', 'lys_L[e]', 'met_L[e]', 'phe_L[e]', 'cys_L[e]', 'ser_L[e]', 'thr_L[e]'};
mask = ~ismember(optionalSettingsBis.medium, toRemove);
optionalSettingsBis.medium = optionalSettingsBis.medium(mask);

[BCRAmodelBis, retainedRxnsBis] = rFastcormics(consistentModel, discretized, rownames, dico, biomassRxn, ...
    consensusProportion, epsilon, optionalSettingsBis); %filling the medium is activated by default

%% Example without optionalSettings
[BCRAmodelTer, retainedRxnsTer] = rFastcormics(consistentModel, discretized, rownames, dico, biomassRxn, ...
    consensusProportion, epsilon);

