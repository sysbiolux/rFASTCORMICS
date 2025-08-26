%%
clear all % This clears all your workspace variables. It helps you to have a fresh start in your script
feature astheightlimit 2000;
%changeCobraSolver('ibm_cplex');
changeCobraSolver('gurobi');

%% loading workspace
load workspace_example.mat

%%  Checking consistency of the original model
A = fastcc(model, 1e-4, 1);

consistent_model = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns), A))); % create a consistent model based on the vector A
fastcc(consistent_model, 1e-4, 1);

consistent_model = changeObjective(consistent_model, "biomass_reaction");

%% Discretization
fpkm = table2array(fpkm); % transform table to array
discretized = discretizeFPKM(fpkm, colnames);

%% Optional settings
load optionalSettings.mat
biomassRxn = 'biomass_reaction';
consensusProportion = 0.9;
epsilon = 1e-4;

%% Running rFASTCORMICS_v2
[BCRAmodel, retainedRxns] = rFastcormics4cobra_v2(consistent_model, discretized, rownames, dico, ...
    consensusProportion, epsilon, optionalSettings, biomassRxn);

%% FBA
fba_control = optimizeCbModel(consistent_model, 'max', 'zero');
fba_BCRA = optimizeCbModel(BCRAmodel, 'max', 'zero');

%% In _silico_ gene deletions

% On the control model
[grRatio_control, grRateKO_control, ...
    grRateWT_control, ~, ~, geneList_control] = singleGeneDeletion(consistent_model, 'FBA', [], 0, 1);

% For the context-specific model:
[grRatio_BCRA, grRateKO_BCRA, ...
    grRateWT_BCRA, ~, ~, geneList_BCRA] = singleGeneDeletion(BCRAmodel, 'FBA', [], 0, 1);
