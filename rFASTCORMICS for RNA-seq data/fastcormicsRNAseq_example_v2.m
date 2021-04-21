%% FASTCORMICS RNA-seq
% _(c) Dr. Maria Pires Pacheco 2016_
% 
% _Example script adapted by Tamara Bintener_
%% Disclaimer
% In order to run FASTCORMICS RNA-seq you will need
%% 
% * Matlab
% * Compatible Cplex version, added to the Matlab path
% * Curve fitting toolbox
% * Cobra Toolbox installed (https://opencobra.github.io/cobratoolbox/latest/installation.html)
%% 
% The example data was downloaded from <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1697009 
% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1697009> (FPKM values 
% of normal tissue). The first two samples (TCGA06067511A32RA36H07 and TCGA06067811A32RA36H07) 
% will serve as an example was already imported, split into 3 variables, and saved 
% under _example_FPKM.mat ._
%% 
% * colnames:   cell array with the sample names (here TCGA06067511A32RA36H07)
% * rownames:  cell aray with the gene IDs
% * fpkm:           fpkm values for the samples, size(fpkm) = lenghth(geneIDs) 
% length(colnames)
%% Setup

addpath(genpath(pwd)) % add all subfolders to path

load('example_FPKM.mat') % load the data

figflag = 1; % set figure flag: 1 to output and save density plot
%% Data discretization

discretized = discretize_FPKM(fpkm, colnames) % no figures
%% FASTCORMICS
%% Prepare FASTCORMICS
% needed:
%% 
% * Cmodel:     a consistent model, here: obtained from running fastcc on Recon204 
% (from https://vmh.uni.lu/#downloadview)
% * dico:   identifier file, which links the geneIDs from the RNA-seq data to 
% the gene identifiers in the model
%% 
% optional
%% 
% * *medium*:     defines metabolites in the growth medium of cells to constrain 
% the model, see example medium_example.mat

load Recon2.v04.mat  % Recon 2.04 model, 
A = fastcc_4_rfastcormics(modelR204, 1e-4,0) % create consistent model by running FASTCC (Vlassis et al., 2014)
Cmodel = removeRxns(modelR204, modelR204.rxns(setdiff(1:numel(modelR204.rxns),A)))

load medium_example % need to define medium for the cells used here
load dico_ML.mat % dictionary to map the rownname identifier to the genes in the model

Cmodel.genes = regexprep(Cmodel.genes,'\.[0-9]+$',''); %removal of gene transcripts

epsilon = 1e-4;
consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one generic model from different samples
%% 
% Set optional settings such as:
%% 
% * unpenalizedSystems:      
% * func: reaction(s) forced to be present in the model
% * not_medium_constrained
% * medium : medium composition 

unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};
unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));
optional_settings.unpenalized = unpenalized;

optional_settings.func = {'DM_atp_c_';'biomass_reaction'}; % forced reactions

not_medium_constrained = 'EX_tag_hs(e)';
optional_settings.not_medium_constrained = not_medium_constrained;

optional_settings.medium = medium_example;% remove field if no constraint is provided
%% Create models
% the cobra toolbox is needed for the model creation
% 
% make sure you have a compatible cplex version for your Matlab version

% initCobraToolbox % initialize the COBRA toolbox if needed
%% 
% single models:

for i = 1:numel(colnames) %for each sample
    [model_out{i}, A_keep{i}] = fastcormics_RNAseq(Cmodel, discretized(:,i), ...
        rownames, dico , 0, consensus_proportion, epsilon, optional_settings);
end
%% 
% generic models:

[model_out_generic, A_keep_generic] = fastcormics_RNAseq(Cmodel, discretized, ...
    rownames, dico , 0, consensus_proportion, epsilon, optional_settings);