function [model, A_final] = fastcormics_2function_20170816(model, data, rownames, dico, already_mapped_tag, consensus_proportion, epsilon, optional_settings)
%(c) Maria Pires Pacheco & Thomas Sauter sep.2015 - System biology group,
% LSRU, University of Luxembourg
%Inputs
%model         cobra model structure containing the fields
%   S                           m x n stoichiometric matrix
%   lb                          n x 1 flux lower bound
%   ub                          n x 1 flux uppper bound
%   rxns                        n x 1 cell array of reaction abbreviations
%   data                        p x r matix containing 1s (expressed), 0s
%                               (unkmown), -1 (not expressed)
%                               r is the number of microarrays used as inputs for the
%                               building of 1 model
%                               p number of dataIDs (i.e. Probe IDs)
%   rownames                    p x 1 cell array dataIDs
%   already_mapped_tag          1, if the data was already to the
%                               model.rxns in this case data p = n
%                               0, if the data has to be mapped  using the
%                               GPR rules of the model
%   consensus porpotion         minimal fraction of arrays that are required
%                               to support the expression / non expression
%                               of a gene/probe ID (default 0.9, a gene
%                               is consider as expressed/non expressed
%                               if it is supported by 90% of the arrays)
%   dico                        t x 2 cell array use to convert dataIDs to the modelIDs
%                               %1st column contains the dataIDs (i.e.data_IDsS)
                                %2nd column contains the modelID
                                % TRICK if the model.gene constains a ".1"
                                %run this command
                                %model.genes=regexprep(model.genes,'\.[0-9]+$','');
                                % before runing FASTCORMICS to get rid of the
                                % "dots"
%% Default options
model_input = model;
if nargin<8
    optional_settings = '';
    if (nargin < 7)
        epsilon = 1e-4;
        if nargin <6
            consensus_proportion = 0.9;
        end
    end
end

%% map expression data to the model
number_of_array_per_model=size(data,2);
%if already_mapped_tag equals to 1,the function map_expression_2_data 
%mappes the probeIDs/gene expression levels  to the reactions via the GPR
%rules
if already_mapped_tag~=1
    mapping = map_expression_2_data(model, data,dico, rownames);
    mapping = sparse(mapping);    
else
    mapping = data;
end

if numel(rownames)~=size(data,1)
    'data and iDs do not correspond';
    return
end
if size(mapping,1)~=numel(model.rxns)
    'when the option already_mapped is used the size of the data has to correpond to the number of reaction of the model';
 return
end

%% Optionally the model can be constrained in function of the medium
%composition
if ~isempty(optional_settings) && isfield(optional_settings, 'medium')
    if isfield(optional_settings, 'not_medium_constrained')
        not_medium_constrained = optional_settings.not_medium_constrained;
    else
        not_medium_constrained = [];
    end
    medium_mets = optional_settings.medium;
    [model] =  constrain_model(model, medium_mets,not_medium_constrained);
end

%% Identifications that are under the control of expressed genes
C =  find(sum(mapping,2)>=(consensus_proportion*number_of_array_per_model)); 
% Additions of the reactions needed for a given function to carry a flux
% to the core set
if ~isempty(optional_settings)&& isfield(optional_settings, 'func') 
    B = find(ismember(model.rxns,optional_settings.func));
    if isempty(B)
        warning('no functions set to be kept')
    end
    A = fastcore_4_fastcormics(B, model, epsilon, C);
    C = union(C,A);% add the ATP and biomass reactions to the core set
end
%% Identification of the inactive reactions set (medium composition and 
%reactions under control of unexpressed genes the or  and removing of
%inactive branches
not_expressed =  find(sum(mapping,2)<=(-consensus_proportion*number_of_array_per_model));
not_expressed = setdiff(not_expressed,A);
model.lb(not_expressed) = 0;
model.ub(not_expressed) = 0;

% Building of a consistent model
[A] = fastcc_4_fastcormics(model,epsilon,0);
model_output = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),A)));

%% Establishment of the reactions in the core set 

[~, ~, IB] = intersect(model_input.rxns(C), model_output.rxns);
C = IB;
%Removal of transporters from the core set

t_keep = [];
if ~isempty(optional_settings)&& isfield(optional_settings, 'unpenalized')
    t_keep = find(ismember(model_output.rxns(C),optional_settings.unpenalized));
    t_keep = C(unique(t_keep));
end
C = setdiff(C,t_keep);
%% building of the context-specific model
A2 = fastcore_4_fastcormics(C, model_output, epsilon, t_keep);
model = removeRxns(model_output,model_output.rxns(setdiff(1:numel(model_output.rxns),A2)));
A_final = find(ismember(model_input.rxns,model_output.rxns(A2)));
end






