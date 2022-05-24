function [mapping] = map_expression_2_data_rFASTCORMICS(model, data, dico, rownames)
%% Maria Pires Pacheco & Thomas Sauter 25.03.2014 - System biology group % University of Luxembourg
% slighty adapted for dico as table input
%% search the dataIds in the dictionnary
if istable(dico)
    dico = table2array(dico);
end


if any(cellfun(@isempty, rownames))
    disp('rownames contains empty entries, please check')
    return
%     rownames(cellfun(@isempty, rownames)) = cellstr('not found');
end


col = find(sum(ismember(dico,rownames)) == max(sum(ismember(dico,rownames)))); % find matching column, usually the one with the highest number of matches
[~,idico,irownames] = intersect(dico(:,col),rownames); % get indices

if isempty(irownames);
    'the dico does match the dataIds';
    return
end

mapped = data(irownames,:); %get data for matching rownames

mapped2(:,1) = rownames(irownames,1); %create dico with same order than input
for i = 1:size(dico,2)
    try;mapped2(:,i+1) = dico(idico,i);end
    try;mapped2(:,i+1) = cellstr(dico(idico,i));end %might need conversion to cell
end

% deal with transcripts (if any)
if any(contains(model.genes,'.')) && ~contains(model.description, 'recon','IgnoreCase',true) % if possible transcripts but not recon
    warning('Does your input model contain transcripts? If not, please remove any dots . in the model.genes field')
    disp('Temporarily removing transcripts...')
    model.genes = regexprep(model.genes, '\.[0-9]*','');
elseif contains(model.description, 'recon','IgnoreCase',true) %if Recon model
    model.genes = regexprep(model.genes, '\.[0-9]*','');
end


mapped_to_genes = zeros(numel(model.genes), size(mapped,2)); %initialise variable

genes_matched = 0;

% maps the discretized expression data to the genes

for i=1:numel(model.genes);
    [match,~] = find(ismember(mapped2,model.genes(i))); %find row of match
    if numel(match)==1;
        mapped_to_genes(i,:) = mapped(match,:);
    elseif isempty(match);
    else
        mapped_to_genes(i,:) = max(mapped(match,:),[],1); % take the highest value if
        %more probeIDs correspond to one modelID
    end
    if ~isempty(match)
        genes_matched =  genes_matched + 1;
    end
end


fprintf('%i of %i genes matched\n', genes_matched, numel(model.genes))

% maps the expression data to the reactions
mapping = zeros(numel(model.rxns), size(mapped,2));
for j= 1:size(mapped_to_genes,2)
    x=mapped_to_genes(:,j);
    for k=1:numel(model.rxns);
        mapping(k,j)= GPRrulesMapper_rFASTCORMICS(cell2mat(model.rules(k)),x);
    end
end

