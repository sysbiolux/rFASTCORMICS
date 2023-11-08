function [mapping, mapped_to_genes] = map_expression_2_data_v2(model, data, dico, rownames)
%% Maria Pires Pacheco & Thomas Sauter 25.03.2014 - System biology group % University of Luxembourg
% slighty adapted for dico as table input
%% search the rownames in the dictionnary
col = find(sum(ismember(table2array(dico),rownames)));
[~,IA,loc] = intersect(table2array(dico(:,col)),rownames);

if isempty(loc);
    'the dico does match the dataIds';
    return
end

mapped2(:,1) = rownames(loc,1);
mapped2(:,2) = table2array(dico(IA,1));
mapped2(:,3) = table2array(dico(IA,2));
mapped = data(loc,:);

mapped_to_genes = zeros(numel(model.genes), size(mapped,2)); %initialise variable

% maps the expression data to the genes
for i=1:numel(model.genes);
    match = find(strcmp(mapped2(:,3),model.genes(i)));
    if numel(match)==1;
        mapped_to_genes(i,:)=mapped(match,:);
    elseif isempty(match);
        
    else
        mapped_to_genes(i,:)=max(mapped(match,:),[],1); % take the highest value if
        %more probeIDs correspond to one modelID
    end
end

% maps the expression data to the reactions
mapping = zeros(numel(model.rxns), size(mapped,2));
for j= 1:size(mapped_to_genes,2)
    x=mapped_to_genes(:,j);
    for k=1:numel(model.rxns);
        mapping(k,j)= GPRrulesMapper(cell2mat(model.rules(k)),x);
    end
end

end


