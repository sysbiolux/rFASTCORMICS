function [mapping] = map_expression_2_data(model, data, dico, rownames)
%% Maria Pires Pacheco & Thomas Sauter 25.03.2014 - System biology group % University of Luxembourg
% slighty adapted for dico as table input
%% search the dataIds in the dictionnary
if istable(dico)
    dico = table2array(dico);
end
% dico = lower(dico);
% rownames = lower(rownames);

col = find(sum(ismember(dico,rownames)));
[~,idico,irownames] = intersect(dico(:,col(1)),rownames);

if isempty(irownames);
    'the dico does match the dataIds';
    return
end

mapped2(:,1) = rownames(irownames,1); %data input
mapped2(:,2) = dico(idico,1); %Probe ID
mapped2(:,3) = dico(idico,2); %ENTREZ
try;mapped2(:,4) = dico(idico,3);end %HGNC
mapped = data(irownames,:);



mapped_to_genes = zeros(numel(model.genes), size(mapped,2)); %initialise variable

genes_matched = 0;

% maps the expression data to the genes
for i=1:numel(model.genes);
    match = find(strcmp(mapped2(:,2),model.genes(i))); %probe ID
    if numel(match)==1;
        mapped_to_genes(i,:)=mapped(match,:);
    elseif isempty(match);
        try; match = find(strcmp(mapped2(:,3),model.genes(i))); %ENTREZ
            if numel(match)==1;
                mapped_to_genes(i,:)=mapped(match,:);
            elseif isempty(match);
                match = find(strcmp(mapped2(:,4),model.genes(i))); %HGNC
                if numel(match)==1;
                    mapped_to_genes(i,:)=mapped(match,:);
                elseif isempty(match);
                else
                    mapped_to_genes(i,:)=max(mapped(match,:),[],1); % take the highest value if
                    %more probeIDs correspond to one modelID
                end
            else
                mapped_to_genes(i,:)=max(mapped(match,:),[],1); % take the highest value if
                %more probeIDs correspond to one modelID
            end
        end
    else
        mapped_to_genes(i,:)=max(mapped(match,:),[],1); % take the highest value if
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

