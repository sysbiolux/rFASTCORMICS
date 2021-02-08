function [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = DrugDeletion(model, method, DrugList)

load GeneDrugRelations.mat

if (nargin < 2)
    method = 'FBA';
end
if (nargin < 3)
    %     geneList = model.genes;
    % else
    %     if (isempty(geneList))
    %         geneList = model.genes;
    %     end
end


if ~isfield(model,'rxnGeneMat')
    %We need the reactionGeneMat during deleteModelGenes, so build it in
    %the beginning
    model = buildRxnGeneMat(model);
end

model.genes = regexprep(model.genes,'\..*','');

% model_out.genes = regexprep(model_out.genes,'\..*','')

% solWT = optimizeCbModel(model,'max','one'); % by default uses the min manhattan distance norm FBA solution.
solWT = optimizeCbModel(model,'max'); 
grRateWT = solWT.f

grRateKO = ones(numel(DrugList),1)*grRateWT;
grRatio = ones(numel(DrugList),1);
hasEffect = true(numel(DrugList),1);
fluxSolution = zeros(length(model.rxns),numel(DrugList));
delRxns = cell(numel(DrugList),1);


showprogress(0,'Single gene deletion analysis in progress ...');
for i = 1:numel(DrugList)
    showprogress(i/numel(DrugList));
    % delete all alternate transcripts
    %             delGenes = model.genes(strmatch(geneList{i},model.genes));
    
    idx = find(ismember(GxD.DrugName, DrugList(i)));
    geneList = GxD.ENTREZ(idx);
    [modelDel,hasEffect(i),constrRxnNames] = deleteModelGenes(model,geneList);

    delRxns{i} = constrRxnNames;
    if (hasEffect(i))
        switch method
            case 'lMOMA'
                solKO = linearMOMA(model,modelDel,'max');
            case 'MOMA'
                solKO = MOMA(model,modelDel,'max',false,true);
            otherwise
                solKO = optimizeCbModel(modelDel, 'max');
        end
        if (solKO.stat == 1 ||solKO.stat == 5 )
            grRateKO(i) = solKO.f;
            fluxSolution(:,i) = solKO.x;
        else
            grRateKO(i) = NaN;
        end
    end
end




grRatio = grRateKO/grRateWT;

