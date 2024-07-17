%% compare consensus models
clearvars -except solverOK, clc, close all force
load('WORKSPACE170724','model_orig','modelM','modelA','model1','AM','AA','A1','colnames','T')
sampleModelsMatrixConsensus=zeros(numel(model_orig.rxns),2);
sampleModelsMatrixConsensus(AM,1)=1;
sampleModelsMatrixConsensus(AA,2)=1;
colnamesC={'Mold','A375'}
Tmedium=T;

sampleModelsMatrix=zeros(numel(model_orig.rxns),1);
sampleModelsMatrix(A1,1)=1;
colnames=colnames(1)

%% Basic model analysis
disp('A375 vs Mold model:')
disp(' ')
disp('Number of reactions:')
disp([numel(modelA.rxns), numel(modelM.rxns)])
disp(' ')
disp('Number of metabolites:')
disp([numel(modelA.mets), numel(modelM.mets)])
disp('Number of genes:')
disp('... removing unused genes ...')
modelA=removeUnusedGenes(modelA);
disp('... removing unused genes ...')
modelM=removeUnusedGenes(modelM);
disp([numel(modelA.genes), numel(modelM.genes)])

%% Finding specific metabolites and genes (could also be done in Venny)
% and could be repeated for rxns
disp('A375 specific metabolites:')
setdiff(modelA.metNames,modelM.metNames)
disp('Mold specific metabolites:')
setdiff(modelM.metNames,modelA.metNames)

disp('A375 specific genes:')
setdiff(modelA.genes,modelM.genes)
disp('Mold specific genes:')
setdiff(modelM.genes,modelA.genes)

%% similarity between two models can be assessed via the Jaccard similarity index
altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; %shorter 10% = 1 bar

% model similarity
numel(intersect(AM,AA))/numel(union(AM,AA))

J = squareform(pdist(sampleModelsMatrixConsensus','jaccard'));
% J(isnan(J))=1;
cgo_J = clustergram(1-J,...
    'RowLabels', colnamesC,...
    'ColumnLabels', colnamesC,...
    'ColumnLabelsRotate',270, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor);
addTitle(cgo_J,{'Model similarity based on Jaccard distance'})

%% Pathway analysis
subSys=vertcat(model_orig.subSystems{:});
Pathways = table(unique(subSys));
[pathways, ~, ub] = unique(subSys);
path_counts = histc(ub, 1:length(pathways));
T = table(pathways, path_counts);
[~, ia, ib] = intersect(Pathways.Var1, T.pathways);
Pathways.consistent(ia) = T.path_counts(ib)

% Pathway information for the consensus models
[pathways, ~, ub] = unique(subSys((sampleModelsMatrixConsensus(:,1))~=0));
path_counts = histc(ub, 1:length(pathways));
T = table(pathways, path_counts);
[~, ia, ib] = intersect(Pathways.Var1, T.pathways);
Pathways.Var2(ia) = T.path_counts(ib);
Pathways.Properties.VariableNames{3} = colnamesC{1};
[pathways, ~, ub] = unique(subSys(find(sampleModelsMatrixConsensus(:,2))));
path_counts = histc(ub, 1:length(pathways));
T = table(pathways, path_counts);
[~, ia, ib] = intersect(Pathways.Var1, T.pathways);
Pathways.Var2(ia) = T.path_counts(ib) ;
Pathways.Properties.VariableNames{4} = colnamesC{2}

% pathway information for the sample-specific models
for i=1:numel(colnames)
    [pathways, ~, ub] = unique(subSys(find(sampleModelsMatrix(:,i))));
    path_counts = histc(ub, 1:length(pathways));
    T = table(pathways, path_counts);
    [~, ia, ib] = intersect(Pathways.Var1, T.pathways);
    Pathways.Var2(ia) = T.path_counts(ib) ;
    Pathways.Properties.VariableNames{4+i} = colnames{i};
end

% pathway activity rates
PathwayActivity = Pathways;
for i=3:size(PathwayActivity,2)
    PathwayActivity(:,i) = array2table(table2array(PathwayActivity(:,i))./table2array(PathwayActivity(:,2)));
end
% comparison of 2 conditions
% pathways with a difference higher than 20% 
diff_idx = find(abs(table2array(PathwayActivity(:,3))- table2array(PathwayActivity(:,4))) > 0.2);

%% plotting
figure
hold on
scatter(table2array(PathwayActivity(:,3)),table2array(PathwayActivity(:,4)),'filled',...
    'MarkerFaceColor',[0.9 0.9 0.9])
scatter(table2array(PathwayActivity(diff_idx,3)),table2array(PathwayActivity(diff_idx,4)),...
    'black')
xlabel(colnamesC(1))
ylabel(colnamesC(2))
title('Pathway presence rate in the consensus models')
line([0 1], [0,1],'Color','k')
line([0 0.8], [0.2,1],'Color','k','LineStyle','--')
line([0.2 1], [0,0.8],'Color','k','LineStyle','--')
legend({'All pathways','>20%'},"Location","best")
text(table2array(PathwayActivity(diff_idx,3)),table2array(PathwayActivity(diff_idx,4)), PathwayActivity.Var1(diff_idx))

%% comparison of multiple samples
[~,I] = sort(sum(abs(table2array(PathwayActivity(:,3:end))-mean(table2array(PathwayActivity(:,3:end)),2)),2),'descend');
cgo = clustergram(table2array(PathwayActivity(I(1:20),3:end)),...
    'RowLabels', regexprep(PathwayActivity.Var1(I(1:20)),'metabolism',''),...
    'ColumnLabels', regexprep(PathwayActivity.Properties.VariableNames(3:end),'_TCGA.*',''),...
    'ColumnLabelsRotate',270, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor);
addTitle(cgo,'Pathway activity for all models');

%% Drug screening A375
load GeneDrugRelations.mat
DrugList=unique(GeneDrugRelations.DrugName)
method='FBA'
m=modelA

idx=find(contains(m.rxns, 'biomass'))
m = changeObjective(m,m.rxns(idx));

% set medium uptake rates to concentrations
[temp,IA,IB]=intersect(m.rxns,Tmedium.medium_RPMI_EX);
tempS=nonzeros(m.S(:,IA));
for counter=1:numel(IA)
    if tempS(counter)>0
        m.ub(IA(counter))=Tmedium.Var3(IB(counter));
    else
        m.lb(IA(counter))=-Tmedium.Var3(IB(counter));
    end
end

[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = DrugDeletion(m, method, DrugList);
tabulate(grRatio)

idx=find(grRatio<0.5);
DrugList(idx)
