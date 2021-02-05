% cd('\\atlas\users\tamara.bintener\PhD Project\DATA\TCGA')

tumour = readtable('C:\Users\tamar\Desktop\Home Office\Book Chapter data\GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt');
normal = readtable('C:\Users\tamar\Desktop\Home Office\Book Chapter data\GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FPKM.txt');

tumour_names = readtable('\\atlas\users\tamara.bintener\PhD Project\DATA\TCGA\GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt','ReadVariableNames',false);
tumour_names.Properties.VariableNames = {'Name' 'Cancer'};
tumour_names.Name = regexprep(tumour_names.Name, '-','_');
tumour_names2 = strcat(strcat(strcat(strcat(tumour_names.Cancer','_'),'C'),'_')',tumour_names.Name);


normal_names = readtable('\\atlas\users\tamara.bintener\PhD Project\DATA\TCGA\GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples','ReadVariableNames',false);
normal_names.Properties.VariableNames = {'Name' 'Cancer'};
normal_names.Name = regexprep(normal_names.Name, '-','_');
normal_names2 = strcat(strcat(strcat(strcat(normal_names.Cancer','_'),'H'),'_')',normal_names.Name);

rownames = tumour.Var1;

tumour.Var1 = [];
normal.Var1 = [];


[ia1, ib1] = ismember(tumour.Properties.VariableNames, tumour_names.Name);
[tumour.Properties.VariableNames', tumour_names.Name(ib1)];

[ia2, ib2] = ismember(normal.Properties.VariableNames, normal_names.Name);
[normal.Properties.VariableNames', normal_names.Name(ib2)];


fpkm = [tumour, normal];
colnames = [tumour_names2(ib1); normal_names2(ib2)];

[fpkm.Properties.VariableNames',colnames];

fpkm.Properties.VariableNames = colnames;


BRCA_C = find(contains(colnames, 'BRCA_C'));
BRCA_H = find(contains(colnames, 'BRCA_H'));


% colnames(BRCA_C(randperm(numel(BRCA_C),10)));
% colnames(BRCA_H(randperm(numel(BRCA_H),10)));



fpkm_BRCA_C = fpkm(:, BRCA_C(randperm(numel(BRCA_C),10)));
fpkm_BRCA_H = fpkm(:, BRCA_H(randperm(numel(BRCA_H),10)));


fpkm_BRCA_C.Properties.RowNames = rownames;
fpkm_BRCA_H.Properties.RowNames = rownames;

writetable(fpkm_BRCA_C,'fpkm_BRCA_cancer.txt','Delimiter','\t','WriteRowNames',true)
writetable(fpkm_BRCA_H,'fpkm_BRCA_control.txt','Delimiter','\t','WriteRowNames',true)
