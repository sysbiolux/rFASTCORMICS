function[model]=structureAnalyseFastbox(model)
% (c) Maria Pires Pacheco and Thomas Sauter 2020
%

%The function checks at the structure of the network and verifies if reversible
%reactions can carry a flux in both directions. Due to neighboring reactions
%numerous reversible reactions can carry a flux in only one direction.
%Ex: if a metabolite is connected to r reactions and r-1 reactions are
% irreversible consuming reactions then the remaining reaction is
% an irreversible producing reaction

% INPUT
%model         cobra model structure containing the fields

%   S                           m x n stoichiometric matrix
%   lb                          n x 1 flux lower bound
%   ub                          n x 1 flux upper bound
%   rxns                        n x 1 cell array of reaction abbreviations
%   c                           n x 1 

% OUTPUT(S)
%model         cobra model structure containing the fields

I_keep=[];
% find irrvesible reactions
I=find(model.lb >=0 & model.ub>0| model.ub<=0 & model.lb<0);
[model]=fixIrr_rFASTCORMICS(model);

while numel(I)> numel(I_keep)
    I=find(model.lb >=0 & model.ub>0| model.ub<=0 & model.lb<0);
    I_keep=I;    
    N=1:numel(model.rxns);
    R=setdiff(N,I);
    %creates a structure matrix all the non-zero entries are changed
    %into ones
    Struct_matrix=model.S~=0;
    Sign_matrix=model.S;
    %creates a sign matrix all the positive non-zero entries are changed
    %into ones and the negatives to minus ones
    Sign_matrix(model.S>0)=1;
    Sign_matrix(model.S<0)=-1;
    
    %all irreversible are producing
    
    I_mat=Sign_matrix(:,I);
    %find metabolites for which the irreversible reactions are all producing reactions
    prodI=find(sum(I_mat,2)==sum(abs(I_mat),2));
    
    % Find metabolites that are associated with only one reversible reaction
    NbRev=sum(Struct_matrix(:,R)~=0,2);
    ONErev=find(NbRev==1);
   % If a metabolite is connected to exactly one reversible reaction 
   % and all other connected reactions are irreversible producers, 
   % then the reversible reaction is reclassified as an irreversible 
   % consuming reaction.

    
    Met_to_change=intersect(prodI,ONErev);
    
    [~,r]=find(model.S(Met_to_change,R)<0);
    
    model.lb(R(r))=0;
    model.rev(R(r))=0;
    model.rxns(R(r));
    
    [model]=fixIrr_rFASTCORMICS(model);
    
    [~,r]=find(model.S(Met_to_change,R)>0);
 
    model.S(:,R(r))=- model.S(:,R(r));
    tmp= model.ub(R(r));
    model.ub(R(r))=-model.lb(R(r));
    model.lb(R(r))=-tmp;
    model.lb(R(r))=0;

    I=find(model.lb >=0 & model.ub>0| model.ub<=0 & model.lb<0);
    R=setdiff(N,I);
    %creates a structure matrix all the non-zero entries are changed
    %into ones
    Struct_matrix=model.S~=0;
    %creates a sign matrix all the positive non-zero entries are changed
    %into ones and the negatives to minus ones
    Sign_matrix=model.S;
    Sign_matrix(model.S>0)=1;
    Sign_matrix(model.S<0)=-1;
    
    
    I_mat=Sign_matrix(:,I);
    %find metabolites for which the irreversible reactions are all producing reactions
    comsI=find(sum(I_mat,2)==-sum(abs(I_mat),2));
    NbRev=sum(Struct_matrix(:,R)~=0,2);
    ONErev=find(NbRev==1);
   % If a metabolite is connected to exactly one reversible reaction 
   % and all other connected reactions are irreversible consumers, 
   % then the reversible reaction is reclassified as an irreversible 
   % producing reaction.
    Met_to_change=intersect(comsI,ONErev);
    [~,r]=find(model.S(Met_to_change,R)>0);
    
    model.lb(R(r))=0;
    I=find(model.lb >=0 & model.ub>0| model.ub<=0 & model.lb<0);
    R=setdiff(N,I);
    
    [~,r]=find(model.S(Met_to_change,R)<0);
    model.S(:,R(r))=- model.S(:,R(r));
    tmp= model.ub(R(r));
    model.ub(R(r))=-model.lb(R(r));
    model.lb(R(r))=-tmp;
    
    model.lb(R(r))=0;
    model.rev(R(r))=0;
    
    [model]=fixIrr_rFASTCORMICS(model);    
    I=find(model.lb >=0 & model.ub>0| model.ub<=0 & model.lb<0);    
    
end