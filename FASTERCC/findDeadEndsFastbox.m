function[ model]=findDeadEndsFastbox(model)
%% Tutorial and sanity check of the  compareModel functions
%
% model                 = cobra model

%% closed reactions
modelOrg=model;
deadEndsI=zeros(numel(model.mets),1);
deadEndsII=zeros(numel(model.mets),1);

%Identify closed reactions and remove them
closed= model.lb==0 & model.ub==0;
if sum(closed)>0
    [model]=removeRxns(model,model.rxns(closed));
end
k=0;
%%

deadendmets=model.mets;
dead_end_mets_II=model.mets;
%%
while (numel(deadendmets)>0 || numel(dead_end_mets_II)>0)

    [model]=structureAnalyseFastbox(model);
    [model] = fixIrr_rFASTCORMICS(model);
    %% this part identify deadends type I and remove them
    k=k+1;
    if k==1 % only to run in the 1st run
       
        Struct_matrix=model.S~=0;
        % find deadends metabolites (metabolites only connected to 1
        % reaction, class I)
        dead_end_mets_I=sum(Struct_matrix,2)==1;
        deadendmets=model.mets(dead_end_mets_I);% name of the deadends
        deadEndsI(ismember(modelOrg.mets,deadendmets))=1;
    end
    %remove the reactions associated to type I deadends
    [~,r]=find(model.S(ismember(model.mets,deadendmets),:));
    model= removeRxns(model, model.rxns(r));
    
    
    
    
    
    %%     Then identify dead end type II (deadends only connected to comsuming or only producing reactions)
    % compute the sign matrix (negative coefficients are converted
    % to -1 and positives coefficients to +1
    if k==1
        Struct_matrix=model.S~=0;
        Sign_matrix=model.S;
        Sign_matrix(model.S>0)=1 ;
        Sign_matrix(model.S<0)=-1 ;
        % all irreversible reactions
        I=(model.lb >=0 & model.ub>0| model.ub<=0 & model.lb<0);
        %metabolites that are connected only to irreversible reactions
        mets_linked_only_Irr=find(sum(Struct_matrix,2)== sum(Struct_matrix(:,I),2));
        %Blocked or exchange reactions
        incont_and_exchange_mets=mets_linked_only_Irr(abs(sum(Sign_matrix(mets_linked_only_Irr,:),2))==sum(Struct_matrix(mets_linked_only_Irr,I),2));
        ExRxns=sum(model.S~=0,1)==1;
        ExMets =find(sum(abs(model.S(:,ExRxns)),2));
        
        dead_end_mets_II=setdiff(incont_and_exchange_mets,ExMets);
        dead_end_mets_II=model.mets(dead_end_mets_II);
        deadEndsII(ismember(modelOrg.mets,deadendmets))=1;
    end
    %remove the reactions associated to type II deadends
    [~,r]=find(model.S(ismember(model.mets,dead_end_mets_II),:));
    model= removeRxns(model, model.rxns(r));
    
    % check if new dead end I appeared
    Struct_matrix=model.S~=0;
    dead_end_mets_I=sum(Struct_matrix,2)==1;
    deadendmets=model.mets(dead_end_mets_I);% name of the deadends
     deadEndsI(ismember(modelOrg.mets,deadendmets))=1;

    
    % check if new deadend II appeared
    Sign_matrix=model.S;
    Sign_matrix(model.S>0)=1 ;
    Sign_matrix(model.S<0)=-1 ;
    % all irreversible reactions
    I=(model.lb >=0 & model.ub>0| model.ub<=0 & model.lb<0);
    % For metabolites that are connected only to irreversible reactions:  
    % - The left-hand value is the total number of reactions connected to the metabolite.  
    % - The right-hand value is the number of irreversible reactions connected to it.  
    % If both values are equal, the metabolite has no reversible reactions.
    mets_linked_only_Irr=find(sum(Struct_matrix,2)== sum(Struct_matrix(:,I),2));
    % Blocked or exchange reactions:  
    % - Metabolites connected to only one reaction may represent exchange reactions.  
    % - Dead-end metabolites are those connected to only one reaction that is not an exchange reaction.

    incont_and_exchange_mets=mets_linked_only_Irr(abs(sum(Sign_matrix(mets_linked_only_Irr,:),2))==sum(Struct_matrix(mets_linked_only_Irr,I),2));
    ExRxns=sum(model.S~=0,1)==1;
    ExMets =find(sum(abs(model.S(:,ExRxns)),2));
    
    dead_end_mets_II=setdiff(incont_and_exchange_mets,ExMets);
    dead_end_mets_II=model.mets(dead_end_mets_II);
    
    deadEndsII(ismember(modelOrg.mets,dead_end_mets_II))=1;

    
end




end
