function [model] = constrainModelRFASTCORMICS(model, mediumMets, notMediumConstrained, biomassReaction, functionToKeep)
% (c) Maria Pires Pacheco 2015
% Close all the exchange reactions which are carbon sources except the
% biomass reaction, the reactions that are supplying the medium or those
% that are in the functionToKeep input

%identification of the exchange reactions ExRxns
if ~isempty(notMediumConstrained)
    notMediumConstrained = find(ismember(model.rxns,notMediumConstrained)); %exchange reactions that we don't want to constrain
    %keeping into memory the bound values of the reaction that we do not
    %want to constrain
    lb = model.lb(notMediumConstrained);
    ub = model.ub(notMediumConstrained);
end
[exRxns, ExOrgaInd] = findEXRxnsRFASTCORMICS(model, biomassReaction, functionToKeep); % Retrieves the exchange reactions to be closed (excluding biomassReaction and functionToKeep), as well as the carbon sources
notConstrained = intersect(findRxnsFromMets(model,mediumMets),exRxns); % Exchange reactions supplied by the medium
notConstrained = find(ismember(model.rxns,notConstrained)); % Exchange reactions not to be constrained because they are supplied by the medium --> all others will be closed
lb2 = model.lb(notConstrained);
ub2 = model.ub(notConstrained);
[~,match]= find(model.S(:,ExOrgaInd)<0);
if ~isempty(match)
    model.lb(ExOrgaInd(match)) = 0; % close the uptake of all carbon sources
end
[~,match]= find(model.S(:,ExOrgaInd)>0);
if ~isempty(match)
    model.ub(ExOrgaInd(match)) = 0; % close the secretion of all the carbon sources
end
if ~isempty(notMediumConstrained)
    % Restore the bounds of the reactions that we do not want to constrain
    model.lb(notMediumConstrained) = lb;
    model.ub(notMediumConstrained) = ub;
end
if  ~isempty(notConstrained)
    model.lb(notConstrained) = lb2;
    model.ub(notConstrained) = ub2;
end
end
