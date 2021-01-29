function [model] = constrain_model(model, medium_mets, not_medium_constrained, function_keep)
% (c) Maria Pires Pacheco 2015

%identification of the exchange reactions ExRxns
if ~isempty(not_medium_constrained)
    not_medium_constrained = find(ismember(model.rxns,not_medium_constrained));
    lb = model.lb(not_medium_constrained);
    ub = model.ub(not_medium_constrained);
end
[exRxns, Ex_orgaInd] = findEX_Rxns_rFASTCORMICS(model, function_keep);
not_constrained = intersect(findRxnsFromMets(model,medium_mets),exRxns);
not_constrained = find(ismember(model.rxns,not_constrained));
lb2 = model.lb(not_constrained);
ub2 = model.ub(not_constrained);
model.lb(Ex_orgaInd) = 0; % close all the carbon sources
if ~isempty(not_medium_constrained)
    model.lb(not_medium_constrained) = lb;
    model.ub(not_medium_constrained) = ub;
end
if  ~isempty(not_constrained);
    model.lb(not_constrained) = lb2;
    model.ub(not_constrained) = ub2;
end
end
