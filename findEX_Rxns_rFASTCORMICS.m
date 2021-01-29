function [exRxns, Ex_orgaInd] = findEX_Rxns(model, function_keep)
exRxnsInd       = [];
exchange_mets   = [];
exRxnsInd       = find(sum(abs(model.S),1)==1);
biomass_id      = find(ismember(model.rxns, function_keep));
if ~isempty(biomass_id)
   exRxnsInd = setdiff(exRxnsInd,biomass_id);
end
for i=1:numel(exRxnsInd);
    exmets  = find(model.S(:,exRxnsInd(i)));
    exchange_mets(end+1) = exmets;
    if model.S(exmets,exRxnsInd(i))==1;
        model.S(exmets,exRxnsInd(i))=-1;
    end
end
exRxns  = model.rxns(exRxnsInd);
model.mets(exchange_mets);
ex_mets_carbon  = (regexp(model.metFormulas(exchange_mets),'C'));
is_organic      = ~cellfun('isempty', ex_mets_carbon);
mets_ex_orga    = exchange_mets(is_organic);
mets_ex_inorga  = setdiff(exchange_mets,mets_ex_orga);
model.mets(mets_ex_orga);
model.mets(mets_ex_inorga);
Ex_orgaInd  = exRxnsInd(is_organic);
Ex_orga     = model.rxns(Ex_orgaInd);
