function [exRxns, Ex_orgaInd] = findEX_Rxns_rFASTCORMICS(model, biomass_rxn, function_keep)
exRxnsInd       = [];
exchange_mets   = [];
exRxnsInd       = find(sum(abs(model.S),1)==1);
biomass_id      = find(ismember(model.rxns, biomass_rxn));

if isempty(biomass_id) && isempty(function_keep)
    warning('No biomass set for the model, please verify if the medium constraints do not affect the biomass production')
elseif isempty(biomass_id) && ~isempty(function_keep)
    warning('No biomass set for the model, please verify if the medium constraints do not affect the biomass production')
    function_id = find(ismember(model.rxns, function_keep));
    exRxnsInd = setdiff(exRxnsInd, function_id);
else
    function_id = find(ismember(model.rxns, function_keep));
    function_id = unique([biomass_id; function_id]);
    exRxnsInd = setdiff(exRxnsInd, function_id);
end

for i=1:numel(exRxnsInd)
    exmets  = find(model.S(:,exRxnsInd(i)));
    exchange_mets(end+1) = exmets;
    if model.S(exmets,exRxnsInd(i))==1
        model.S(exmets,exRxnsInd(i))=-1;
    end
end
exRxns  = model.rxns(exRxnsInd);
model.mets(exchange_mets);

ex_mets_X  = (regexp(model.metFormulas(exchange_mets),'X')) ;
ex_mets_Y  = (regexp(model.metFormulas(exchange_mets),'Y')) ;

if ~isempty(model.metFormulas(exchange_mets(~cellfun('isempty', ex_mets_X))))
    disp('Warning metaboltes with X in their Formulas these inputs are close if no carbon is present')
end
if ~isempty(model.metFormulas(exchange_mets(~cellfun('isempty', ex_mets_Y))))
    disp('Warning metaboltes with Y in their Formulas these inputs are close if no carbon is present')
end
ex_mets_carbon  = (regexp(model.metFormulas(exchange_mets),'C')) ;
ex_mets_R  = (regexp(model.metFormulas(exchange_mets),'R')) ;
ex_know_inorganic=(ismember(model.metFormulas(exchange_mets),'Ca') |ismember(model.metFormulas(exchange_mets),'Cl')  |ismember(model.metFormulas(exchange_mets),'Co')  | ismember(model.metFormulas(exchange_mets),'Cu'))  
model.mets(exchange_mets(ex_know_inorganic))

is_organic      = ~cellfun('isempty', ex_mets_carbon)| ~cellfun('isempty', ex_mets_R)& ~ex_know_inorganic;

mets_ex_orga    = exchange_mets(is_organic);
mets_ex_inorga  = setdiff(exchange_mets,mets_ex_orga);
model.mets(mets_ex_orga);
model.mets(mets_ex_inorga);
Ex_orgaInd  = exRxnsInd(is_organic);
Ex_orga     = model.rxns(Ex_orgaInd);
end
