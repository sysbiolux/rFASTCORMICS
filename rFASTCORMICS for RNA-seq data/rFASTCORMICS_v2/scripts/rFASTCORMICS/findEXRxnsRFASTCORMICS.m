function [exRxns, ExOrgaInd] = findEXRxnsRFASTCORMICS(model, biomassReaction, functionToKeep)
% (c) Maria Pires Pacheco 2015
% Close all the exchange reactions which are carbon sources except the
% biomass reaction, the reactions that are supplying the medium or those
% that are in the functionToKeep input

exchangeMets = []; % exchange metabolites
exRxnsInd = find(sum(abs(model.S),1)==1); % Exchange reactions only have one non-zero (+1 / -1) element in the corresponding column of the stoichiometric matrix.
biomassId = find(ismember(model.rxns, biomassReaction)); % index of the biomass reaction(s) in the objective function

if isempty(biomassId) && isempty(functionToKeep) % if no biomass and no functionToKeep -> exRxnsInd = all exchange rxns
    warning('No biomass set for the model, please verify if the medium constraints do not affect the biomass production')
elseif isempty(biomassId) && ~isempty(functionToKeep) % if biomass is empty but functionToKeep is provided
    warning('No biomass set for the model, please verify if the medium constraints do not affect the biomass production')
    functionId = find(ismember(model.rxns, functionToKeep)); % ids of the reactions to keep (those listed in functionToKeep)
    exRxnsInd = setdiff(exRxnsInd, functionId); % ids of exchange reactions that are not in the reactions to keep
else % if both biomass and functionToKeep are provided
    % what about the case where biomass is provided but functionToKeep is empty?
    functionId = find(ismember(model.rxns, functionToKeep)); % ids of the reactions to keep (those listed in functionToKeep)
    functionId = unique([biomassId; functionId]); % keep both the biomass reaction ids and the functionToKeep reactions
    exRxnsInd = setdiff(exRxnsInd, functionId); % exchange reactions are those not in the ids to keep
end
% we remove from exchange reactions those included in biomass or functionToKeep

for i = 1:numel(exRxnsInd)
    exmets = find(model.S(:,exRxnsInd(i)));
    exchangeMets(end+1) = exmets; % progressively add the exchange metabolites
    if model.S(exmets,exRxnsInd(i)) == 1
        model.S(exmets,exRxnsInd(i)) = -1; % in case the exchange reaction (import) was defined as reversible, we force it to irreversible import (RECON convention?)
    end % what is this loop for???
end

exRxns  = model.rxns(exRxnsInd);
model.mets(exchangeMets); % not used

exMetsX  = (regexp(model.metFormulas(exchangeMets),'X')); % why X and Y metabolites?
exMetsY  = (regexp(model.metFormulas(exchangeMets),'Y')); % used to check whether exchange reactions involving X or Y metabolites are closed or not. Useful? Should be passed as argument?

if ~isempty(model.metFormulas(exchangeMets(~cellfun('isempty', exMetsX))))
    disp('Warning metabolites with X in their Formulas these inputs are closed')
end
if ~isempty(model.metFormulas(exchangeMets(~cellfun('isempty', exMetsY))))
    disp('Warning metabolites with Y in their Formulas these inputs are closed')
end

exMetsCarbon  = (regexp(model.metFormulas(exchangeMets),'C')); % exchange reactions involving carbon
exMetsR  = (regexp(model.metFormulas(exchangeMets),'R')); % exchange reactions involving R - what is R?
exKnowInorganic = (ismember(model.metFormulas(exchangeMets),'Ca') | ismember(model.metFormulas(exchangeMets),'Cl') | ismember(model.metFormulas(exchangeMets),'Co') | ismember(model.metFormulas(exchangeMets),'Cu'));
model.mets(exchangeMets(exKnowInorganic)); % print exchange reactions that are inorganic (non-nutritive, e.g., O2)

isOrganic = (~cellfun('isempty', exMetsCarbon) | ~cellfun('isempty', exMetsX) | ~cellfun('isempty', exMetsY) | ~cellfun('isempty', exMetsR)) & ~exKnowInorganic; % list of organic exchange reactions

% not used
metsExOrganic    = exchangeMets(isOrganic);
metsExInorganic  = setdiff(exchangeMets, metsExOrganic);
model.mets(metsExOrganic);
model.mets(metsExInorganic);

ExOrgaInd  = exRxnsInd(isOrganic); % output = indices of organic exchange reactions
%Ex_orga     = model.rxns(ExOrgaInd); % not used
end
