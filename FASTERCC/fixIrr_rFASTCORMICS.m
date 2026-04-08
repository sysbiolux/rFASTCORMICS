function [model] = fixIrr_rFASTCORMICS(model)
% [model] = fixIrr_rFASTCORMICS(model)
% Fix irreversible reactions to have(lb>=0 and ub>0)
%
% INPUT
% model         cobra model structure containing the fields
% OUTPUT
% model         cobra model structure containing the fields
%
% (c) Maria Pacheco, Evelyn Gonzalez, Thomas Sauter, 2026

model.rev = zeros(numel(model.rxns),1);
model.rev(model.lb <0 & model.ub> 0) = 1;
Irr=(model.lb >=0 & model.ub>0| model.ub<=0 & model.lb<0);
model.rev(Irr) = 0;
FakeIrr= model.ub <=0 & model.lb<0;
tmp=model.ub(FakeIrr);

model.S(:,FakeIrr) = -model.S(:,FakeIrr);
model.ub(FakeIrr) = -model.lb(FakeIrr);
model.lb(FakeIrr) = -tmp;
end
