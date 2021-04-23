function[model] = fixIrr(model)
% Maria Pires Pacheco
model.rev = zeros(numel(model.rxns),1);
model.rev(model.lb <0 & model.ub> 0) = 1;
Irr=(model.lb >=0 & model.ub>0| model.ub<=0 & model.lb<0);
model.rev(Irr) = 0;
FakeIrr= model.ub <=0 & model.lb<0;
model.S(:, FakeIrr) = -model.S(:,FakeIrr);
model.ub(FakeIrr) = -model.lb(FakeIrr);
model.lb(FakeIrr) = zeros(sum(FakeIrr),1);
end