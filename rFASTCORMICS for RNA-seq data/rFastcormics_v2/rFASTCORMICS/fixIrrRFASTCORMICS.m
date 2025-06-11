function[model] = fixIrrRFASTCORMICS(model)
% The function converts irreversible backwards reactions into irreversible
% forward reactions

% USAGE:
%
%   [model] = fixIrrRFASTCORMICS(model)

%INPUTS:
%    model:             (the following fields are required - others can be supplied)
%                         * S  - `m x 1` Stoichiometric matrix
%                         * lb - `n x 1` Lower bounds
%                         * ub - `n x 1` Upper bounds
%                         * rxns   - `n x 1` cell array of reaction
%                         abbreviations

% OUTPUT:                 
% model:                 model with corrected reversibilties

% .. Authors:
%       - Maria Pires Pacheco, Thomas Sauter, 2016, University of Luxembourg
%       - Maria Pires Pacheco, Thomas Sauter, 2022, adaptation of the code to the Cobra toolbox


model.rev = zeros(numel(model.rxns),1); %on prend toutes les reactions du modèle et on assigne 0
model.rev(model.lb <0 & model.ub> 0) = 1; %on assigne un 1 aux réactions réversibles 
Irr =(model.lb >=0 & model.ub >0| model.ub <=0 & model.lb <0); %les reactions irréversibles sont soient celles avec un flux que négatif ou que positif
model.rev(Irr) = 0; %on assigne 0 à ces réactions
FakeIrr= model.ub <=0 & model.lb<0; %réactions qui peuvent avoir un flux que négatif sont définies comme les fakeIrr
model.S(:, FakeIrr) = -model.S(:,FakeIrr); %ces réactions là on change leur signe dans la matrice S
model.ub(FakeIrr) = -model.lb(FakeIrr); %on change le signe de ces réactions
model.lb(FakeIrr) = zeros(sum(FakeIrr),1);
end