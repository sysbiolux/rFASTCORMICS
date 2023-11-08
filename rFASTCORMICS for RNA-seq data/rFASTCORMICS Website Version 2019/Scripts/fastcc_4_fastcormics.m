function A = fastcc_4_fastcormics( model, epsilon, printLevel )
% [A,V] = fastcc(model,epsilon,printLevel)
%
% The FASTCC algorithm for testing the consistency of a stoichiometric model
% Output A is the consistent part of the model
%
% INPUT
% model         cobra model structure containing the fields
%   S           m x n stoichiometric matrix
%   lb          n x 1 flux lower bound
%   ub          n x 1 flux uppper bound
%   rxns        n x 1 cell array of reaction abbreviations
%
% epsilon
% printLevel    0 = silent, 1 = summary
%
%
% OUTPUT
% A             n x 1 boolean vector indicating the flux
% (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013

tic
%number of reactions
N = (1:numel(model.rxns));
%reactions assumed to be irreversible in forward direction
I = find(model.rev==0);
not_expressed=find(model.lb==0 & model.ub==0);
A = [];
% start with I
% J is the set of irreversible reactions
J = intersect( N, I );
J= setdiff(J,not_expressed);
if printLevel==1
    fprintf('|J|=%d  ', numel(J));
end
%V is the flux vector that approximately maximizes the cardinality
V = LP7_4_rFASTCORMICS( J, model, epsilon );
Supp = find( abs(V) >= 0.99*epsilon );
%A is the set of reactions in v with absoulte value greater than epsilon
A = Supp;
if printLevel==1;
    fprintf('|A|=%d\n', numel(A));
end

%incI is the set of irreversible reactions that are flux inconsistent
incI = setdiff( J, A );
if ~isempty( incI ) & printLevel==1
    fprintf('\n(inconsistent subset of I detected)\n');
end

%J is the set of reactions with absolute value less than epsilon in V

J = setdiff( setdiff( N, A ), incI);
J= setdiff(J,not_expressed);

if printLevel==1
    fprintf('|J|=%d  ', numel(J));
end
% reversible reactions have to be tried for flux consistency in both
% directions
flipped = false;
singleton = false;
while ~isempty( J )
    if singleton
        Ji = J(1);
        V = LP3_4_rFASTCORMICS( Ji, model ) ;
    else
        Ji = J;
        V = LP7_4_rFASTCORMICS( Ji, model, epsilon ) ;
    end
    Supp = find( abs(V) >= 0.99*epsilon );
    A = union( A, Supp);
    if printLevel==1;
        fprintf('|A|=%d\n', numel(A));
    end
    if ~isempty( intersect( J, A ))
        J = setdiff( J, A );
        if printLevel==1
            
            fprintf('|J|=%d  ', numel(J));
        end
        flipped = false;
    else
        JiRev = setdiff( Ji, I );
        if flipped || isempty( JiRev )
            flipped = false;
            if singleton
                J = setdiff( J, Ji );
                if printLevel==1
                    
                    fprintf('\n(inconsistent reversible reaction detected)\n');
                end
            else
                singleton = true;
            end
        else
            model.S(:,JiRev) = -model.S(:,JiRev);
            tmp = model.ub(JiRev);
            model.ub(JiRev) = -model.lb(JiRev);
            model.lb(JiRev) = -tmp;
            flipped = true;
            if printLevel==1
                
                fprintf('(flip)  ');
            end
        end
    end
end

if numel(A) == numel(N)
    fprintf('\nThe input model is consistent.\n');
end
toc
end
