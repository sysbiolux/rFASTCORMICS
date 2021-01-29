function [res] = GPRrulesMapper (exp, x)
    if (~( strcmp (exp , '') ) & ~isempty(exp))% verifies if there is an association rule for the each gene of the gene ID association list (no association is displayed with '' and if there is an expression value 
       
        res = eval(exp); % the association rules are stored like ‘x(25)’ meaning that the object res will take the expression value in the 25th row
       
       
    else
       res = 0; %if the association rule or the expression is missing the object takes the value 0
    end
end
 
function res = and(a,b)% in case of the reaction is associated to more genes and the association rule is given by the Boolean AND, the object takes the lowest expression value 
    res = min(a,b);
    res;
end 
 
%if the Boolean rule is an OR, the object Res takes the highest expression value
function res = or(a,b)
    res = max(a,b);
end
