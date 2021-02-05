function V = LP9_4_fastcormics( K, P, model, epsilon,NonP )
%
%addpath('C:\Users\Maria\Documents\IBM\ILOG\CPLEX_Studio126')
% V = LP9( K, P, model, epsilon )
%
% CPLEX implementation of LP-9 for input sets K, P (see FASTCORE paper)

% (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
%     LCSB / LSRU, University of Luxembourg

% Maria added a reactions set that are not penalized 
scalingfactor = 1e5;

V = [];
if isempty(P) || isempty(K)
    return;
end

np = numel(P);
nk = numel(K);
[m,n] = size(model.S);
not_penalized= find(ismember(P,NonP));

% x = [v;z]

% objective
f = [zeros(1,n), ones(1,np)];
f(n+not_penalized)=0;

% equalities
Aeq = [model.S, sparse(m,np)];
beq = zeros(m,1);

% inequalities
Ip = sparse(np,n); Ip(sub2ind(size(Ip),(1:np)',P(:))) = 1;
Ik = sparse(nk,n); Ik(sub2ind(size(Ik),(1:nk)',K(:))) = 1;
Aineq = sparse([[Ip, -speye(np)]; ...
                [-Ip, -speye(np)]; ...
                [-Ik, sparse(nk,np)]]);
bineq = [zeros(2*np,1); -ones(nk,1)*epsilon*scalingfactor];

% bounds
lb = [model.lb; zeros(np,1)] * scalingfactor;
ub = [model.ub; max(abs(model.ub(P)),abs(model.lb(P)))] * scalingfactor;

x = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub);

V = x(1:n);
end
