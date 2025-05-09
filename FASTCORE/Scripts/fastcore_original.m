function A = fastcore( C, model, epsilon ) 
%
% A = fastcore( C, model, epsilon )
%
% The FASTCORE algorithm for context-specific metabolic network reconstruction
% Input C is the core set, and output A is the reconstruction

% (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
%     LCSB / LSRU, University of Luxembourg

tic
model_org = model;

%% Check input model
%model.rev field
if ~isfield(model,'rev')
    disp('creating model.rev')
    model.rev = zeros(numel(model.lb),1);
    model.rev(model.lb<0) = 1;
end

% unnest subsystems
if length(model.subSystems{1}) == 1
    disp('unnesting subsystems')
    model.subSystems = vertcat(model.subSystems{:});
end

%fix irreversible reactions
model = fixIrr_FASTCORE(model);

%fix rules
if ~isfield(model, 'rules')
    disp('creating model.rules')
    model = generateRules_FASTCORE(model);
end

%% FASTCORE

N = 1:numel(model.rxns);
I = find(model.rev==0);

A = [];
flipped = false;
singleton = false;  

% start with I
J = intersect( C, I ); fprintf('|J|=%d  ', length(J));
P = setdiff( N, C);
Supp = findSparseMode_fastcore( J, P, singleton, model, epsilon );
if ~isempty( setdiff( J, Supp ) ) 
  fprintf ('Error: Inconsistent irreversible core reactions.\n');
  return;
end
A = Supp;  fprintf('|A|=%d\n', length(A));
J = setdiff( C, A ); fprintf('|J|=%d  ', length(J));

% main loop     
while ~isempty( J )
    P = setdiff( P, A);
    Supp = findSparseMode_fastcore( J, P, singleton, model, epsilon );
    A = union( A, Supp );   fprintf('|A|=%d\n', length(A)); 
    if ~isempty( intersect( J, A ))
        J = setdiff( J, A );     fprintf('|J|=%d  ', length(J));
        flipped = false;
    else
        if singleton
            JiRev = setdiff(J(1),I);
        else
            JiRev = setdiff(J,I);
        end
        if flipped || isempty( JiRev )
            if singleton
                fprintf('\nError: Global network is not consistent.\n');
                return
            else
              flipped = false;
              singleton = true;
            end
        else
            model.S(:,JiRev) = -model.S(:,JiRev);
            tmp = model.ub(JiRev);
            model.ub(JiRev) = -model.lb(JiRev);
            model.lb(JiRev) = -tmp;
            flipped = true;  fprintf('(flip)  ');
        end
    end
end
fprintf('|A|=%d\n', length(A));

toc

