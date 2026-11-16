% toy_spectra_fastcc.m -- a toy model where a reaction-level epsilon threshold removes a consumer but keeps its producer.
%
% Network (all reactions irreversible, default bounds 0..1000, uptake v1 <= 1e-4, epsilon = 1e-4):
%   v1:  -> A[e]                v2:  A[e] -> A[c]
%   B branch:  v3: A[c] -> B[c]     v4: B[c] -> 0.1 B[e]     v5: B[e] ->
%   C branch:  v6: A[c] -> C[c]     v7: C[c] -> C[e]         v8: C[e] ->
% At most 1e-4 enters (v1). v4 can carry 1e-4, but v5 = 0.1 * v4 can carry at most 1e-5 < epsilon, so v5 never reaches
% epsilon. A method that drops v5 but keeps v3 and v4 returns a network in which B[e] is produced and never consumed:
% a dead end, so the "consistent" network is not consistent, and the method run again on it removes v3 and v4.
% The consistent set at epsilon (closed: running the method again changes nothing) is {v1, v2, v6, v7, v8}.
%
% The script runs every method it finds on the path (SPECTRA spectraCC, FASTCC from the COBRA Toolbox, FasterCC) on
% the toy, then again on its own output, and prints for each: the reactions kept, the dead-end metabolites left in the
% output, and whether the second run keeps the same set (a consistent output must be a fixed point).
% It also prints FVA of the toy (largest flux of each reaction).
%
% Requirements: COBRA Toolbox (createModel, addReaction, removeRxns, fluxVariability, fastcc) with an LP solver set by
% changeCobraSolver; for SPECTRA its src/ folder (github.com/bisect-group/SPECTRA); for FasterCC the folder
% ../fastercc_cplex (CPLEX MATLAB API). A method that is not on the path is skipped and reported as "not run".
% The same toy without MATLAB: python3 check_toys_python.py (FVA and exact Poolman only).

clear
epsilon = 1e-4;
here = fileparts(mfilename('fullpath'));
if exist('fasterCCV10_SnapshotRemove', 'file') ~= 2 && isfolder(fullfile(here, '..', 'fastercc_cplex'))
    addpath(fullfile(here, '..', 'fastercc_cplex'));
end

model = toyModelSpectra();
fprintf('\nToy model: %d reactions, %d metabolites, v1 <= %g, epsilon = %g\n', numel(model.rxns), numel(model.mets), model.ub(1), epsilon);

%% FVA: the largest flux each reaction can carry (no objective)
try
    m0 = model; m0.c(:) = 0;
    [minF, maxF] = fluxVariability(m0, 0, 'max', m0.rxns);
    fprintf('\nFVA (min / max flux):\n');
    for j = 1:numel(model.rxns)
        fprintf('  %-3s  %10.3g  %10.3g   %s\n', model.rxns{j}, minF(j), maxF(j), ...
            ternary(max(abs([minF(j) maxF(j)])) >= epsilon, 'reaches epsilon', 'below epsilon'));
    end
catch e
    fprintf('\nFVA not run: %s\n', e.message);
end

%% the three methods, each on the toy and again on its own output
methods = {'SPECTRA', 'FASTCC', 'FasterCC'};
fprintf('\n%-9s | %-34s | %-22s | %-34s | %s\n', 'method', 'first run keeps', 'dead ends in output', 'run on its own output keeps', 'fixed point');
fprintf('%s\n', repmat('-', 1, 125));
for k = 1:numel(methods)
    [ids1, err] = runMethod(methods{k}, model, epsilon);
    if ~isempty(err), fprintf('%-9s | not run: %s\n', methods{k}, err); continue, end
    sub = removeRxns(model, setdiff(model.rxns, ids1));
    de = deadEndMets(sub);
    [ids2, err2] = runMethod(methods{k}, sub, epsilon);
    if ~isempty(err2), ids2 = {}; end
    fp = isempty(setxor(ids1, ids2));
    fprintf('%-9s | %-34s | %-22s | %-34s | %s\n', methods{k}, strjoin(sort(ids1), ','), strjoin(de, ','), ...
        strjoin(sort(ids2), ','), ternary(fp, 'yes', 'NO'));
end
fprintf('\nConsistent set at epsilon (closed): v1,v2,v6,v7,v8.  Every reaction can carry some flux (FVA), but v5 only 1e-5.\n');

%% ------------------------------------------------------------------------------------------------------------------
function model = toyModelSpectra()
model = createModel();
model = addReaction(model, 'v1', 'reactionFormula', ' -> A[e]', 'geneRule', 'Gene1');
model = addReaction(model, 'v2', 'reactionFormula', 'A[e] -> A[c]', 'geneRule', 'Gene2');
model = addReaction(model, 'v3', 'reactionFormula', 'A[c] -> B[c]', 'geneRule', 'Gene4');
model = addReaction(model, 'v4', 'reactionFormula', 'B[c] -> 0.1 B[e]', 'geneRule', 'Gene4');
model = addReaction(model, 'v5', 'reactionFormula', 'B[e] ->', 'geneRule', 'Gene4');
model = addReaction(model, 'v6', 'reactionFormula', 'A[c] -> C[c]', 'geneRule', 'Gene4');
model = addReaction(model, 'v7', 'reactionFormula', 'C[c] -> C[e]', 'geneRule', 'Gene4');
model = addReaction(model, 'v8', 'reactionFormula', 'C[e] ->', 'geneRule', 'Gene4');
model.ub(1) = 1e-4;
end

function [ids, err] = runMethod(name, model, epsilon)
% reaction IDs each method calls consistent; err = reason when the method is not available or fails
ids = {}; err = '';
try
    switch name
        case 'SPECTRA'
            if exist('spectraCC', 'file') ~= 2, err = 'spectraCC not on the path'; return, end
            rng(1, 'twister');                                   % spectraCC draws random LP weights
            A = spectraCC(model, epsilon, 'stoichiometry');
            ids = model.rxns(A(:));
        case 'FASTCC'
            if exist('fastcc', 'file') ~= 2, err = 'fastcc (COBRA) not on the path'; return, end
            A = fastcc(model, epsilon, 0);
            ids = model.rxns(A(:));
        case 'FasterCC'
            if exist('fasterCCV10_SnapshotRemove', 'file') ~= 2, err = 'fastercc_cplex not on the path'; return, end
            A = fasterCCV10_SnapshotRemove(model, epsilon, 0);
            ids = model.rxns(A(:));
    end
    ids = unique(ids(:))';
catch e
    err = e.message;
end
end

function de = deadEndMets(model)
% metabolites that can only be produced or only be consumed, given the reaction directions of the model
de = {};
for i = 1:numel(model.mets)
    j = find(model.S(i, :));
    if isempty(j), continue, end
    rev = model.lb(j) < 0 & model.ub(j) > 0;
    if any(rev), continue, end
    sgn = sign(full(model.S(i, j))) .* sign(model.ub(j)' + model.lb(j)');   % direction of use
    if all(sgn > 0) || all(sgn < 0), de{end + 1} = model.mets{i}; end %#ok<AGROW>
end
if isempty(de), de = {'none'}; end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
