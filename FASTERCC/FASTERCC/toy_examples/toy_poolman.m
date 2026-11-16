% toy_poolman.m -- the Poolman toy: what exact Poolman removes before any LP, and what it must not remove.
%
% Network: 5 metabolites (A-E), 7 irreversible reactions, bounds 0..10
%   R1:  -> A      R2: A -> B      R3: B ->         (the pathway: carries flux)
%   R4:  A -> C    R5: C -> D      R6: D -> C       (R4 feeds the pool C + D; nothing ever leaves it)
%   R7:  B -> E                                     (E is a dead end: produced, never consumed)
% Truth (FVA, exact):
%   R4 is blocked by stoichiometry alone: adding the balances of C and D gives d(C+D)/dt = v4, so at steady state v4 = 0.
%      Exact Poolman finds it (its row is zero in every null-space basis; certificate y = C + D, S'y = e_R4).
%   R7 is blocked: E has no consumer (a type I dead end), also a zero row.
%   R5 and R6 form a free loop (v5 = v6 balances C and D), so they CAN carry flux and are consistent for FASTCC and FVA.
% FasterCC's exact Poolman removes R4 and R7 before any LP and keeps the loop.
%
% The script prints the exact Poolman result (fastercc_cplex/poolmanExactJava.m, Java engine), FVA, and the sets kept by
% FASTCC, FasterCC and SPECTRA when they are on the path (each also run again on its own output).

clear
epsilon = 1e-4;
here = fileparts(mfilename('fullpath'));
if exist('fasterCCV10_SnapshotRemove', 'file') ~= 2 && isfolder(fullfile(here, '..', 'fastercc_cplex'))
    addpath(fullfile(here, '..', 'fastercc_cplex'));
end

model = createModel();
model = addReaction(model, 'R1', 'reactionFormula', ' -> A[c]');
model = addReaction(model, 'R2', 'reactionFormula', 'A[c] -> B[c]');
model = addReaction(model, 'R3', 'reactionFormula', 'B[c] ->');
model = addReaction(model, 'R4', 'reactionFormula', 'A[c] -> C[c]');
model = addReaction(model, 'R5', 'reactionFormula', 'C[c] -> D[c]');
model = addReaction(model, 'R6', 'reactionFormula', 'D[c] -> C[c]');
model = addReaction(model, 'R7', 'reactionFormula', 'B[c] -> E[c]');
model.lb(:) = 0; model.ub(:) = 10;
fprintf('\nPoolman toy: %d reactions, %d metabolites, bounds 0..10\n', numel(model.rxns), numel(model.mets));
disp(array2table(full(model.S), 'RowNames', model.mets, 'VariableNames', model.rxns));

%% exact Poolman (no LP): reactions whose row is zero in every null-space basis of S
if exist('poolmanExactJava', 'file') == 2
    [blk, info] = poolmanExactJava(model.S, struct('verify', true));
    fprintf('Exact Poolman blocks: %s   (rank S = %d, proofs verified: %d)\n', strjoin(model.rxns(blk), ', '), info.rankS, info.verified);
else
    fprintf('Exact Poolman not run: fastercc_cplex is not on the path\n');
end

%% FVA
try
    m0 = model; m0.c(:) = 0;
    [minF, maxF] = fluxVariability(m0, 0, 'max', m0.rxns);
    fprintf('\nFVA (min / max flux):\n');
    for j = 1:numel(model.rxns), fprintf('  %-3s  %6.3g  %6.3g\n', model.rxns{j}, minF(j), maxF(j)); end
catch e
    fprintf('\nFVA not run: %s\n', e.message);
end

%% the methods (same helper as toy_spectra_fastcc.m)
fprintf('\n%-9s | %-26s | %-26s | %s\n', 'method', 'first run keeps', 'run on its own output keeps', 'fixed point');
for name = {'FASTCC', 'FasterCC', 'SPECTRA'}
    [ids1, err] = run1(name{1}, model, epsilon);
    if ~isempty(err), fprintf('%-9s | not run: %s\n', name{1}, err); continue, end
    sub = removeRxns(model, setdiff(model.rxns, ids1));
    ids2 = run1(name{1}, sub, epsilon);
    fprintf('%-9s | %-26s | %-26s | %d\n', name{1}, strjoin(sort(ids1), ','), strjoin(sort(ids2), ','), isempty(setxor(ids1, ids2)));
end
fprintf('\nExpected (FVA): R1, R2, R3, R5, R6 consistent; R4 and R7 blocked.\n');

function [ids, err] = run1(name, model, epsilon)
ids = {}; err = '';
try
    switch name
        case 'SPECTRA'
            if exist('spectraCC', 'file') ~= 2, err = 'spectraCC not on the path'; return, end
            rng(1, 'twister'); A = spectraCC(model, epsilon, 'stoichiometry');
        case 'FASTCC'
            if exist('fastcc', 'file') ~= 2, err = 'fastcc (COBRA) not on the path'; return, end
            A = fastcc(model, epsilon, 0);
        case 'FasterCC'
            if exist('fasterCCV10_SnapshotRemove', 'file') ~= 2, err = 'fastercc_cplex not on the path'; return, end
            A = fasterCCV10_SnapshotRemove(model, epsilon, 0);
    end
    ids = unique(model.rxns(A(:)))';
catch e
    err = e.message;
end
end
