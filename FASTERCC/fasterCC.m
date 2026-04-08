function [A, model] = fasterCC(model, epsilon, printLevel, revopt)
% [A, V] = fasterCC(model, epsilon, printLevel, revopt)
%
% The fasterCC algorithm for testing the consistency of a stoichiometric model.
% Output A is the consistent part of the model.
%
% INPUT
% model         COBRA model structure containing the fields:
%   S           m x n stoichiometric matrix
%   lb          n x 1 flux lower bound
%   ub          n x 1 flux upper bound
%   rxns        n x 1 cell array of reaction abbreviations
%
% epsilon       numerical tolerance threshold
% printLevel    0 = silent, 1 = summary
% revopt        0 = reversibility check is not run, 1 = reversibility check is run
%
% OUTPUT
% A             n x 1 boolean vector indicating the consistent flux
% model         COBRA model structure containing the consistent reactions
%
% (c) Maria Pacheco, Evelyn Gonzalez, Thomas Sauter, 2026

model_org = model;

if nargin < 4
    revopt = 0;
end
if nargin < 3
    printLevel = 0;
end
if nargin < 2
    epsilon = 1e-4;
end
%%

%% Save original bounds
orig_lb = model.lb;
orig_ub = model.ub;


%% Track constrained reactions
changed_lb = find(model.lb > 0);
changed_ub = find(model.ub < 0);

%% Relax bounds
model.lb(changed_lb) = 0;
model.ub(changed_ub) = 0;

%{
Part 1: Structural preprocessing. Prior to any flux testing, type I and type II
dead-ends are removed to reduce model size and avoid the costly one-by-one step.
The reversibility of reactions is then corrected to maximise the number of
reactions included in the first batch LP7 call. Note that findDeadEndsFastbox
internally calls fixIrr_rFASTCORMICS and structureAnalyseFastbox to maximise
dead-end identification.
%}

% Remove type I and II dead ends
[model] = findDeadEndsFastbox(model);
% Take advantage of the structure to turn some reversible reactions into
% irreversible reactions
[model] = structureAnalyseFastbox(model);
% Fix irreversible reactions to have (lb >= 0 and ub > 0)
model = fixIrr_rFASTCORMICS(model);

%% Run the consistency checking based on fastcc (Vlassis et al. 2014)
%{
Part 2: Flux testing of irreversible reactions. LP7 (Vlassis et al. 2014) is
applied to the full set of irreversible reactions (J) in a single batch.
Reactions that fail to carry flux are flagged as blocked (incI) and removed
from J. Those that carry flux, whether irreversible or reversible, are considered
flux-consistent and added to the active set (A). Reversible reactions carrying
a negative flux are flipped.
%}

% Number of reactions
N = 1:numel(model.rxns); % indices of the working model

% Reactions assumed to be irreversible in forward direction
I = find(model.lb >= 0 & model.ub > 0);
not_expressed = find(model.lb == 0 & model.ub == 0); % completely shut reactions

% Start with I
% J is the set of irreversible reactions that were not tested yet
J = intersect(N, I);
J = setdiff(J, not_expressed);

if printLevel == 1
    fprintf('|J| = %d ', numel(J));
end

% V is the flux vector that approximately maximizes the cardinality
V = LP7_rFASTCORMICS(J, model, epsilon);

% Flip reversible reactions with negative flux
Vrev = V < -epsilon & model.lb < 0 & model.ub > 0;
% Flip columns of S for those reactions and swap/negate their bounds
model.S(:, Vrev) = -model.S(:, Vrev);
temp = model.ub(Vrev);
model.ub(Vrev) = -model.lb(Vrev);
model.lb(Vrev) = -temp;

Supp = find(abs(V) >= 0.99 * epsilon);
% A is the set of reactions in V with absolute value greater than epsilon
A = Supp;
if printLevel == 1
    fprintf('|A| = %d\n', numel(A));
end

% At this step, all irreversible reactions should be in A. Those not in A are blocked.
incI = setdiff(J, A);

% Create a temporary model tmp in which dead-end metabolites (connected to
% only 1 reaction) and the associated reactions are removed
tmp = model;
[tmp] = removeRxns(tmp, tmp.rxns(incI));
%%
%{
Part 3: Iterative pruning. Blocked irreversible reactions (incI) are removed and a
temporary model is built. Dead-end metabolites (degree-1 nodes) and their associated
reactions are then pruned iteratively until no further dead-ends remain. Since the
temporary model is a subset of the original, the active set (A), the candidate set (J),
and the blocked set (incI) are remapped to the indices of the temporary model before
proceeding.
%}

r = numel(model.rxns);
while ~isempty(r)
    Struct_matrix = tmp.S ~= 0;
    dead_end_mets_I = sum(Struct_matrix, 2) == 1; % finds metabolites that participate in exactly one reaction in tmp (degree-1 nodes)
    deadendmets = tmp.mets(dead_end_mets_I);
    [~, r] = find(tmp.S(ismember(tmp.mets, deadendmets), :));
    tmp = removeRxns(tmp, tmp.rxns(r));
end

% Blocked reactions
block = setdiff(model.rxns, tmp.rxns);

% Find the index of the blocked reactions in the working model
blockID = find(ismember(model.rxns, block));
tmp = removeRxns(tmp, block);

if ~isempty(incI) && printLevel == 1
    fprintf('\n(inconsistent subset of I detected)\n');
end

% At this step, all irreversible reactions should be in A. Those not in A are blocked.
% Next step is to check reversible reactions.
J = setdiff(setdiff(N, A), incI);
J = setdiff(J, not_expressed);
J = setdiff(J, blockID);
J = find(ismember(tmp.rxns, model.rxns(J)));
A = find(ismember(tmp.rxns, model.rxns(A)));

model = tmp;
clear N V

if printLevel == 1
    fprintf('|J| = %d ', numel(J));
end
%%
%{
Part 4: As reversible exchange reactions are often part of the reactions tested in the
costly one-by-one step, to speed up the flux testing we first test exchange reactions
in the forward direction. Reactions carrying flux are added to A. Reversible reactions
with a negative flux are flipped. The set of reversible reactions that did not carry
flux are flipped and retested.
%}

% Testing exchange reactions
Ex = find(sum(model.S ~= 0, 1) == 1); % identify exchange-like reactions
V = LP7_rFASTCORMICS(Ex, model, epsilon); % test the consistency of all exchange reactions
Supp = find(abs(V) >= 0.99 * epsilon); % reactions that carry flux with an absolute positive value are not blocked
% Add supported exchanges to A, remove from J
% A is the set of reactions in V with absolute value greater than epsilon
A = union(A, Supp);
J = setdiff(J, A);

% If a reversible exchange came out negative, flip its orientation
% Flip reversible reactions that are negative
Vrev = V < 0 & model.lb < 0 & model.ub > 0;
model.S(:, Vrev) = -model.S(:, Vrev);
temp = model.ub(Vrev);
model.ub(Vrev) = -model.lb(Vrev);
model.lb(Vrev) = -temp;

% Flip exchange reactions that did not pass the consistency test
Ex = intersect(Ex, J);
model.S(:, Ex) = -model.S(:, Ex);
temp = model.ub(Ex);
model.ub(Ex) = -model.lb(Ex);
model.lb(Ex) = -temp;

% Recompute exchanges (still degree-1 after sign flip) and run LP7 again
V = LP7_rFASTCORMICS(Ex, model, epsilon);
Supp = find(abs(V) >= 0.99 * epsilon);
A = union(A, Supp);
J = setdiff(J, A);

flipped = false;
singleton = false;
%%
% We start in batch mode.
%{
Part 5: Main consistency loop.
The remaining reversible reactions are first tested in batch mode using LP7. Reactions
that fail to carry flux are reoriented and retested. If no progress is made, fasterCC
switches to a one-by-one test (LP3, see Vlassis et al. 2014). The key difference from
the original FASTCC implementation is that whenever a blocked reaction is identified,
findDeadEndsFastbox is called to detect and remove any newly arising dead-ends, keeping
the model as small as possible throughout the loop.
%}

while ~isempty(J)

    if singleton
        % Start of the one-by-one step
        Ji = J(1);
        V = LP3_rFASTCORMICS(Ji, model);

        % Flip reversible reactions with negative flux
        Vrev = V < 0 & model.rev == 1;
        % Optimal flux is negative and reactions are marked reversible; those can be re-oriented.
        model.S(:, Vrev) = -model.S(:, Vrev);
        temp = model.ub(Vrev);
        model.ub(Vrev) = -model.lb(Vrev);
        model.lb(Vrev) = -temp;

    else
        % Batch mode: try to support as many reactions as possible in one go.
        Ji = J;
        V = LP7_rFASTCORMICS(Ji, model, epsilon);

        Vrev = V < 0 & model.lb < 0 & model.ub > 0;
        % Flip reversible reactions that have a negative flux
        model.S(:, Vrev) = -model.S(:, Vrev);
        temp = model.ub(Vrev);
        model.ub(Vrev) = -model.lb(Vrev);
        model.lb(Vrev) = -temp;

    end

    Supp = find(abs(V) >= 0.99 * epsilon);
    A = union(A, Supp);
    if printLevel == 1
        fprintf('|A| = %d\n', numel(A));
    end

    if ~isempty(intersect(J, A))
        J = setdiff(J, A);
        if printLevel == 1
            fprintf('|J| = %d ', numel(J));
        end
        flipped = false;
    else
        JiRev = setdiff(Ji, I);
        if flipped || isempty(JiRev)
            flipped = false;
            if singleton
                J = setdiff(J, Ji);
                % Create a temporary model tmp in which dead-end metabolites
                % (connected to only 1 reaction) and the associated reactions are removed
                tmp = model;
                tmp = removeRxns(tmp, tmp.rxns(Ji));
                % Iteratively prune Class-I dead ends
                r = numel(model.rxns);
                while ~isempty(r)
                    Struct_matrix = tmp.S ~= 0;
                    dead_end_mets_I = sum(Struct_matrix, 2) == 1;
                    deadendmets = tmp.mets(dead_end_mets_I);
                    [~, r] = find(tmp.S(ismember(tmp.mets, deadendmets), :));
                    tmp = removeRxns(tmp, tmp.rxns(r));
                end

                % Blocked reactions
                block = setdiff(model.rxns, tmp.rxns);
                % Find the index of the blocked reactions in the original model
                blockID = find(ismember(model.rxns, block));
                J = setdiff(J, blockID);
                blockID = find(ismember(model.rxns, block));
                J = setdiff(J, blockID);
                J = find(ismember(tmp.rxns, model.rxns(J)));
                A = find(ismember(tmp.rxns, model.rxns(A)));
                model = tmp; % update the current model
                if printLevel == 1
                    fprintf('\n(inconsistent reversible reaction detected)\n');
                end
            else
                singleton = true;
            end
        else
            model.S(:, JiRev) = -model.S(:, JiRev);
            tmp = model.ub(JiRev);
            model.ub(JiRev) = -model.lb(JiRev);
            model.lb(JiRev) = -tmp;
            if singleton
                model.lb(JiRev) = 0;
                [tmp] = findDeadEndsFastbox(model);
                [tmp] = structureAnalyseFastbox(tmp);
                block = setdiff(model.rxns, tmp.rxns);
                blockID = find(ismember(model.rxns, block));
                J = setdiff(J, blockID);
                J = find(ismember(tmp.rxns, model.rxns(J)));
                A = find(ismember(tmp.rxns, model.rxns(A)));
                model = tmp;
                clear N V
            end
            flipped = true;
            if printLevel == 1
                fprintf('(flip) ');
            end
        end
    end
end

if numel(A) == numel(model_org.rxns)
    fprintf('\nThe input model is consistent.\n');
end

model = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns), A)));
A = find(ismember(model_org.rxns, model.rxns));
if ~isempty(setdiff(model_org.rxns(A), model.rxns))
    error('fasterCC: reactions in A are missing from the trimmed model. Index mismatch after removeRxns.');
end
if ~isempty(setdiff(model.rxns, model_org.rxns(A)))
    error('fasterCC: trimmed model contains reactions not present in A. Index mismatch after removeRxns.');
end
model_keep = model;
if revopt == 1
%{
Part 6: Optional full reversibility test. The reversible reactions remaining in
the consistent set are tested in both the forward and backward directions. In
each direction, LP7 is first applied in batch mode. Reactions that fail to carry
flux are then tested individually using the one-by-one step. Reactions confirmed
as unidirectional are tightened accordingly by setting their lower or upper bound
to zero.
%}

    % Test the reversible reactions in the forward direction
    rev = find(model.lb < 0 & model.ub > 0);
    V = LP7_rFASTCORMICS(rev, model, epsilon);

    SuppRev = find(abs(V) >= 0.99 * epsilon);
    NewIrr_test = setdiff(rev, SuppRev);
    % Flip all reversible reactions to test backward direction
    IrrTestkeep = zeros(numel(model.rxns), 1);

    for i = 1:numel(NewIrr_test)
        model = model_keep;
        IrrTest = NewIrr_test(i);
        V = LP7_rFASTCORMICS(IrrTest, model, epsilon);
        SuppRev = find(abs(V) >= 0.99 * epsilon);
        if ~isempty(setdiff(IrrTest, SuppRev))
            model.ub(IrrTest) = 0;
            IrrTestkeep(IrrTest) = 1;
        end
    end

    model.S(:, rev) = -model.S(:, rev);
    temp = model.ub(rev);
    model.ub(rev) = -model.lb(rev);
    model.lb(rev) = -temp;

    V = LP7_rFASTCORMICS(rev, model, epsilon);
    SuppRev = find(abs(V) >= 0.99 * epsilon);
    NewIrr_testback = setdiff(rev, SuppRev);
    % Flip all reversible reactions to test backward direction
    IrrTestkeepback = zeros(numel(model.rxns), 1);
    model_org_rev = model;

    for i = 1:numel(NewIrr_testback)
        model = model_org_rev;
        IrrTest = NewIrr_testback(i);
        V = LP7_rFASTCORMICS(IrrTest, model, epsilon);
        SuppRev = find(abs(V) >= 0.99 * epsilon);
        if ~isempty(setdiff(IrrTest, SuppRev))
            model.ub(IrrTest) = 0;
            IrrTestkeepback(IrrTest) = 1;
        end
    end

    model = model_keep;
    model.ub(IrrTestkeep ~= 0) = 0;
    model.lb(IrrTestkeepback ~= 0) = 0;
end
[model]=restoreBounds(model, A,model_org, changed_lb, changed_ub, orig_lb, orig_ub);

end
