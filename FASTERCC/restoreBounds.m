function model = restaureBounds(model, A, modelOrg, changed_lb, changed_ub, orig_lb, orig_ub)
% RESTAUREBOUNDS Restore reaction bounds from original model,
% accounting for possible reaction direction flips.
%
% INPUTS:
% model        - constrained model
% A            - mapping: model index -> original model index
% modelOrg     - original model
% changed_lb   - indices (in original model) where lb was changed
% changed_ub   - indices (in original model) where ub was changed
% orig_lb      - original lower bounds
% orig_ub      - original upper bounds

%% 1. Match reactions and metabolites
[~, IA, IB]   = intersect(modelOrg.rxns, model.rxns);
[~, IAA, IBB] = intersect(modelOrg.mets, model.mets);

%% 2. Extract aligned stoichiometric matrices
S_orig = modelOrg.S(IAA, IA);
S_cons = model.S(IBB, IB);

%% 3. Detect reaction direction consistency
SameDirr = find(all(S_orig == S_cons, 1));     % identical columns
OppDirr  = find(all(S_orig == -S_cons, 1));    % exact sign flip

% (Optional safety: reactions that are neither identical nor flipped)
% OtherDirr = setdiff(1:length(IA), union(SameDirr, OppDirr));

%% 4. Identify which reactions were modified (map via A)
changed_all = union(changed_lb, changed_ub);

% Find indices in *model* corresponding to changed reactions in original
changed_idx_model = find(ismember(A, changed_all));

%% 5. Restore bounds for SAME direction reactions
idx_same = intersect(SameDirr, changed_idx_model);

model.lb(idx_same) = orig_lb(A(idx_same));
model.ub(idx_same) = orig_ub(A(idx_same));

%% 6. Restore bounds for OPPOSITE direction reactions (flip bounds)
idx_opp = intersect(OppDirr, changed_idx_model);

model.lb(idx_opp) = -orig_ub(A(idx_opp));
model.ub(idx_opp) = -orig_lb(A(idx_opp));

end