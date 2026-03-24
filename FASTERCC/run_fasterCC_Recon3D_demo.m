%% run_fasterCC_Recon3D_demo.m
% =========================================================
% fasterCC demo on Recon3D
%
% Goals:
%   1) demonstrate how to run fasterCC on Recon3D
%   2) compare results with and without reversibility check (revopt)
%   3) inspect changes in the model (reactions and bounds)
%   4) verify correct restoration of bounds (restoreBounds)
%   5) use FVA as a sanity check of the final model
%
% Required functions in path:
%   Core functions:
%       fasterCC.m
%       fixIrr_rFASTCORMICS.m
%       restoreBounds.m
%       findDeadEndsFastbox.m
%       structureAnalyseFastbox.m
%
%   LP helper functions:
%       LP3_rFASTCORMICS.m
%       LP7_rFASTCORMICS.m
%       LP9_rFASTCORMICS.m
%
%   COBRA Toolbox functions:
%       removeRxns
%       fluxVariability
%       changeCobraSolver
%       initCobraToolbox   (if needed)
%
%   Model file:
%       Recon3DModel.mat
% =========================================================

clear; clc;

%% 0. User settings
modelFile  = 'Recon3D.mat';   % adjust if needed 
epsilon    = 1e-4;
printLevel = 1;
solverName = 'ibm_cplex';          % adjust if needed
tolFVA     = epsilon;

outBase = 'fastercc_vs_fva_results_Recon3D';
if ~exist(outBase,'dir')
    mkdir(outBase);
end
modelTag = 'Recon3D';
modelOutDir = fullfile(outBase, modelTag);
if ~exist(modelOutDir, 'dir')
    mkdir(modelOutDir);
end

fprintf('=========================================================\n');
fprintf('              fasterCC DEMO ON Recon3D\n');
fprintf('=========================================================\n\n');

%% 1. Check required functions
fprintf('1) Checking required functions in the MATLAB path...\n');

requiredCore = { ...
    'fasterCC', ...
    'fixIrr_rFASTCORMICS', ...
    'restoreBounds', ...
    'findDeadEndsFastbox', ...
    'structureAnalyseFastbox', ...
    'LP3_rFASTCORMICS', ...
    'LP7_rFASTCORMICS', ...
    'LP9_rFASTCORMICS'};

requiredCobra = { ...
    'removeRxns', ...
    'fluxVariability', ...
    'changeCobraSolver'};

missingFlag = false;

fprintf('   Custom fasterCC functions:\n');
for i = 1:numel(requiredCore)
    p = which(requiredCore{i});
    if isempty(p)
        fprintf('      [MISSING] %s\n', requiredCore{i});
        missingFlag = true;
    else
        fprintf('      [FOUND]   %s\n', requiredCore{i});
    end
end

fprintf('\n   COBRA functions:\n');
for i = 1:numel(requiredCobra)
    p = which(requiredCobra{i});
    if isempty(p)
        fprintf('      [MISSING] %s\n', requiredCobra{i});
        missingFlag = true;
    else
        fprintf('      [FOUND]   %s\n', requiredCobra{i});
    end
end

if missingFlag
    error('Some required functions are missing from the MATLAB path. Please fix the path and rerun.');
end
fprintf('\n');

%% 2. Initialize COBRA / solver
fprintf('2) Initializing COBRA / solver...\n');

try
    changeCobraSolver(solverName, 'LP', 0, -1);
    fprintf('   Solver set to: %s\n', solverName);
catch
    fprintf('   Warning: could not set solver automatically.\n');
    fprintf('   Please make sure a COBRA LP solver is already configured.\n');
end
fprintf('\n');

%% 3. Load Recon3D model
fprintf('3) Loading Recon3D model...\n');

tmp = load(modelFile);
fn = fieldnames(tmp);
model = [];
modelName = '';

for k = 1:numel(fn)
    cand = tmp.(fn{k});
    if isstruct(cand) && isfield(cand,'S') && isfield(cand,'rxns') && ...
            isfield(cand,'mets') && isfield(cand,'lb') && isfield(cand,'ub')
        model = cand;
        modelName = fn{k};
        break;
    end
end

if isempty(model)
    error('No COBRA model structure found in %s', modelFile);
end

model_org = model;

fprintf('   Model variable detected: %s\n', modelName);
fprintf('   Number of reactions   : %d\n', numel(model_org.rxns));
fprintf('   Number of metabolites : %d\n', numel(model_org.mets));

if ~isfield(model_org, 'rev') || numel(model_org.rev) ~= numel(model_org.rxns)
    model_org.rev = double(model_org.lb < 0 & model_org.ub > 0);
end
fprintf('   Reactions marked reversible: %d\n', nnz(model_org.rev == 1));
fprintf('\n');

%% 4. Run fasterCC with revopt = 0
fprintf('4) Running fasterCC with revopt = 0...\n');

tStart0 = tic;
[A0, model_cc0] = fasterCC(model_org, epsilon, printLevel, 0);
time_fasterCC_0 = toc(tStart0);

fprintf('   ===== Results for revopt = 0 =====\n');
fprintf('   Kept reactions   : %d\n', numel(model_cc0.rxns));
fprintf('   Removed reactions: %d\n', numel(model_org.rxns) - numel(model_cc0.rxns));
fprintf('   Runtime (s)      : %.2f\n', time_fasterCC_0);

if numel(model_cc0.rxns) == numel(model_org.rxns)
    fprintf('   Interpretation: Recon3D appears fully consistent in this run.\n');
else
    fprintf('   Interpretation: some reactions were removed as inconsistent.\n');
end
fprintf('\n');

%% 5. Run fasterCC with revopt = 1
fprintf('5) Running fasterCC with revopt = 1...\n');

tStart1 = tic;
[A1, model_cc1] = fasterCC(model_org, epsilon, printLevel, 1);
time_fasterCC_1 = toc(tStart1);

fprintf('   ===== Results for revopt = 1 =====\n');
fprintf('   Kept reactions   : %d\n', numel(model_cc1.rxns));
fprintf('   Removed reactions: %d\n', numel(model_org.rxns) - numel(model_cc1.rxns));
fprintf('   Runtime (s)      : %.2f\n', time_fasterCC_1);

if numel(model_cc1.rxns) == numel(model_org.rxns)
    fprintf('   Interpretation: Recon3D remains fully consistent after reversibility checking.\n');
else
    fprintf('   Interpretation: some reactions were removed during consistency analysis.\n');
end
fprintf('\n');

%% 6. Compare revopt = 0 and revopt = 1
fprintf('6) Comparing revopt = 0 and revopt = 1...\n');

only_in_0 = setdiff(model_cc0.rxns, model_cc1.rxns);
only_in_1 = setdiff(model_cc1.rxns, model_cc0.rxns);

fprintf('   Reactions only in revopt = 0 model: %d\n', numel(only_in_0));
fprintf('   Reactions only in revopt = 1 model: %d\n', numel(only_in_1));

[sharedRxns01, i0, i1] = intersect(model_cc0.rxns, model_cc1.rxns);
sameLB01 = sum(abs(model_cc0.lb(i0) - model_cc1.lb(i1)) < epsilon);
sameUB01 = sum(abs(model_cc0.ub(i0) - model_cc1.ub(i1)) < epsilon);

fprintf('   Shared reactions: %d\n', numel(sharedRxns01));
fprintf('   Same lower bounds on shared reactions: %d / %d\n', sameLB01, numel(sharedRxns01));
fprintf('   Same upper bounds on shared reactions: %d / %d\n', sameUB01, numel(sharedRxns01));

if isempty(only_in_0) && isempty(only_in_1)
    fprintf('   Interpretation: both runs keep the same reaction set.\n');
else
    fprintf('   Interpretation: reaction sets differ between the two runs.\n');
end

if sameLB01 == numel(sharedRxns01) && sameUB01 == numel(sharedRxns01)
    fprintf('   Interpretation: revopt does not change bounds on shared reactions.\n');
else
    fprintf('   Interpretation: revopt changes some bounds, which can be expected when reversibility is tightened.\n');
end
fprintf('\n');

rev0 = find(model_cc0.lb < -epsilon & model_cc0.ub > epsilon);
rev1 = find(model_cc1.lb < -epsilon & model_cc1.ub > epsilon);

fprintf('Reversible reactions in revopt = 0: %d\n', numel(rev0));
fprintf('Reversible reactions in revopt = 1: %d\n', numel(rev1));
fprintf('Reactions no longer reversible after revopt = 1: %d\n', numel(rev0) - numel(rev1));



%% 7. Verify restoreBounds correctness (quality check)
fprintf('7) Checking restoreBounds correctness...\n');

orig_lb = model_org.lb;
orig_ub = model_org.ub;
changed_lb = find(model_org.lb > 0);
changed_ub = find(model_org.ub < 0);

fprintf('   Reactions with temporarily changed lower bound: %d\n', numel(changed_lb));
fprintf('   Reactions with temporarily changed upper bound: %d\n', numel(changed_ub));

fprintf('\n   Validation output for revopt = 0:\n');
check_restoreBounds(model_cc0, A0, model_org, changed_lb, changed_ub, orig_lb, orig_ub);

fprintf('\n   Validation output for revopt = 1:\n');
check_restoreBounds(model_cc1, A1, model_org, changed_lb, changed_ub, orig_lb, orig_ub);

fprintf('\n   Interpretation of restoreBounds check:\n');
fprintf('   - Four 1''s and a final 0 mean the restoration is correct.\n');
fprintf('   - Any 0 in the first four checks means bounds were not restored correctly.\n');
fprintf('   - A final number > 0 means some modified reactions were not classified.\n\n');

%% 8. Inspect bound changes relative to the original model
fprintf('8) Detailed reporting of bound changes by step...\n');

nRxns = numel(model_org.rxns);

% INITIAL STATE

rev_init = find(model_org.lb < -epsilon & model_org.ub > epsilon);
nRev_init = numel(rev_init);

backward_init = find(model_org.ub <= 0 & model_org.lb < 0);
nBackward_init = numel(backward_init);

fprintf('\n--- INITIAL MODEL ---\n');
fprintf('Total reactions: %d\n', nRxns);
fprintf('Reversible reactions: %d\n', nRev_init);
fprintf('Backward irreversible reactions: %d\n', nBackward_init);

% STEP 1 — BACKWARD NORMALIZATION

[tfBack, locBack] = ismember(model_org.rxns(backward_init), model_cc1.rxns);

backInFinal = backward_init(tfBack);
locFinal = locBack(tfBack);

lbChangedBack = abs(model_org.lb(backInFinal) - model_cc1.lb(locFinal)) > epsilon;
ubChangedBack = abs(model_org.ub(backInFinal) - model_cc1.ub(locFinal)) > epsilon;

nBackwardNormalized = nnz(lbChangedBack | ubChangedBack);

fprintf('\n--- STEP 1: backward ? forward normalization ---\n');
fprintf('Backward irreversible reactions normalized: %d / %d\n', ...
    nBackwardNormalized, nBackward_init);

% STEP 2 — revopt = 0 effect

rev0 = find(model_cc0.lb < -epsilon & model_cc0.ub > epsilon);
nRev0 = numel(rev0);

rev0_ids = string(model_cc0.rxns(rev0));
rev_init_ids = string(model_org.rxns(rev_init));

lostRev0 = setdiff(rev_init_ids, rev0_ids);
nLostRev0 = numel(lostRev0);

fprintf('\n--- STEP 2: revopt = 0 ---\n');
fprintf('Reversible before: %d\n', nRev_init);
fprintf('Reversible after revopt=0: %d\n', nRev0);
fprintf('Reversible ? irreversible (revopt=0): %d\n', nLostRev0);

% STEP 3 — revopt = 1 effect

rev1 = find(model_cc1.lb < -epsilon & model_cc1.ub > epsilon);
nRev1 = numel(rev1);

rev1_ids = string(model_cc1.rxns(rev1));

lostRev1_extra = setdiff(rev0_ids, rev1_ids);
nLostRev1_extra = numel(lostRev1_extra);

fprintf('\n--- STEP 3: revopt = 1 (extra tightening) ---\n');
fprintf('Reversible after revopt=0: %d\n', nRev0);
fprintf('Reversible after revopt=1: %d\n', nRev1);
fprintf('Additional reversible ? irreversible: %d\n', nLostRev1_extra);

% FINAL — TOTAL BOUND CHANGES

[shared1, iOrg1, iCC1] = intersect(model_org.rxns, model_cc1.rxns);

lbChanged = abs(model_org.lb(iOrg1) - model_cc1.lb(iCC1)) > epsilon;
ubChanged = abs(model_org.ub(iOrg1) - model_cc1.ub(iCC1)) > epsilon;

nTotalChanged = nnz(lbChanged | ubChanged);

fprintf('\n--- FINAL vs ORIGINAL ---\n');
fprintf('Total reactions with changed bounds: %d\n', nTotalChanged);

% BREAKDOWN 

fprintf('\n--- BREAKDOWN OF CHANGES ---\n');
fprintf('Backward normalization        : %d\n', nBackwardNormalized);
fprintf('revopt=0 reversible loss      : %d\n', nLostRev0);
fprintf('revopt=1 additional tightening: %d\n', nLostRev1_extra);

fprintf('Note: total changes (%d) may include overlaps and flips.\n', nTotalChanged);
fprintf('=====================================================\n\n');

%% 9. FVA-based comparison with fasterCC on Recon3D
% fprintf('9) FVA-based comparison with fasterCC on Recon3D...\n');
% 
% nRxns = numel(model_org.rxns);
% nMets = numel(model_org.mets);
% 
% revBeforeIdx = find(model_org.lb < -epsilon & model_org.ub > epsilon);
% revBeforeIDs = model_org.rxns(revBeforeIdx);
% nRevBefore = numel(revBeforeIdx);
% 
% fprintf('   Input model metabolites : %d\n', nMets);
% fprintf('   Input model reactions   : %d\n', nRxns);
% fprintf('   Reversible before fasterCC: %d\n', nRevBefore);
% 
% fprintf('\n   Running FVA on the full Recon3D model...\n');
% tStartFVA = tic;
% [lbFVA_full, ubFVA_full] = fluxVariability(model_org, 0, 'max');
% time_FVA_org = toc(tStartFVA);
% fprintf('   FVA finished in %.2f s\n', time_FVA_org);
% 
% %% FVA classification on reactions reversible before fasterCC
% lbFVA_rev = lbFVA_full(revBeforeIdx);
% ubFVA_rev = ubFVA_full(revBeforeIdx);
% 
% revFVA_mask          = (lbFVA_rev < -epsilon & ubFVA_rev > epsilon);
% irreversibleFwd_mask = (lbFVA_rev >= -epsilon & ubFVA_rev > epsilon);
% irreversibleBwd_mask = (ubFVA_rev <=  epsilon & lbFVA_rev < -epsilon);
% blocked_mask         = (lbFVA_rev >= -epsilon & ubFVA_rev <= epsilon);
% 
% revFVA_IDs          = revBeforeIDs(revFVA_mask);
% irreversibleFwd_IDs = revBeforeIDs(irreversibleFwd_mask);
% irreversibleBwd_IDs = revBeforeIDs(irreversibleBwd_mask);
% blocked_IDs         = revBeforeIDs(blocked_mask);
% 
% nRevFVA   = sum(revFVA_mask);
% nFwdOnly  = sum(irreversibleFwd_mask);
% nBwdOnly  = sum(irreversibleBwd_mask);
% nBlocked  = sum(blocked_mask);
% 
% fprintf('\n   FVA classification of reactions that were reversible before fasterCC:\n');
% fprintf('      Reversible by FVA   : %d\n', nRevFVA);
% fprintf('      Forward-only by FVA : %d\n', nFwdOnly);
% fprintf('      Backward-only by FVA: %d\n', nBwdOnly);
% fprintf('      Blocked by FVA      : %d\n', nBlocked);
% 
% %% Get reversible reactions after fasterCC
% nRxnsClean = numel(model_cc1.rxns);
% revAfterIdx = find(model_cc1.lb < -epsilon & model_cc1.ub > epsilon);
% revAfterIDs = model_cc1.rxns(revAfterIdx);
% nRevAfter = numel(revAfterIdx);
% 
% fprintf('\n   After fasterCC (revopt = 1):\n');
% fprintf('      Clean model reactions      : %d\n', nRxnsClean);
% fprintf('      Reversible after fasterCC  : %d\n', nRevAfter);
% 
% %% Compare reversible IDs
% revAfterStr = string(revAfterIDs(:));
% revFVAStr   = string(revFVA_IDs(:));
% 
% sharedRev = intersect(revAfterStr, revFVAStr);
% onlyFasterCC = setdiff(revAfterStr, revFVAStr);
% onlyFVA = setdiff(revFVAStr, revAfterStr);
% 
% nSharedRev = numel(sharedRev);
% nOnlyFasterCC = numel(onlyFasterCC);
% nOnlyFVA = numel(onlyFVA);
% 
% fprintf('\n   Comparison of reversible reactions: fasterCC vs FVA\n');
% fprintf('      Shared reversible reactions : %d\n', nSharedRev);
% fprintf('      Only reversible in fasterCC : %d\n', nOnlyFasterCC);
% fprintf('      Only reversible in FVA      : %d\n', nOnlyFVA);
% 
% %% Build reaction-level comparison table
% modelRxns = string(model_org.rxns(:));
% 
% revBeforeLogical = false(nRxns,1);
% revBeforeLogical(revBeforeIdx) = true;
% 
% revFVALogical = false(nRxns,1);
% revFVALogical(revBeforeIdx(revFVA_mask)) = true;
% 
% revAfterLogical = ismember(modelRxns, string(revAfterIDs(:)));
% inCleanModelLogical = ismember(modelRxns, string(model_cc1.rxns(:)));
% 
% Tcompare = table();
% Tcompare.rxnID = modelRxns;
% Tcompare.lb = model_org.lb(:);
% Tcompare.ub = model_org.ub(:);
% Tcompare.revBefore = revBeforeLogical;
% Tcompare.FVA_min = lbFVA_full(:);
% Tcompare.FVA_max = ubFVA_full(:);
% Tcompare.revByFVA = revFVALogical;
% Tcompare.inCleanModel = inCleanModelLogical;
% Tcompare.revAfterFasterCC = revAfterLogical;
% 
% %% Build summary table
% summaryRecon3D = table( ...
%     string(modelTag), ...
%     nMets, ...
%     nRxns, ...
%     nRevBefore, ...
%     nRxnsClean, ...
%     nRevAfter, ...
%     nRevFVA, ...
%     nFwdOnly, ...
%     nBwdOnly, ...
%     nBlocked, ...
%     nSharedRev, ...
%     nOnlyFasterCC, ...
%     nOnlyFVA, ...
%     time_fasterCC_0, ...
%     time_fasterCC_1, ...
%     time_FVA_org, ...
%     'VariableNames', { ...
%         'Model', ...
%         'InputModelMets', ...
%         'InputModelSize', ...
%         'ReversibleBeforeFasterCC', ...
%         'CleanModelSize', ...
%         'ReversibleAfterFasterCC', ...
%         'ReversibleByFVA', ...
%         'FVA_FwdOnly', ...
%         'FVA_BwdOnly', ...
%         'FVA_Blocked', ...
%         'SharedReversible_FasterCC_FVA', ...
%         'OnlyFasterCCReversible', ...
%         'OnlyFVAReversible', ...
%         'fasterCC_revopt0_sec', ...
%         'fasterCC_revopt1_sec', ...
%         'FVA_fullmodel_sec' ...
%     });
% 
% fprintf('\n   Summary table:\n');
% disp(summaryRecon3D);
% 
% %% Save outputs
% save(fullfile(modelOutDir, 'Recon3D_fastercc_vs_fva.mat'), ...
%     'epsilon', 'model_org', 'model_fix', 'model_cc0', 'model_cc1', ...
%     'A0', 'A1', 'Tcompare', 'summaryRecon3D', ...
%     'revBeforeIDs', 'revFVA_IDs', 'revAfterIDs', ...
%     'irreversibleFwd_IDs', 'irreversibleBwd_IDs', 'blocked_IDs', ...
%     'sharedRev', 'onlyFasterCC', 'onlyFVA', ...
%     'lbFVA_full', 'ubFVA_full', ...
%     'time_fasterCC_0', 'time_fasterCC_1', 'time_FVA_org', '-v7.3');
% 
% writetable(Tcompare, fullfile(modelOutDir, 'Recon3D_reaction_comparison.csv'));
% writetable(summaryRecon3D, fullfile(modelOutDir, 'Recon3D_summary.csv'));
% 
% fprintf('\n   Saved outputs in: %s\n\n', modelOutDir);
% 
% %% 11. Compact summary for presentation
% fprintf('==================== SUMMARY ====================\n');
% fprintf('Original model reactions : %d\n', numel(model_org.rxns));
% fprintf('revopt = 0 kept          : %d\n', numel(model_cc0.rxns));
% fprintf('revopt = 1 kept          : %d\n', numel(model_cc1.rxns));
% fprintf('Only in revopt = 0       : %d\n', numel(only_in_0));
% fprintf('Only in revopt = 1       : %d\n', numel(only_in_1));
% fprintf('Changed bounds revopt=1  : %d\n', numel(idxShow));
% fprintf('fasterCC time revopt=0 s : %.2f\n', time_fasterCC_0);
% fprintf('fasterCC time revopt=1 s : %.2f\n', time_fasterCC_1);
% fprintf('FVA full model time s    : %.2f\n', time_FVA_org);
% fprintf('Shared reversible        : %d\n', nSharedRev);
% fprintf('Only fasterCC reversible : %d\n', nOnlyFasterCC);
% fprintf('Only FVA reversible      : %d\n', nOnlyFVA);
% fprintf('=================================================\n');

%% =========================================================
% Local function: validate restoreBounds
%% =========================================================
function check_restoreBounds(model_cons, A, model_org, changed_lb, changed_ub, orig_lb, orig_ub)

changed_all = union(changed_lb, changed_ub);

[tf_met, loc_met] = ismember(model_cons.mets, model_org.mets);
assert(all(tf_met), 'Metabolite mismatch between final and original model.');

S_cons = model_cons.S;
S_orig = model_org.S(loc_met, A);

sameDir = find(all(S_orig == S_cons, 1));
oppDir  = find(all(S_orig == -S_cons, 1));

changed_idx_model = find(ismember(A, changed_all));

idx_same = intersect(sameDir, changed_idx_model);
idx_opp  = intersect(oppDir, changed_idx_model);

ok_same_lb = isequal(model_cons.lb(idx_same), orig_lb(A(idx_same)));
ok_same_ub = isequal(model_cons.ub(idx_same), orig_ub(A(idx_same)));
ok_opp_lb  = isequal(model_cons.lb(idx_opp), -orig_ub(A(idx_opp)));
ok_opp_ub  = isequal(model_cons.ub(idx_opp), -orig_lb(A(idx_opp)));

idx_other = setdiff(changed_idx_model, union(idx_same, idx_opp));

fprintf('      Same-dir LB ok: %d\n', ok_same_lb);
fprintf('      Same-dir UB ok: %d\n', ok_same_ub);
fprintf('      Opp-dir  LB ok: %d\n', ok_opp_lb);
fprintf('      Opp-dir  UB ok: %d\n', ok_opp_ub);
fprintf('      Unclassified changed reactions: %d\n', numel(idx_other));
end