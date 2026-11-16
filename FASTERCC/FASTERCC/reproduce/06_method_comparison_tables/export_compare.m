function export_compare(what, outDir, limit, folders)
% export_compare  READ ONLY on all inputs. Exports, as plain TSV, what FasterCC (NEW), FASTCC and SPECTRA kept, with the
% input bounds of every reaction, so the comparison (counts, reversible/irreversible before and after, reaction IDs
% that differ) can be done outside MATLAB. Writes only into outDir; refuses to overwrite its own files.
%
%   export_compare('real',  OUT)        58 real models in bench2026/results_real/{NEW,FASTCC,SPECTRA}
%       OUT/real_rxns.tsv     one row per INPUT reaction per model:
%                             model rxn lb_in ub_in kept_NEW kept_FASTCC kept_SPECTRA lb_NEWout ub_NEWout
%                             (kept_* = 1/0, or NA when that method has no finished result; *_NEWout = bounds in
%                             FasterCC's returned resModel, NA if not kept)
%       OUT/real_runs.tsv     one row per model and method: n, finished, time, jobid, date, version (which), FVA check,
%                             normalization, input size, and ID checks (kept IDs not found in the input)
%   export_compare('bench', OUT)        1,140 CC tasks in results_full_cplex, results_full_cplex_spectra,
%                                       results_full_gurobi_spectra (folder CC)
%       OUT/bench_tasks.tsv   one row per task, folder and method: n, counts of the kept set by the input bounds
%                             (rev / fwd / bwd / closed), nRev reported by the method, time, error text, and the input
%                             counts (rev / fwd / bwd / closed; the input = base cluster model with removedIdx closed)
%                             plus a check that the computed input rev count equals res.nRev
%       OUT/bench_diffs.tsv   every reaction kept by one method and not by another, within the same result file
%   third argument (optional) = number of models / tasks per folder to export (dry run), e.g. 2
%   fourth argument (optional, bench only) = other result folders, e.g. {'results_spec15_cplex_spectra','results_spec15_gurobi_spectra'}
if nargin < 3 || isempty(limit), limit = Inf; end
if nargin < 4, folders = {}; end
B = '/mnt/scratch/users/egonzalez/fastbox2copy/bench2026';
if ~isfolder(outDir), mkdir(outDir); end
switch what
    case 'real',  doReal(B, outDir, limit);
    case 'bench', doBench(B, outDir, limit, folders);
    otherwise, error('first argument: real or bench');
end
end

% ---------------------------------------------------------------------------------------------------------------
function doReal(B, outDir, limit)
f1 = fullfile(outDir, 'real_rxns.tsv'); f2 = fullfile(outDir, 'real_runs.tsv');
if isfile(f1) || isfile(f2), error('%s or %s exists: not overwriting', f1, f2); end
M = {'NEW', 'FASTCC', 'SPECTRA'};
d = dir(fullfile(B, 'results_real', 'NEW', '*.mat')); names = sort({d.name});
names = names(1:min(numel(names), limit));
h1 = fopen(f1, 'w'); h2 = fopen(f2, 'w');
fprintf(h1, 'model\trxn\tlb_in\tub_in\tkept_NEW\tkept_FASTCC\tkept_SPECTRA\tlb_NEWout\tub_NEWout\n');
fprintf(h2, 'model\tmethod\tfinished\tn\ttime_s\tjobid\tdate\tnFVA\tnFalseBlocked\tnFalseConsistent\tcorrect\tnormalization\tnRxns_input\tids_not_in_input\terr\twhich\n');
for k = 1:numel(names)
    m = erase(names{k}, '.mat'); fprintf('[%d/%d] %s\n', k, numel(names), m);
    R = struct(); inputFile = ''; inputVar = '';
    for j = 1:3
        f = fullfile(B, 'results_real', M{j}, names{k});
        if isfile(f)
            S = load(f); R.(M{j}) = S;
            if isempty(inputFile) && isfield(S.res, 'file'), inputFile = char(S.res.file); inputVar = char(S.res.variable); end
        end
    end
    if isempty(inputFile), fprintf('   no result file with an input path: skipped\n'); continue, end
    I = load(inputFile, inputVar); model = I.(inputVar);
    rx = string(model.rxns(:)); lb = model.lb(:); ub = model.ub(:); n = numel(rx);
    kept = nan(n, 3); lbo = nan(n, 1); ubo = nan(n, 1);
    for j = 1:3
        if ~isfield(R, M{j}), fprintf(h2, '%s\t%s\tNA\t\t\t\t\t\t\t\t\t\t%d\t\tno result file\t\n', m, M{j}, n); continue, end
        r = R.(M{j}).res; fin = isfield(r, 'finished') && r.finished;
        ids = strings(0, 1); if isfield(r, 'ids'), ids = string(r.ids(:)); end
        if fin, kept(:, j) = ismember(rx, ids); end
        nf = sum(~ismember(ids, rx));
        fprintf(h2, '%s\t%s\t%d\t%d\t%g\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\n', m, M{j}, fin, numel(ids), ...
            num(r, 'time'), str(r, 'jobid'), str(r, 'date'), numstr(r, 'nFVA'), numstr(r, 'nFalseBlocked'), ...
            numstr(r, 'nFalseConsistent'), numstr(r, 'correct'), oneline(str(r, 'normalization')), n, nf, oneline(str(r, 'err')), oneline(str(r, 'which')));
        if j == 1 && fin && isfield(R.NEW, 'resModel')
            rm = R.NEW.resModel; [tf, loc] = ismember(rx, string(rm.rxns(:)));
            lbo(tf) = rm.lb(loc(tf)); ubo(tf) = rm.ub(loc(tf));
        end
    end
    for i = 1:n
        fprintf(h1, '%s\t%s\t%.17g\t%.17g\t%s\t%s\t%s\t%s\t%s\n', m, rx(i), lb(i), ub(i), na(kept(i, 1)), na(kept(i, 2)), ...
            na(kept(i, 3)), na(lbo(i)), na(ubo(i)));
    end
end
fclose(h1); fclose(h2); fprintf('written %s and %s\n', f1, f2);
end

% ---------------------------------------------------------------------------------------------------------------
function doBench(B, outDir, limit, folders)
f1 = fullfile(outDir, 'bench_tasks.tsv'); f2 = fullfile(outDir, 'bench_diffs.tsv');
if isfile(f1) || isfile(f2), error('%s or %s exists: not overwriting', f1, f2); end
F = {'results_full_cplex', 'results_full_cplex_spectra', 'results_full_gurobi_spectra'};
if ~isempty(folders), F = folders; end
M = {'NEW', 'FASTCC', 'SPECTRA'};
h1 = fopen(f1, 'w'); h2 = fopen(f2, 'w');
fprintf(h1, ['folder\ttask\tcluster\tlevel\trep\tsolver\tmethod\tn\tkept_rev\tkept_fwd\tkept_bwd\tkept_closed\tnRev_reported\t' ...
             'time_s\terr\tin_total\tin_rev\tin_fwd\tin_bwd\tin_closed\tnRev_res\tin_rev_check\tjobid\tdate\n']);
fprintf(h2, 'folder\ttask\tmethodA\tmethodB\tside\trxn\n');
% list tasks per cluster so each base model is loaded once
T = cell(0, 3);
for a = 1:numel(F)
    d = dir(fullfile(B, F{a}, 'CC', 'CC_c*.mat')); nm = sort({d.name}); nm = nm(1:min(numel(nm), limit));
    for k = 1:numel(nm), T(end + 1, :) = {F{a}, nm{k}, sscanf(nm{k}, 'CC_c%d')}; end %#ok<AGROW>
end
for c = unique(cell2mat(T(:, 3)))'
    fprintf('cluster %d: loading base model\n', c);
    S = load(fullfile(B, 'models', sprintf('expanded_input_cluster%d.mat', c)), 'model_in');
    rx0 = string(S.model_in.rxns(:)); lb0 = S.model_in.lb(:); ub0 = S.model_in.ub(:); clear S
    for t = find(cell2mat(T(:, 3)) == c)'
        r = load(fullfile(B, T{t, 1}, 'CC', T{t, 2})); r = r.res;
        lb = lb0; ub = ub0; lb(r.removedIdx) = 0; ub(r.removedIdx) = 0;
        [cr, cf, cb, cz] = cls(lb, ub); nres = num(r, 'nRev');
        common = sprintf('%d\t%d\t%d\t%d\t%d\t%g\t%s\t%s\t%s', numel(lb), cr, cf, cb, cz, nres, okstr(cr == nres), str(r, 'jobid'), str(r, 'date'));
        present = {};
        for j = 1:3
            fld = ['in_' M{j}];
            if ~isfield(r, fld), continue, end
            in = logical(r.(fld)(:)); present{end + 1} = M{j}; %#ok<AGROW>
            if numel(in) ~= numel(lb), error('%s %s: in_%s has %d entries, input has %d', T{t, 1}, T{t, 2}, M{j}, numel(in), numel(lb)); end
            [kr, kf, kb, kz] = cls(lb(in), ub(in));
            fprintf(h1, '%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%g\t%s\t%s\n', T{t, 1}, r.task, r.cluster, r.level, r.rep, ...
                str(r, 'solver'), M{j}, nnz(in), kr, kf, kb, kz, numstr(r, ['nRev_' M{j}]), num(r, ['time_' M{j}]), ...
                oneline(str(r, ['err_' M{j}])), common);
        end
        for p = 1:numel(present)
            for q = p + 1:numel(present)
                A = logical(r.(['in_' present{p}])(:)); Bq = logical(r.(['in_' present{q}])(:));
                for i = find(A & ~Bq)', fprintf(h2, '%s\t%d\t%s\t%s\tonly_A\t%s\n', T{t, 1}, r.task, present{p}, present{q}, rx0(i)); end
                for i = find(~A & Bq)', fprintf(h2, '%s\t%d\t%s\t%s\tonly_B\t%s\n', T{t, 1}, r.task, present{p}, present{q}, rx0(i)); end
            end
        end
    end
end
fclose(h1); fclose(h2); fprintf('written %s and %s\n', f1, f2);
end

% ---------------------------------------------------------------------------------------------------------------
function [r, f, b, z] = cls(lb, ub)
z = nnz(lb == 0 & ub == 0); r = nnz(lb < 0 & ub > 0); f = nnz(lb >= 0 & ub > 0); b = nnz(lb < 0 & ub <= 0);
end
function s = str(r, k), if isfield(r, k), s = char(join(string(r.(k)), ' ')); else, s = ''; end, end
function x = num(r, k), if isfield(r, k) && isnumeric(r.(k)) && ~isempty(r.(k)), x = double(r.(k)(1)); else, x = NaN; end, end
function s = numstr(r, k), if isfield(r, k) && ~isempty(r.(k)), s = num2str(double(r.(k)(1))); else, s = 'NA'; end, end
function s = na(x), if isnan(x), s = 'NA'; else, s = sprintf('%.17g', x); end, end
function s = okstr(b), if b, s = 'OK'; else, s = 'DIFF'; end, end
function s = oneline(s), s = regexprep(s, '[\t\r\n]+', ' '); end
