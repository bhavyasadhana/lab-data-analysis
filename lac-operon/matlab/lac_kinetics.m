%% lac operon induction kinetics from the shared course table
%
%  1. Reads every shared table in data/shared (xlsx or csv) and cleans it:
%     column-name aliases, decimal commas, spelling variants, empty and
%     duplicate rows.
%  2. Where raw readings exist (data/lab_notebook/*.csv) it recomputes the
%     beta-galactosidase activity  U = A420*1000 / (t * V * A600)  and
%     cross-checks the shared table (wrong labels, swapped rows,
%     calculation errors), correcting it from the raw data.
%  3. Flags A420 > 1.0 (outside linear range), "> 3" (lower bound only) and
%     A420 < 0.02 (near the detection limit).
%  4. Pools groups, computes glycerol/glucose and IPTG/no-IPTG ratios and
%     draws a multi-curve kinetics plot. Tables open in their own windows.
%
%  Base MATLAB, R2020b or newer. Run from the matlab/ folder: >> lac_kinetics

clear; close all force;

P.tReaction  = 10;                 % min, ONPG reaction
P.a600Dil    = 5;                  % 200 µl + 800 µl Z buffer
P.volGlucose = 0.020;              % ml suspension in enzyme assay
P.volGlycerol= 0.010;
P.a420Max    = 1.0;                % linear range limit (script)
P.a420Min    = 0.020;              % near detection limit

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', 'data');
outDir  = fullfile(here, '..', 'results', 'matlab');
if ~exist(outDir, 'dir'), mkdir(outDir); end

shared = readShared(fullfile(dataDir, 'shared'));
nb     = readNotebooks(fullfile(dataDir, 'lab_notebook'), P);
[data, issues] = reconcile(shared, nb, P);
data   = addBaseline(data);
data   = sortrows(data, {'group','strain','csource','inducer','time_min'});
ratios = keyRatios(data);
pooled = pooledKinetics(data);

writetable(data,   fullfile(outDir, 'clean_data.csv'));
writetable(issues, fullfile(outDir, 'data_issues.csv'));
writetable(ratios, fullfile(outDir, 'key_ratios.csv'));
writetable(pooled, fullfile(outDir, 'pooled_kinetics.csv'));
plotKinetics(data, outDir);

fprintf('Read %d rows from data/shared, groups: %s\n', height(shared), mat2str(unique(shared.group)'));
disp('Problems found in the shared table'); disp(issues);
disp('Key ratios'); disp(ratios);
fprintf('Results written to %s\n', outDir);

showTable(issues, 'Problems found in the shared table');
showTable(data,   'Cleaned data (labels and activity from raw readings)');
showTable(ratios, 'Key ratios');
showTable(pooled, 'Pooled kinetics (reliable points only)');


%% ======================================================================
%  Reading and cleaning
%% ======================================================================
function T = readAsText(f)
    opts = detectImportOptions(f, 'VariableNamingRule', 'preserve');
    opts = setvartype(opts, 'string');
    T = readtable(f, opts);
end

function name = aliasColumn(n)
    n = lower(strtrim(n));
    map = containers.Map( ...
        {'gruppe','grp','nr','nr.','probe','sample','no.','stamm','c-source','c_source', ...
         'c-quelle','kohlenstoffquelle','medium','zeit','t','t (min)','time (min)','time_min', ...
         'induktor','activity','aktivität','aktivitaet','miller','units','beta-gal'}, ...
        {'group','group','no','no','no','no','no','strain','csource','csource', ...
         'csource','csource','csource','time','time','time','time','time', ...
         'inducer','betagal','betagal','betagal','betagal','betagal','betagal'});
    if isKey(map, n), name = map(n); else, name = matlab.lang.makeValidName(n); end
end

function [v, lower] = toNumber(s)
    % '1347,83', '1.347,83', '>3', '' -> numbers; lower = true for '>x'
    s = strrep(strtrim(string(s)), " ", "");
    lower = startsWith(s, ">");
    s = erase(s, [">", "<", "="]);
    v = nan(size(s));
    for i = 1:numel(s)
        x = s(i);
        if ismissing(x) || x == "", continue; end
        if contains(x, ",") && (~contains(x, ".") || max(strfind(x, ",")) > max(strfind(x, ".")))
            x = replace(replace(x, ".", ""), ",", ".");
        end
        v(i) = str2double(erase(x, ","));
    end
    lower(ismissing(lower)) = false;
end

function c = normCsource(s)
    s = lower(strtrim(string(s)));
    c = s;
    c(ismember(s, ["glu","glc","glucose","glukose"]))   = "glucose";
    c(ismember(s, ["gly","glyc","glycerol","glycerin"])) = "glycerol";
    c(ismember(s, ["lac","lactose","laktose"]))         = "lactose";
end

function c = normInducer(s)
    s = lower(strtrim(string(s)));
    s(ismissing(s)) = "";
    c = s;
    c(ismember(s, ["","-","–","ohne","none","no","kein","0","nan","w/o"])) = "none";
    c(contains(s, "iptg")) = "IPTG";
end

function c = normStrain(s)
    raw = strtrim(string(s));
    k = lower(regexprep(raw, '[\s_]', ''));
    c = raw;
    c(startsWith(k, "laci") & (contains(k, "neg") | endsWith(k, "-"))) = "lacI-";
    c(ismember(k, ["wt","wildtyp","wildtype","wildtypus"])) = "WT";
    c(startsWith(k, "crp") & (contains(k, "neg") | endsWith(k, "-"))) = "crp-";
end

function S = readShared(folder)
    files = [dir(fullfile(folder, '*.xlsx')); dir(fullfile(folder, '*.xls')); dir(fullfile(folder, '*.csv'))];
    if isempty(files), error('No shared tables found in %s', folder); end
    parts = {};
    for i = 1:numel(files)
        T = readAsText(fullfile(folder, files(i).name));
        T.Properties.VariableNames = cellfun(@aliasColumn, T.Properties.VariableNames, 'UniformOutput', false);
        A = T{:, :};                                          % all columns are strings
        T = T(~all(ismissing(A) | strtrim(A) == "", 2), :);   % drop empty rows
        R = table;
        R.group    = toNumber(T.group);
        R.no       = toNumber(T.no);
        R.strain   = normStrain(T.strain);
        R.csource  = normCsource(T.csource);
        R.inducer  = normInducer(T.inducer);
        R.time_min = toNumber(T.time);
        R.betagal_shared = toNumber(T.betagal);
        R.source_file = repmat(string(files(i).name), height(R), 1);
        parts{end+1} = R; %#ok<AGROW>
    end
    S = vertcat(parts{:});
    [~, ia] = unique(S(:, {'group','no','strain','csource','inducer','time_min'}), 'rows', 'stable');
    S = S(sort(ia), :);
end

function N = readNotebooks(folder, P)
    files = dir(fullfile(folder, '*.csv'));
    N = table;
    for i = 1:numel(files)
        T = readAsText(fullfile(folder, files(i).name));
        R = table;
        R.group    = toNumber(T.group);
        R.no       = toNumber(T.no);
        R.strain   = normStrain(T.strain);
        R.csource  = normCsource(T.csource);
        R.inducer  = normInducer(T.inducer);
        R.time_min = toNumber(T.time_min);
        R.a600     = toNumber(T.a600_diluted) * P.a600Dil;
        [a420, overflow] = toNumber(T.a420);
        R.a420 = a420; R.a420_overflow = overflow;
        N = [N; R]; %#ok<AGROW>
    end
    if isempty(N), return; end
    vol = nan(height(N), 1);
    vol(N.csource == "glucose")  = P.volGlucose;
    vol(N.csource == "glycerol") = P.volGlycerol;
    N.activity_U = N.a420 * 1000 ./ (P.tReaction .* vol .* N.a600);
end

function q = qualityFlag(a420, overflow, P)
    q = repmat("ok", numel(a420), 1);
    q(isnan(a420)) = "not verified (no raw data)";
    q(a420 < P.a420Min) = "A420 < 0.02 (near detection limit)";
    q(a420 > P.a420Max) = "A420 > 1.0 (outside linear range)";
    q(overflow) = "lower bound (A420 > 3)";
end

%% ======================================================================
%  Cross-check shared table against raw readings
%% ======================================================================
function [D, I] = reconcile(S, N, P)
    I = table('Size', [0 4], 'VariableTypes', {'double','double','string','string'}, ...
              'VariableNames', {'group','no','issue','detail'});
    lbl = @(T, i) sprintf('%s/%s/%g min', T.csource(i), T.inducer(i), T.time_min(i));
    same = @(A, i, B, j) A.strain(i) == B.strain(j) && A.csource(i) == B.csource(j) && ...
                         A.inducer(i) == B.inducer(j) && A.time_min(i) == B.time_min(j);
    for i = 1:height(N)
        j = find(S.group == N.group(i) & S.no == N.no(i), 1);
        if isempty(j)
            I(end+1, :) = {N.group(i), N.no(i), "missing from shared table", ""}; %#ok<AGROW>
            continue
        end
        if ~same(S, j, N, i)
            % partner = raw row whose real labels equal this shared row's labels
            swap = "";
            for k = find(N.group == N.group(i) & N.no ~= N.no(i))'
                if same(N, k, S, j)
                    jk = find(S.group == N.group(k) & S.no == N.no(k), 1);
                    if ~isempty(jk) && same(S, jk, N, i)
                        swap = sprintf('; swapped with no. %d', N.no(k));
                    end
                end
            end
            I(end+1, :) = {N.group(i), N.no(i), "condition label wrong in shared table", ...
                string(sprintf('shared says %s, raw data says %s%s', lbl(S, j), lbl(N, i), swap))}; %#ok<AGROW>
        end
        if isfinite(N.activity_U(i)) && abs(S.betagal_shared(j) - N.activity_U(i)) > 0.01 * abs(N.activity_U(i))
            I(end+1, :) = {N.group(i), N.no(i), "activity differs from recalculation", ...
                string(sprintf('shared %.2f U, recalculated %.2f U', S.betagal_shared(j), N.activity_U(i)))}; %#ok<AGROW>
        end
    end

    cols = {'group','no','strain','csource','inducer','time_min','a600','a420','activity_U'};
    D = N(:, cols);
    D.overflow = N.a420_overflow;
    D.source = repmat("raw readings", height(D), 1);
    % shared rows without raw data
    if isempty(N), keep = true(height(S), 1);
    else, keep = ~ismember([S.group S.no], [N.group N.no], 'rows'); end
    U = S(keep, {'group','no','strain','csource','inducer','time_min'});
    U.a600 = nan(height(U), 1); U.a420 = nan(height(U), 1);
    U.activity_U = S.betagal_shared(keep);
    U.overflow = false(height(U), 1);
    U.source = repmat("shared table only", height(U), 1);
    D = [D; U];
    D.quality = qualityFlag(D.a420, D.overflow, P);
    D.overflow = [];
end

%% ======================================================================
%  Kinetics summary
%% ======================================================================
function D = addBaseline(D)
    % IPTG is added right after t = 0, so the uninduced t0 of the same
    % group / strain / carbon source is also t0 of the IPTG curve.
    add = D([], :);
    t0 = find(D.time_min == 0 & D.inducer == "none");
    for i = t0'
        hasIPTG = any(D.group == D.group(i) & D.strain == D.strain(i) & D.csource == D.csource(i) & D.inducer == "IPTG");
        hasT0   = any(D.group == D.group(i) & D.strain == D.strain(i) & D.csource == D.csource(i) & D.inducer == "IPTG" & D.time_min == 0);
        if hasIPTG && ~hasT0
            r = D(i, :); r.inducer = "IPTG"; r.source = "t0 before IPTG (shared)";
            add = [add; r]; %#ok<AGROW>
        end
    end
    D = [D; add];
end

function R = keyRatios(D)
    D = D(D.source ~= "t0 before IPTG (shared)", :);
    R = table('Size', [0 6], 'VariableTypes', {'double','string','string','double','double','string'}, ...
              'VariableNames', {'group','strain','comparison','time_min','ratio','reliability'});
    G = unique(D(:, {'group','strain'}), 'rows');
    for g = 1:height(G)
        d = D(D.group == G.group(g) & D.strain == G.strain(g), :);
        [v, q] = ratioOf(d, "glycerol", "none", 0, "glucose", "none", 0);
        R(end+1, :) = {G.group(g), G.strain(g), "glycerol / glucose (t = 0, catabolite repression)", 0, v, q}; %#ok<AGROW>
        for cs = ["glucose", "glycerol"]
            for t = unique(d.time_min(d.time_min > 0))'
                [v, q] = ratioOf(d, cs, "IPTG", t, cs, "none", t);
                R(end+1, :) = {G.group(g), G.strain(g), "IPTG / no IPTG on " + cs, t, v, q}; %#ok<AGROW>
            end
        end
    end
end

function [v, q] = ratioOf(d, cs1, in1, t1, cs2, in2, t2)
    a = find(d.csource == cs1 & d.inducer == in1 & d.time_min == t1, 1);
    b = find(d.csource == cs2 & d.inducer == in2 & d.time_min == t2, 1);
    if isempty(a) || isempty(b), v = NaN; q = "missing"; return; end
    v = d.activity_U(a) / d.activity_U(b);
    qa = d.quality(a); qb = d.quality(b);
    if qa == "ok" && qb == "ok",                     q = "reliable";
    elseif startsWith(qa, "lower bound") && qb == "ok", q = "lower limit (≥)";
    elseif qa ~= "ok",                              q = "uncertain (numerator: " + qa + ")";
    else,                                           q = "uncertain (denominator: " + qb + ")";
    end
end

function T = pooledKinetics(D)
    ok = D(D.quality == "ok", :);
    [G, strain, csource, inducer, time_min] = findgroups(ok.strain, ok.csource, ok.inducer, ok.time_min);
    n      = splitapply(@numel, ok.activity_U, G);
    mean_U = splitapply(@mean,  ok.activity_U, G);
    sd_U   = splitapply(@(x) std(x, 0), ok.activity_U, G);
    sd_U(n < 2) = NaN;
    T = table(strain, csource, inducer, time_min, n, mean_U, sd_U);
end

%% ======================================================================
%  Plot and display
%% ======================================================================
function plotKinetics(D, outDir)
    strains = unique(D.strain);
    f = figure('Name', 'lac operon kinetics', 'NumberTitle', 'off', ...
               'Position', [80 80 640 * numel(strains) 460]);
    tiledlayout(f, 1, numel(strains));
    col = containers.Map({'glucose','glycerol','lactose'}, {[0.12 0.37 0.75], [0.76 0.25 0.05], [0.18 0.52 0.35]});
    for s = 1:numel(strains)
        ax = nexttile; hold(ax, 'on');
        d = D(D.strain == strains(s), :);
        conds = unique(d(:, {'csource','inducer'}), 'rows');
        h = gobjects(0); names = strings(0);
        for c = 1:height(conds)
            k = d.csource == conds.csource(c) & d.inducer == conds.inducer(c);
            [t, ia] = sort(d.time_min(k)); U = d.activity_U(k); U = U(ia); Q = d.quality(k); Q = Q(ia);
            if isKey(col, char(conds.csource(c))), cc = col(char(conds.csource(c))); else, cc = [0.5 0.5 0.5]; end
            ls = '--'; lab = conds.csource(c) + ", no IPTG";
            if conds.inducer(c) == "IPTG", ls = '-'; lab = conds.csource(c) + ", + IPTG"; end
            h(end+1) = plot(ax, t, U, ls, 'Color', cc, 'LineWidth', 1.8); names(end+1) = lab; %#ok<AGROW>
            for i = 1:numel(t)
                switch Q(i)
                    case "ok",                                   mk = 'o'; fc = cc;
                    case "A420 < 0.02 (near detection limit)",   mk = 'o'; fc = 'w';
                    case "A420 > 1.0 (outside linear range)",    mk = 's'; fc = 'w';
                    case "lower bound (A420 > 3)",               mk = '^'; fc = cc;
                    otherwise,                                   mk = 'd'; fc = 'w';
                end
                plot(ax, t(i), U(i), mk, 'MarkerEdgeColor', cc, 'MarkerFaceColor', fc, 'MarkerSize', 8, 'LineWidth', 1.4);
            end
        end
        % marker legend entries
        mk = {'o','o','s','^'}; fc = {'k','w','w','k'};
        lab = ["ok", "A420 < 0.02 (near detection limit)", "A420 > 1.0 (outside linear range)", "lower bound (A420 > 3)"];
        for m = 1:4
            if any(d.quality == lab(m))
                h(end+1) = plot(ax, NaN, NaN, mk{m}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', fc{m}); %#ok<AGROW>
                names(end+1) = lab(m); %#ok<AGROW>
            end
        end
        set(ax, 'YScale', 'log', 'XTick', unique(d.time_min));
        grid(ax, 'on');
        xlabel(ax, 'time after IPTG addition (min)');
        ylabel(ax, '\beta-galactosidase activity (U)');
        title(ax, sprintf('%s strain (groups: %s)', strains(s), strjoin(string(unique(d.group))', ', ')));
        legend(ax, h, cellstr(names), 'Location', 'best', 'FontSize', 8);
    end
    exportgraphics(f, fullfile(outDir, 'lac_kinetics.png'), 'Resolution', 200);
end

function showTable(T, name)
    for v = T.Properties.VariableNames
        if isnumeric(T.(v{1})), T.(v{1}) = round(T.(v{1}), 2); end
    end
    if height(T) == 0, T = table("none found", 'VariableNames', {'result'}); end
    fig = uifigure('Name', name, 'Position', [120 120 1150 min(700, 70 + 24 * height(T))]);
    g = uigridlayout(fig, [1 1]);
    uitable(g, 'Data', T);
end
