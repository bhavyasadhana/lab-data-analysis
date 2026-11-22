%% IL-2 sandwich ELISA: standard-curve fitting and back-calculation
%
% Fits three calibration models to the IL-2 standards (0-600 pg/ml):
%   linear           OD = m*c + b
%   semi-log linear  OD = m*log10(c) + b         (standards > 0 only)
%   4PL              OD = d + (a - d) ./ (1 + (c/C).^B)
% then back-calculates the samples, checks replicate CVs and standard
% recovery, and estimates LOD / LOQ.
%
% Base MATLAB only (R2016b or newer): the 4PL is fitted with fminsearch,
% so no Optimization or Statistics Toolbox is required.
%
% Run from the matlab/ folder:   >> elisa_curve

clear; close all force;                       % also closes table windows from the last run

dilution = 100;                               % 10 µl + 1000 µl (strictly 1:101)
here     = fileparts(mfilename('fullpath'));
dataDir  = fullfile(here, '..', 'data');
outDir   = fullfile(here, '..', 'results', 'matlab');
if ~exist(outDir, 'dir'), mkdir(outDir); end

stdAll = readtable(fullfile(dataDir, 'standards.csv'));
smpAll = readtable(fullfile(dataDir, 'samples.csv'), 'TextType', 'string');
stdAll.person = string(stdAll.person);

persons = unique(stdAll.person);
sampleRows = {}; recRows = {}; fitRows = {};

for k = 1:numel(persons)
    p   = persons(k);
    S   = stdAll(stdAll.person == p, :);
    P   = smpAll(smpAll.person == p, :);

    % ---- assemble replicate points (blank only if both wells filled) ----
    blankRow = S(S.conc_pgml == 0, :);
    haveBlank = ~isempty(blankRow) && all(~isnan([blankRow.od1 blankRow.od2]));
    U = S;
    if ~haveBlank, U = S(S.conc_pgml > 0, :); end
    conc = [U.conc_pgml; U.conc_pgml];
    od   = [U.od1; U.od2];

    % ---- fit models ------------------------------------------------------
    M = fitModels(conc, od);
    L = detectionLimits(S, M);

    % ---- samples -----------------------------------------------------------
    for i = 1:height(P)
        meanOD = mean([P.od1(i) P.od2(i)]);
        c4 = M.fpl.inv(meanOD);
        if isnan(c4),            status = "outside asymptotes";
        elseif c4 > 600,         status = "above top standard";
        elseif c4 < L.lod,       status = "below LOD";
        elseif c4 < L.loq,       status = "below LOQ";
        else,                    status = "in range";
        end
        sampleRows(end+1, :) = {p, P.sample(i), meanOD, cvPct(P.od1(i), P.od2(i)), ...
            c4, M.fpl.inv(P.od1(i)), M.fpl.inv(P.od2(i)), M.semi.inv(meanOD), ...
            M.lin.inv(meanOD), c4 * dilution / 1000, ...
            mean([P.software_pgml_1(i) P.software_pgml_2(i)]), status}; %#ok<SAGROW>
    end

    % ---- standard recovery -------------------------------------------------
    R = S(S.conc_pgml > 0, :);
    for i = 1:height(R)
        meanOD = mean([R.od1(i) R.od2(i)]);
        c = R.conc_pgml(i);
        recRows(end+1, :) = {p, c, meanOD, cvPct(R.od1(i), R.od2(i)), ...
            M.fpl.inv(meanOD)/c*100, M.semi.inv(meanOD)/c*100, M.lin.inv(meanOD)/c*100}; %#ok<SAGROW>
    end

    % ---- fit summary -------------------------------------------------------
    names = {'4pl', 'semilog', 'linear'}; mods = {M.fpl, M.semi, M.lin};
    for j = 1:3
        fitRows(end+1, :) = {p, string(names{j}), mods{j}.r2, mods{j}.rmse, ...
            mods{j}.info, L.lod, L.loq, L.lod4pl, L.method}; %#ok<SAGROW>
    end

    plotCurve(p, conc, od, M, P, L, outDir);
end

samples = cell2table(sampleRows, 'VariableNames', {'person','sample','mean_od', ...
    'od_cv_pct','pgml_4pl','pgml_4pl_rep1','pgml_4pl_rep2','pgml_semilog', ...
    'pgml_linear','supernatant_ngml_4pl','software_pgml_mean','status'});
recovery = cell2table(recRows, 'VariableNames', {'person','conc_pgml','mean_od', ...
    'od_cv_pct','recovery_4pl_pct','recovery_semilog_pct','recovery_linear_pct'});
fits = cell2table(fitRows, 'VariableNames', {'person','model','r2','rmse_od', ...
    'parameters','lod_pgml','loq_pgml','lod_4pl_pgml','lod_method'});

writetable(samples,  fullfile(outDir, 'sample_results.csv'));
writetable(recovery, fullfile(outDir, 'standard_recovery.csv'));
writetable(fits,     fullfile(outDir, 'fit_summary.csv'));

disp('Fit summary'); disp(fits(:, {'person','model','r2','rmse_od','lod_pgml','loq_pgml'}));
disp('Samples (in-well pg/ml unless noted)');
disp(samples(:, {'person','sample','mean_od','od_cv_pct','pgml_4pl','pgml_semilog', ...
    'pgml_linear','supernatant_ngml_4pl','software_pgml_mean','status'}));
fprintf('Results written to %s\n', outDir);

% ---- pop-up windows with the result tables ------------------------------
showTable(samples, 'IL-2 sample results (in-well pg/ml, supernatant ng/ml)');
showTable(recovery, 'Standard recovery (%)');
showTable(fits, 'Fit quality and detection limits');


%% ======================================================================
%  Local functions
%% ======================================================================
function y = fourPL(x, p)
    % p = [a B C d]; a = OD at zero, d = plateau, C = EC50, B = slope
    y = p(4) + (p(1) - p(4)) ./ (1 + (x ./ p(3)).^p(2));
end

function x = fourPLinv(y, p)
    lo = min(p(1), p(4)); hi = max(p(1), p(4));
    x = p(3) .* ((p(1) - p(4)) ./ (y - p(4)) - 1).^(1 ./ p(2));
    x(y <= lo | y >= hi) = NaN;
    x = real(x);
end

function r2 = rsq(y, yhat)
    r2 = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
end

function M = fitModels(conc, od)
    % Linear
    c1 = polyfit(conc, od, 1);
    yhat = polyval(c1, conc);
    M.lin.fwd  = @(c) polyval(c1, c);
    M.lin.inv  = @(y) (y - c1(2)) / c1(1);
    M.lin.r2   = rsq(od, yhat);
    M.lin.rmse = sqrt(sum((od - yhat).^2) / (numel(od) - 2));
    M.lin.info = sprintf('m=%.5g, b=%.4g', c1(1), c1(2));

    % Semi-log linear
    pos = conc > 0;
    cs = polyfit(log10(conc(pos)), od(pos), 1);
    yhat = polyval(cs, log10(conc(pos)));
    M.semi.fwd  = @(c) polyval(cs, log10(c));
    M.semi.inv  = @(y) 10.^((y - cs(2)) / cs(1));
    M.semi.r2   = rsq(od(pos), yhat);
    M.semi.rmse = sqrt(sum((od(pos) - yhat).^2) / (nnz(pos) - 2));
    M.semi.info = sprintf('m=%.4g, b=%.4g', cs(1), cs(2));

    % 4PL: minimise SSE with fminsearch. B and C are fitted on a log scale
    % so they stay positive.
    q0  = [min(od), log(1), mean(log(conc(pos))), max(od) * 1.15];
    unpack = @(q) [q(1), exp(q(2)), exp(q(3)), q(4)];
    sse = @(q) sum((od - fourPL(conc, unpack(q))).^2);
    opts = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5, 'TolX', 1e-10, 'TolFun', 1e-12);
    q = fminsearch(sse, q0, opts);
    q = fminsearch(sse, q, opts);            % restart once to escape simplex stalls
    p = unpack(q);
    yhat = fourPL(conc, p);
    M.fpl.p    = p;
    M.fpl.fwd  = @(c) fourPL(c, p);
    M.fpl.inv  = @(y) fourPLinv(y, p);
    M.fpl.r2   = rsq(od, yhat);
    M.fpl.rmse = sqrt(sum((od - yhat).^2) / (numel(od) - 4));
    M.fpl.info = sprintf('a=%.4g, B=%.4g, C=%.4g pg/ml, d=%.4g', p);
end

function L = detectionLimits(S, M)
    blank = S{S.conc_pgml == 0, {'od1', 'od2'}};
    blank = blank(~isnan(blank));
    if numel(blank) >= 2
        L.lod = M.fpl.inv(mean(blank) + 3 * std(blank));
        L.loq = M.fpl.inv(mean(blank) + 10 * std(blank));
        L.method = "blank + 3 SD / 10 SD (4PL)";
    else
        low = S(S.conc_pgml > 0 & S.conc_pgml <= 50, :);
        x = [low.conc_pgml; low.conc_pgml];
        y = [low.od1; low.od2];
        c = polyfit(x, y, 1);
        s = sqrt(sum((y - polyval(c, x)).^2) / (numel(x) - 2));
        L.lod = 3.3 * s / c(1);
        L.loq = 10 * s / c(1);
        L.method = "ICH Q2: 3.3*s/slope, 10*s/slope (standards <= 50 pg/ml)";
    end
    L.lod4pl = M.fpl.inv(M.fpl.p(1) + 3 * M.fpl.rmse);
end

function v = cvPct(a, b)
    v = std([a b]) / mean([a b]) * 100;
end

function plotCurve(person, conc, od, M, P, L, outDir)
    f = figure('Name', sprintf('IL-2 standard curve - Person %s', person), ...
        'NumberTitle', 'off', 'Position', [100 100 760 500]);
    xx = logspace(log10(2.2), log10(800), 300);
    hold on;
    patch([2.2 L.lod L.lod 2.2], [0 0 3 3], [0.5 0.5 0.5], ...
        'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', sprintf('< LOD (%.1f pg/ml)', L.lod));
    plot(xx, M.fpl.fwd(xx),  '-',  'LineWidth', 2.2, 'Color', [0.12 0.37 0.75], 'DisplayName', '4PL');
    plot(xx, M.semi.fwd(xx), '--', 'LineWidth', 1.6, 'Color', [0.76 0.25 0.05], 'DisplayName', 'Semi-log linear');
    plot(xx, M.lin.fwd(xx),  ':',  'LineWidth', 1.6, 'Color', [0.42 0.48 0.55], 'DisplayName', 'Linear');
    pos = conc > 0;
    plot(conc(pos), od(pos), 's', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 7, 'DisplayName', 'Standards');
    for i = 1:height(P)
        y = mean([P.od1(i) P.od2(i)]);
        x = M.fpl.inv(y);
        if isfinite(x)
            plot([2.2 x x], [y y 0], ':', 'Color', [0.54 0.42 0], 'HandleVisibility', 'off');
            plot(x, y, 'o', 'MarkerFaceColor', [0.95 0.76 0], 'MarkerEdgeColor', [0.54 0.42 0], ...
                'MarkerSize', 8, 'HandleVisibility', 'off');
            text(x * 1.08, y + 0.08, P.sample(i));
        end
    end
    hold off;
    set(gca, 'XScale', 'log', 'XLim', [2.2 800], 'XTick', [3.1 12.5 50 300 600]);
    ylim([0 max(od) * 1.15]);
    xlabel('IL-2 (pg/ml, log scale)'); ylabel('OD 450 nm');
    title(sprintf('IL-2 standard curve, Person %s  (4PL R^2 = %.3f)', person, M.fpl.r2));
    grid on; legend('Location', 'northwest');
    exportgraphics_or_print(f, fullfile(outDir, sprintf('curve_person_%s.png', person)));
    % figure stays open so it shows up in MATLAB
end

function showTable(T, name)
    % Opens a window with the table, numbers rounded to 3 decimals.
    for v = T.Properties.VariableNames
        if isnumeric(T.(v{1})), T.(v{1}) = round(T.(v{1}), 3); end
    end
    fig = uifigure('Name', name, 'Position', [120 120 1100 60 + 26 * height(T)]);
    g = uigridlayout(fig, [1 1]);
    uitable(g, 'Data', T);
end

function exportgraphics_or_print(f, file)
    if exist('exportgraphics', 'file')      % R2020a+
        exportgraphics(f, file, 'Resolution', 200);
    else
        print(f, file, '-dpng', '-r200');
    end
end
