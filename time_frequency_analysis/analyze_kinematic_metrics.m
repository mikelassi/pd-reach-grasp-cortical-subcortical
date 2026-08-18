%% ============================================================
%  Kinematic Metrics — Group-Level Analysis
%  Replicates Vissani et al. 2021 (PD cohort, n = 8)
%
%  Averaging strategy:
%    Level 1: all trials → per-subject mean  (blocks pooled)
%    Level 2: per-subject means → group mean ± SEM
%
%  Reference values (Vissani et al. 2021, PD group, n = 8):
%    Peak wrist velocity  reach : 0.53 ± 0.07 m/s
%    Peak wrist velocity  pull  : 0.57 ± 0.06 m/s
%    Time to peak vel.    reach : 0.50 ± 0.03 s
%    Time to peak vel.    pull  : 0.55 ± 0.03 s
%    Radius of curvature  reach : 0.36 ± 0.04 m
%    Radius of curvature  pull  : 0.38 ± 0.04 m
%    PHA                        : 35.65 ± 3.87 %
%    PCI                        : 0.67 ± 0.08
%    C-score (PD)               : 0.73 ± 0.17     (paper p.2, main text)
%
%  IMPORTANT — these ± values are SEM, not SD.
%    The paper states on p.1 that "all data is reported as average ± standard
%    deviation", but that is demonstrably wrong for the kinematics.  It reports
%    C-score = 0.73 ± 0.17 (n = 8) while also stating (p.2) that those values
%    range from 0.026 to 1.404.  For 8 values spanning that range with mean 0.73
%    the MINIMUM possible SD is 0.368 — more than double the reported 0.17 — so
%    0.17 cannot be an SD.  Read as SEM it implies SD ≈ 0.48, consistent both
%    with the stated range and with the Supp. Fig. 9 scatter (≈ 0.46).
%    Corroboration: our between-subject SDs are 2.1–3.1x the paper's reported ±
%    across 7 of 9 metrics (sqrt(8) = 2.83); read as SEM our dispersion matches
%    theirs closely (PHA 3.89 vs 3.87, radius 0.036 vs 0.04, t-to-peak 0.033 vs 0.03).
%
%  Author: Michael Lassi
%% ============================================================

close all; clear; clc;

%% ========== PARAMETERS ==========
CSV_PATH   = 'G:\Projects\Parkinson_ReachGrasp\Reprocessing\kinematic_metrics.csv';
RESULTS_DIR = 'G:\Projects\Parkinson_ReachGrasp\Reprocessing\RESULTS_kinematics';

if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end

%% ========== REFERENCE VALUES (Vissani et al. 2021, PD group) ==========
% Format: [mean, SEM]  — see header note: the paper's ± values are SEM, not SD.
% Paper PD group n = 8, so the implied between-subject SD is SEM * sqrt(N_PAPER).
N_PAPER = 8;
REF.peak_vel_reach  = [0.53, 0.07];
REF.peak_vel_pull   = [0.57, 0.06];
REF.t_peak_reach    = [0.50, 0.03];
REF.t_peak_pull     = [0.55, 0.03];
REF.radius_reach    = [0.36, 0.04];
REF.radius_pull     = [0.38, 0.04];
REF.PHA             = [35.65, 3.87];
REF.PCI             = [0.67,  0.08];
REF.C_score         = [0.73,  0.17];   % paper p.2 main text (was 0.65 ± 0.46, eyeballed from Supp. Fig. 9)

%% ========== LOAD CSV ==========
T = readtable(CSV_PATH);
fprintf('Loaded %d trial rows from %d subjects.\n', height(T), numel(unique(T.subject)));
disp(T(1:3,:))

subjects = unique(T.subject);
n_subj   = numel(subjects);

%% ========== LEVEL 1: Per-subject means (pool all trials across blocks) ==========
metrics = {'reach_peak_vel_m_s','reach_time_to_peak_vel_s', ...
           'pull_peak_vel_m_s', 'pull_time_to_peak_vel_s', ...
           'reach_median_radius_m','pull_median_radius_m', ...
           'PHA_pct','PCI','C_score','C_tilde', ...
           'dur_reach_s','dur_grasp_s','dur_pull_s', ...
           'n_mvmt_reach','n_mvmt_pull', ...
           'NJ_reach','NJ_pull', ...
           'SI_reach','SI_pull'};

subj_means = array2table(NaN(n_subj, numel(metrics)), 'VariableNames', metrics);
subj_means.subject = subjects;
subj_trials = zeros(n_subj, 1);

for s = 1:n_subj
    rows = strcmp(T.subject, subjects{s});
    subj_trials(s) = sum(rows);
    for m = 1:numel(metrics)
        vals = T.(metrics{m})(rows);
        subj_means.(metrics{m})(s) = mean(vals, 'omitnan');
    end
end

fprintf('\n--- Per-subject trial counts ---\n');
for s = 1:n_subj
    fprintf('  %s : %d trials\n', subjects{s}, subj_trials(s));
end

%% ========== LEVEL 2: Group statistics ==========
grp_mean = varfun(@(x) mean(x,'omitnan'), subj_means, 'InputVariables', metrics);
grp_sem  = varfun(@(x) std(x, 0, 'omitnan') / sqrt(sum(~isnan(x))), ...
                  subj_means, 'InputVariables', metrics);
grp_std  = varfun(@(x) std(x, 0, 'omitnan'), subj_means, 'InputVariables', metrics);

grp_mean.Properties.VariableNames = metrics;
grp_sem.Properties.VariableNames  = metrics;
grp_std.Properties.VariableNames  = metrics;

%% ========== COMPARISON TABLE ==========
metric_labels = { ...
    'Peak vel reach (m/s)',         'reach_peak_vel_m_s',      REF.peak_vel_reach; ...
    'Peak vel pull (m/s)',          'pull_peak_vel_m_s',        REF.peak_vel_pull; ...
    'Time to peak reach (s)',       'reach_time_to_peak_vel_s', REF.t_peak_reach; ...
    'Time to peak pull (s)',        'pull_time_to_peak_vel_s',  REF.t_peak_pull; ...
    'Radius of curvature reach (m)','reach_median_radius_m',   REF.radius_reach; ...
    'Radius of curvature pull (m)', 'pull_median_radius_m',    REF.radius_pull; ...
    'PHA (%)',                      'PHA_pct',                  REF.PHA; ...
    'PCI',                          'PCI',                      REF.PCI; ...
    'C-score',                      'C_score',                  REF.C_score; ...
};

fprintf('\n%s\n', repmat('=',1,88));
fprintf('%-32s  %16s  %18s  %7s %7s\n', ...
    'Metric', 'Ours (mean±SEM)', 'Vissani (mean±SEM)', 'z', 'Agree?');
fprintf('%s\n', repmat('-',1,88));
for r = 1:size(metric_labels,1)
    label   = metric_labels{r,1};
    mname   = metric_labels{r,2};
    ref_mu  = metric_labels{r,3}(1);
    ref_sem = metric_labels{r,3}(2);
    our_mu  = grp_mean.(mname);
    our_sem = grp_sem.(mname);
    % Two-sample z on the group means, combining both standard errors.
    % This is the appropriate test now that the paper's ± is known to be SEM.
    % The previous check (|diff| <= 2*ref) was in effect ±2 SEM ≈ ±0.71 SD,
    % i.e. STRICTER than the 2 SD it claimed to apply — not looser.
    z     = (our_mu - ref_mu) / sqrt(our_sem^2 + ref_sem^2);
    agree = abs(z) < 1.96;
    flag  = '  ✓';
    if ~agree, flag = '  ✗'; end
    fprintf('%-32s  %7.3f ± %6.3f  %8.3f ± %6.3f  %7.2f %s\n', ...
        label, our_mu, our_sem, ref_mu, ref_sem, z, flag);
end
fprintf('%s\n', repmat('=',1,88));
fprintf(['z = two-sample z on group means (combined SEM); |z| < 1.96 = no ' ...
    'significant difference.\nBoxplots show between-subject spread; the shaded ' ...
    'band is the paper''s implied SD (SEM*sqrt(%d)).\n\n'], N_PAPER);

% Also print C-score on its own (C-tilde has no published reference)
fprintf('C-score  : %.3f ± %.3f (mean ± SEM across subjects)\n', ...
    grp_mean.C_score, grp_sem.C_score);
fprintf('C-tilde  : %.3f ± %.3f\n', grp_mean.C_tilde, grp_sem.C_tilde);

%% ========== FIGURES ==========
% Shared style
box_col  = [0.4 0.6 0.9];   % blue for PD
dot_col  = [0.1 0.1 0.7];
fs_ax    = 13;
fs_title = 14;

% Helper: draw a single boxplot panel with reference overlay
function ax = one_box(y_data, y_label, ref_mu, ref_sem, fig_title, box_col, dot_col)
    ax = gca;
    bp = boxchart(ones(size(y_data)), y_data, ...
        'BoxFaceColor', box_col, 'WhiskerLineColor', box_col, ...
        'MarkerStyle', 'o', 'MarkerColor', dot_col);
    hold on;
    % Jittered individual subject points
    jitter = 0.08 * (rand(numel(y_data),1) - 0.5);
    scatter(1 + jitter, y_data, 50, dot_col, 'filled', 'MarkerFaceAlpha', 0.7);
    % Reference band (Vissani et al. 2021).  The boxplot and dots above show
    % the BETWEEN-SUBJECT spread of our per-subject means, so the reference
    % band has to be on that same scale.  The paper's ± is SEM, so convert to
    % the implied SD before shading — shading the raw SEM draws a band about
    % 2.83x too narrow and makes the agreement look far worse than it is.
    if ~isnan(ref_mu)
        ref_sd = ref_sem * sqrt(8);          % paper PD group n = 8
        xpatch = [0.55, 1.45, 1.45, 0.55];
        ypatch = [ref_mu-ref_sd, ref_mu-ref_sd, ref_mu+ref_sd, ref_mu+ref_sd];
        patch(xpatch, ypatch, [0.9 0.5 0.5], 'FaceAlpha', 0.25, 'EdgeColor','none');
        yline(ref_mu, '--r', 'LineWidth', 1.5);
    end
    ylabel(y_label, 'FontSize', 13);
    title(fig_title, 'FontSize', 14);
    xlim([0.5, 1.5]);
    ax.XTick = [];
    grid on; box on;
    hold off;
end

% ---- Figure A: Supplementary Fig. 3a — Peak velocity ----
figure('Name','Peak Wrist Velocity','Position',[100 100 600 480]);
tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');

nexttile;
one_box(subj_means.reach_peak_vel_m_s, 'Velocity (m/s)', ...
    REF.peak_vel_reach(1), REF.peak_vel_reach(2), 'Reach', box_col, dot_col);

nexttile;
one_box(subj_means.pull_peak_vel_m_s, 'Velocity (m/s)', ...
    REF.peak_vel_pull(1), REF.peak_vel_pull(2), 'Pull', box_col, dot_col);

sgtitle('Peak Wrist Velocity  (PD, n=8) — cf. Supp. Fig. 3a', 'FontSize', fs_title);
print(fullfile(RESULTS_DIR,'fig_peak_velocity'), '-dpng', '-r300');

% ---- Figure B: Supplementary Fig. 3b — Time to peak velocity ----
figure('Name','Time to Peak Velocity','Position',[100 100 600 480]);
tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');

nexttile;
one_box(subj_means.reach_time_to_peak_vel_s, 'Time (s)', ...
    REF.t_peak_reach(1), REF.t_peak_reach(2), 'Reach', box_col, dot_col);

nexttile;
one_box(subj_means.pull_time_to_peak_vel_s, 'Time (s)', ...
    REF.t_peak_pull(1), REF.t_peak_pull(2), 'Pull', box_col, dot_col);

sgtitle('Time to Peak Wrist Velocity  (PD, n=8) — cf. Supp. Fig. 3b', 'FontSize', fs_title);
print(fullfile(RESULTS_DIR,'fig_time_to_peak'), '-dpng', '-r300');

% ---- Figure C: Fig. 1c — Radius of curvature ----
figure('Name','Radius of Curvature','Position',[100 100 600 480]);
tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');

nexttile;
one_box(subj_means.reach_median_radius_m, 'Radius (m)', ...
    REF.radius_reach(1), REF.radius_reach(2), 'Reach', box_col, dot_col);

nexttile;
one_box(subj_means.pull_median_radius_m, 'Radius (m)', ...
    REF.radius_pull(1), REF.radius_pull(2), 'Pull', box_col, dot_col);

sgtitle('Median Radius of Curvature  (PD, n=8) — cf. Fig. 1c', 'FontSize', fs_title);
print(fullfile(RESULTS_DIR,'fig_radius_of_curvature'), '-dpng', '-r300');

% ---- Figure D: Supplementary Fig. 3c — PHA ----
figure('Name','Peak Hand Aperture','Position',[100 100 380 480]);
one_box(subj_means.PHA_pct, 'PHA (%)', ...
    REF.PHA(1), REF.PHA(2), 'Peak Hand Aperture  (PD, n=8) — cf. Supp. Fig. 3c', ...
    box_col, dot_col);
print(fullfile(RESULTS_DIR,'fig_PHA'), '-dpng', '-r300');

% ---- Figure E: Supplementary Fig. 3d — PCI ----
figure('Name','Pre-Shape Coordination Index','Position',[100 100 380 480]);
one_box(subj_means.PCI, 'PCI', ...
    REF.PCI(1), REF.PCI(2), 'Pre-Shape Coord. Index  (PD, n=8) — cf. Supp. Fig. 3d', ...
    box_col, dot_col);
print(fullfile(RESULTS_DIR,'fig_PCI'), '-dpng', '-r300');

% ---- Figure F: Fig. 1d — C-score ----
% Reference: 0.73 ± 0.17 (SEM), paper p.2 main text
figure('Name','C-score','Position',[100 100 380 480]);
one_box(subj_means.C_score, 'C-score  [a.u.]', ...
    REF.C_score(1), REF.C_score(2), 'C-score  (PD, n=8) — cf. Fig. 1d', box_col, dot_col);
print(fullfile(RESULTS_DIR,'fig_Cscore'), '-dpng', '-r300');

% ---- Figure G: Phase durations ----
figure('Name','Phase Durations','Position',[100 100 700 480]);
tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
nexttile; one_box(subj_means.dur_reach_s, 'Duration (s)', NaN, NaN, 'Reach (A→C)', box_col, dot_col);
nexttile; one_box(subj_means.dur_grasp_s, 'Duration (s)', NaN, NaN, 'Grasp (C→D)', box_col, dot_col);
nexttile; one_box(subj_means.dur_pull_s,  'Duration (s)', NaN, NaN, 'Pull (D→F)',  box_col, dot_col);
sgtitle('Phase Durations  (PD, n=8)', 'FontSize', fs_title);
print(fullfile(RESULTS_DIR,'fig_phase_durations'), '-dpng', '-r300');

% ---- Figure H: Movement units ----
figure('Name','Movement Units','Position',[100 100 600 480]);
tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile; one_box(subj_means.n_mvmt_reach, 'N sub-peaks', NaN, NaN, 'Reach', box_col, dot_col);
nexttile; one_box(subj_means.n_mvmt_pull,  'N sub-peaks', NaN, NaN, 'Pull',  box_col, dot_col);
sgtitle('Movement Units — velocity sub-peaks  (PD, n=8)', 'FontSize', fs_title);
print(fullfile(RESULTS_DIR,'fig_movement_units'), '-dpng', '-r300');

% ---- Figure I: Normalized dimensionless jerk ----
figure('Name','Normalized Jerk','Position',[100 100 600 480]);
tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile; one_box(subj_means.NJ_reach, 'NJ (a.u.)', NaN, NaN, 'Reach', box_col, dot_col);
nexttile; one_box(subj_means.NJ_pull,  'NJ (a.u.)', NaN, NaN, 'Pull',  box_col, dot_col);
sgtitle('Normalized Dimensionless Jerk  (PD, n=8)  [more negative = smoother]', 'FontSize', fs_title);
print(fullfile(RESULTS_DIR,'fig_norm_jerk'), '-dpng', '-r300');

% ---- Figure J: Straightness index ----
figure('Name','Straightness Index','Position',[100 100 600 480]);
tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile; one_box(subj_means.SI_reach, 'SI (0–1)', NaN, NaN, 'Reach', box_col, dot_col);
nexttile; one_box(subj_means.SI_pull,  'SI (0–1)', NaN, NaN, 'Pull',  box_col, dot_col);
sgtitle('Wrist Path Straightness Index  (PD, n=8)  [1 = perfectly straight]', 'FontSize', fs_title);
print(fullfile(RESULTS_DIR,'fig_straightness'), '-dpng', '-r300');

% Print summary of additional metrics
fprintf('\n--- Additional metrics (no Vissani reference) ---\n');
fprintf('Phase durations (s):\n');
fprintf('  Reach : %.3f ± %.3f\n', grp_mean.dur_reach_s,  grp_sem.dur_reach_s);
fprintf('  Grasp : %.3f ± %.3f\n', grp_mean.dur_grasp_s,  grp_sem.dur_grasp_s);
fprintf('  Pull  : %.3f ± %.3f\n', grp_mean.dur_pull_s,   grp_sem.dur_pull_s);
fprintf('Movement units (N sub-peaks):\n');
fprintf('  Reach : %.2f ± %.2f\n', grp_mean.n_mvmt_reach, grp_sem.n_mvmt_reach);
fprintf('  Pull  : %.2f ± %.2f\n', grp_mean.n_mvmt_pull,  grp_sem.n_mvmt_pull);
fprintf('Normalized dimensionless jerk:\n');
fprintf('  Reach : %.1f ± %.1f\n', grp_mean.NJ_reach,     grp_sem.NJ_reach);
fprintf('  Pull  : %.1f ± %.1f\n', grp_mean.NJ_pull,      grp_sem.NJ_pull);
fprintf('Path straightness index (0–1):\n');
fprintf('  Reach : %.3f ± %.3f\n', grp_mean.SI_reach,     grp_sem.SI_reach);
fprintf('  Pull  : %.3f ± %.3f\n', grp_mean.SI_pull,      grp_sem.SI_pull);

% ---- Figure K: Overview — all metrics per subject (heatmap) ----
metric_short = {'PkVel\newlinereach','PkVel\newlinepull', ...
                'TPk\newlinereach','TPk\newlinepull', ...
                'Rad\newlinereach','Rad\newlinepull', ...
                'PHA','PCI','C̃', ...
                'Dur\newlinereach','Dur\newlinegrasp','Dur\newlinepull', ...
                'MvmtU\newlinereach','MvmtU\newlinepull', ...
                'NJ\newlinereach','NJ\newlinepull', ...
                'SI\newlinereach','SI\newlinepull'};
plot_metrics = {'reach_peak_vel_m_s','pull_peak_vel_m_s', ...
                'reach_time_to_peak_vel_s','pull_time_to_peak_vel_s', ...
                'reach_median_radius_m','pull_median_radius_m', ...
                'PHA_pct','PCI','C_tilde', ...
                'dur_reach_s','dur_grasp_s','dur_pull_s', ...
                'n_mvmt_reach','n_mvmt_pull', ...
                'NJ_reach','NJ_pull', ...
                'SI_reach','SI_pull'};

mat = zeros(n_subj, numel(plot_metrics));
for m = 1:numel(plot_metrics)
    col = subj_means.(plot_metrics{m});
    mn  = min(col); mx = max(col);
    if mx > mn
        mat(:, m) = (col - mn) / (mx - mn);   % 0–1 scaled per metric
    end
end

figure('Name','Per-Subject Overview','Position',[100 100 1400 380]);
imagesc(mat);
colormap(parula); colorbar;
xticks(1:numel(metric_short)); xticklabels(metric_short);
yticks(1:n_subj);              yticklabels(subjects);
title('Per-subject metric values  (min–max scaled per column) — all 18 metrics', 'FontSize', fs_title);
xlabel('Metric'); ylabel('Subject');
set(gca,'FontSize', fs_ax, 'TickLabelInterpreter','tex');
print(fullfile(RESULTS_DIR,'fig_overview_heatmap'), '-dpng', '-r300');

fprintf('\nAll figures saved to %s\n', RESULTS_DIR);

%% ========== SAVE SUBJECT-LEVEL MEANS ==========
writetable(subj_means, fullfile(RESULTS_DIR, 'kinematic_subj_means.csv'));
fprintf('Subject-level means saved.\n');
