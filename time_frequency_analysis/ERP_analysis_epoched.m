%% ERP_analysis_epoched.m
% =========================================================================
% Event-Related Potential (ERP) analysis of STN LFP epochs.
%
% Input:
%   Epochs_allSubjects.mat  (from Epoching.m)
%   Contralateral-hemisphere channel, fixed-length epochs:
%     MovOnset : [-0.5, +1.0] s around event A
%     Pull     : [-0.3, +0.875] s around event D
%
% What the ERP captures for LFP:
%   The trial-average retains only the *phase-locked* component of the signal.
%   Beta oscillations cancel in the average (not phase-locked to the event).
%   What remains is mainly:
%     - Slow potential shifts < 4 Hz (subcortical analog of the readiness
%       potential / movement-related cortical potential)
%     - Sharp transients locked to movement events
%   Both broadband (unfiltered) and LP-filtered (< 4 Hz) views are shown.
%
% Outputs (saved to RESULTS_DIR/ERP/):
%   - ERP_grandAvg.png       grand average both epoch types
%   - ERP_perSubject.png     per-subject grid
%   - ERP_allSubjects.mat    ERP matrices + metadata
%
% Author: Michael Lassi
% =========================================================================

clear; clc;

%% ===== PARAMETERS =====
MAT_FILE    = 'H:\Parkinson_ReachGrasp\Reprocessing\RESULTS_final\Epochs\Epochs_allSubjects.mat';
RESULTS_DIR = 'H:\Parkinson_ReachGrasp\Reprocessing\RESULTS_final\ERP';

LP_CUTOFF_HZ = 4;    % low-pass cutoff for slow-potential view (Hz)
LP_ORDER     = 4;    % Butterworth filter order

% Baseline windows (s) — pre-event periods used for mean subtraction
BASELINE_ONSET = [-0.5, 0];   % rest window in onset epoch
BASELINE_PULL  = [-0.3, 0];   % late-grasp window in pull epoch

% Colours
COL_ONSET  = [0.2 0.5 0.9];   % blue  — movement onset
COL_PULL   = [0.8 0.3 0.1];   % red   — pull onset
COL_SHADE  = 0.35;             % alpha for SEM shading

% Biosig ships a stub butter.m in "maybe-missing/" that shadows the
% Signal Processing Toolbox's butter — remove it before filtering.
biosig_stubs = fullfile(fileparts(which('eeglab')), ...
    'plugins', 'Biosig3.8.4', 'biosig', 'maybe-missing');
if exist(biosig_stubs, 'dir'), rmpath(biosig_stubs); end

%% ===== LOAD =====
fprintf('Loading epochs...\n');
load(MAT_FILE, 'EPOCHS');
SUBJECTS = EPOCHS.params.subjects;
n_subj   = numel(SUBJECTS);

if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end

%% ===== COMPUTE ERPs PER SUBJECT =====
% Storage
erp_onset_subj = cell(1, n_subj);   % per-subject ERPs [1 × n_samp]
erp_pull_subj  = cell(1, n_subj);
erp_onset_lp   = cell(1, n_subj);   % low-pass filtered
erp_pull_lp    = cell(1, n_subj);
t_onset = [];  t_pull = [];
fs = 400;

for s = 1:n_subj
    subj = SUBJECTS{s};
    if ~isfield(EPOCHS, subj), continue; end

    fs      = EPOCHS.(subj).fs;
    t_on    = EPOCHS.(subj).onset.times;   % [1 × 600]
    t_pu    = EPOCHS.(subj).pull.times;    % [1 × 470]
    data_on = EPOCHS.(subj).onset.data;   % [600 × n_trials]
    data_pu = EPOCHS.(subj).pull.data;    % [470 × n_trials]

    if isempty(t_onset), t_onset = t_on; t_pull = t_pu; end

    % --- Baseline correction per trial ---
    bl_on_mask = t_on >= BASELINE_ONSET(1) & t_on < BASELINE_ONSET(2);
    bl_pu_mask = t_pu >= BASELINE_PULL(1)  & t_pu < BASELINE_PULL(2);

    data_on_bc = data_on - mean(data_on(bl_on_mask, :), 1);  % subtract trial-mean
    data_pu_bc = data_pu - mean(data_pu(bl_pu_mask, :), 1);

    % --- Trial average → ERP ---
    erp_onset_subj{s} = mean(data_on_bc, 2)';   % [1 × 600]
    erp_pull_subj{s}  = mean(data_pu_bc, 2)';

    % --- LP filter for slow-potential view ---
    [b, a] = butter(LP_ORDER, LP_CUTOFF_HZ / (fs/2), 'low');
    erp_onset_lp{s} = filtfilt(b, a, erp_onset_subj{s});
    erp_pull_lp{s}  = filtfilt(b, a, erp_pull_subj{s});
end

%% ===== GRAND AVERAGE =====
onset_mat = cell2mat(erp_onset_subj');   % [n_subj × 600]
pull_mat  = cell2mat(erp_pull_subj');

onset_lp_mat = cell2mat(erp_onset_lp');
pull_lp_mat  = cell2mat(erp_pull_lp');

ga_onset     = mean(onset_mat, 1);
ga_pull      = mean(pull_mat,  1);
ga_onset_lp  = mean(onset_lp_mat, 1);
ga_pull_lp   = mean(pull_lp_mat,  1);

sem_onset    = std(onset_mat,    0, 1) / sqrt(n_subj);
sem_pull     = std(pull_mat,     0, 1) / sqrt(n_subj);
sem_onset_lp = std(onset_lp_mat, 0, 1) / sqrt(n_subj);
sem_pull_lp  = std(pull_lp_mat,  0, 1) / sqrt(n_subj);

%% ===== FIGURE 1: Grand Average ERP =====
figure('Name','Grand Average ERP','Position',[50 50 1300 550],'Visible','off');

% --- Helper for shaded SEM ---
shade = @(ax,t,mn,se,col) fill(ax, [t, fliplr(t)], [mn+se, fliplr(mn-se)], ...
    col, 'EdgeColor','none','FaceAlpha',COL_SHADE,'HandleVisibility','off');

subplot(1,2,1); hold on; grid on;
shade(gca, t_onset, ga_onset,    sem_onset,    COL_ONSET);
shade(gca, t_onset, ga_onset_lp, sem_onset_lp, [0 0 0]);
plot(t_onset, ga_onset,    'Color', COL_ONSET, 'LineWidth', 2,   'DisplayName', 'Broadband ERP');
plot(t_onset, ga_onset_lp, 'Color', [0 0 0],   'LineWidth', 2.2, 'DisplayName', sprintf('LP < %d Hz', LP_CUTOFF_HZ));
xline(0, '--k', 'Mov onset (A)', 'LabelVerticalAlignment','bottom','LineWidth',1.5);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Amplitude (z-score re baseline)', 'FontSize', 12);
title('Movement Onset ERP (grand avg ± SEM)', 'FontSize', 13);
legend('Location','best','FontSize',10);
set(gca,'FontSize',11); xlim([t_onset(1) t_onset(end)]);

subplot(1,2,2); hold on; grid on;
shade(gca, t_pull, ga_pull,    sem_pull,    COL_PULL);
shade(gca, t_pull, ga_pull_lp, sem_pull_lp, [0 0 0]);
plot(t_pull, ga_pull,    'Color', COL_PULL, 'LineWidth', 2,   'DisplayName', 'Broadband ERP');
plot(t_pull, ga_pull_lp, 'Color', [0 0 0],  'LineWidth', 2.2, 'DisplayName', sprintf('LP < %d Hz', LP_CUTOFF_HZ));
xline(0, '--k', 'Pull onset (D)', 'LabelVerticalAlignment','bottom','LineWidth',1.5);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Amplitude (z-score re baseline)', 'FontSize', 12);
title('Pull Onset ERP (grand avg ± SEM)', 'FontSize', 13);
legend('Location','best','FontSize',10);
set(gca,'FontSize',11); xlim([t_pull(1) t_pull(end)]);

sgtitle('Grand Average LFP ERP — STN Contralateral Hemisphere  (n=8)', 'FontSize', 14);
set(gcf,'Units','pixels','Position',[50 50 1300 550]);
print(fullfile(RESULTS_DIR, 'ERP_grandAvg'), '-dpng', '-r300');
savefig(fullfile(RESULTS_DIR, 'ERP_grandAvg.fig'));

%% ===== FIGURE 2: Per-Subject ERPs — Movement Onset =====
figure('Name','Per-Subject ERP — MovOnset','Position',[50 50 1400 700],'Visible','off');
for s = 1:n_subj
    subj = SUBJECTS{s};
    subplot(2,4,s); hold on; grid on;
    if isempty(erp_onset_subj{s}), title([subj ' — MISSING']); continue; end
    plot(t_onset, erp_onset_subj{s}, 'Color', [COL_ONSET 0.5], 'LineWidth', 1.2, ...
        'DisplayName', 'Broadband');
    plot(t_onset, erp_onset_lp{s},   'Color', [0 0 0],         'LineWidth', 2.0, ...
        'DisplayName', sprintf('LP<%dHz',LP_CUTOFF_HZ));
    xline(0,'--k','LineWidth',1.2);
    title(subj,'FontSize',11);
    xlabel('Time (s)','FontSize',9); ylabel('\muV','FontSize',9);
    set(gca,'FontSize',9); xlim([t_onset(1) t_onset(end)]);
    if s == 1, legend('Location','best','FontSize',8); end
end
sgtitle('Per-Subject LFP ERP — Movement Onset', 'FontSize', 13);
set(gcf,'Units','pixels','Position',[50 50 1400 700]);
print(fullfile(RESULTS_DIR, 'ERP_perSubject_MovOnset'), '-dpng', '-r300');

%% ===== FIGURE 3: Per-Subject ERPs — Pull =====
figure('Name','Per-Subject ERP — Pull','Position',[50 50 1400 700],'Visible','off');
for s = 1:n_subj
    subj = SUBJECTS{s};
    subplot(2,4,s); hold on; grid on;
    if isempty(erp_pull_subj{s}), title([subj ' — MISSING']); continue; end
    plot(t_pull, erp_pull_subj{s}, 'Color', [COL_PULL 0.5], 'LineWidth', 1.2, ...
        'DisplayName', 'Broadband');
    plot(t_pull, erp_pull_lp{s},   'Color', [0 0 0],        'LineWidth', 2.0, ...
        'DisplayName', sprintf('LP<%dHz',LP_CUTOFF_HZ));
    xline(0,'--k','LineWidth',1.2);
    title(subj,'FontSize',11);
    xlabel('Time (s)','FontSize',9); ylabel('\muV','FontSize',9);
    set(gca,'FontSize',9); xlim([t_pull(1) t_pull(end)]);
    if s == 1, legend('Location','best','FontSize',8); end
end
sgtitle('Per-Subject LFP ERP — Pull Onset', 'FontSize', 13);
set(gcf,'Units','pixels','Position',[50 50 1400 700]);
print(fullfile(RESULTS_DIR, 'ERP_perSubject_Pull'), '-dpng', '-r300');

%% ===== SAVE =====
ERP.subjects       = SUBJECTS;
ERP.fs             = fs;
ERP.t_onset        = t_onset;
ERP.t_pull         = t_pull;
ERP.grand_avg_onset      = ga_onset;
ERP.grand_avg_pull       = ga_pull;
ERP.grand_avg_onset_lp   = ga_onset_lp;
ERP.grand_avg_pull_lp    = ga_pull_lp;
ERP.sem_onset            = sem_onset;
ERP.sem_pull             = sem_pull;
ERP.per_subj_onset       = onset_mat;   % [n_subj × samp]
ERP.per_subj_pull        = pull_mat;
ERP.per_subj_onset_lp    = onset_lp_mat;
ERP.per_subj_pull_lp     = pull_lp_mat;
ERP.params.lp_cutoff_hz  = LP_CUTOFF_HZ;
ERP.params.baseline_onset = BASELINE_ONSET;
ERP.params.baseline_pull  = BASELINE_PULL;

save(fullfile(RESULTS_DIR, 'ERP_allSubjects.mat'), 'ERP', '-v7.3');
fprintf('Saved ERP results to %s\n', RESULTS_DIR);
