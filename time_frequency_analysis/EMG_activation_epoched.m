%% EMG_activation_epoched.m
% =========================================================================
% Event-locked average EMG activation analysis.
%
% Input:
%   H:\...\<subj>\01_Extracted\EMG_KIN\*.set  (raw EMG with events)
%
% What it does:
%   For each subject:
%     1. Load continuous EMG recordings block by block.
%     2. Select contralateral muscles (IOD, Triceps, Deltoid).
%     3. Preprocess: demean → BP 20–180 Hz → notch → rectify → LP envelope.
%     4. Normalise envelope to 95th percentile of the continuous signal.
%     5. Epoch around event A (movement onset) and event D (pull onset),
%        same windows as Epoching.m.
%     6. Subtract pre-event baseline mean per trial.
%     7. Average across trials → per-subject EMG activation.
%   Grand average + SEM across subjects.
%   Figures: grand-average overlay (all muscles), per-subject grid per muscle.
%
% Outputs (saved to RESULTS_DIR/EMG_Activation/):
%   EMG_grandAvg_activation.png / .fig
%   EMG_perSubject_MovOnset_<muscle>.png  (one per muscle)
%   EMG_perSubject_Pull_<muscle>.png
%   EMG_activation.mat
%
% Author: Michael Lassi
% =========================================================================

clear; clc;

%% ===== PARAMETERS =====
SUBJECTS    = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
BASE_PATH   = 'H:\Parkinson_ReachGrasp\Reprocessing';
INPUT_SUB   = fullfile('01_Extracted', 'EMG_KIN');
LFP_SUB     = fullfile('Preprocessed', 'LFP');   % source of kinematic events (*_wEv.set)
RESULTS_DIR = fullfile(BASE_PATH, 'RESULTS_final', 'EMG_Activation');

% Epoch windows — identical to Epoching.m
PRE_ONSET_S  = 0.5;
POST_ONSET_S = 1.0;
PRE_PULL_S   = 0.3;
POST_PULL_S  = 0.875;

% Subjects with dominant / non-dominant EMG sets
DOMINANT_SUBJECTS = {'wue02','wue03','wue10','wue11'};

% Target muscles — must be partial-match strings for channel labels
TARGET_MUSCLES = {'iod', 'triceps', 'deltoid'};
MUSCLE_LABELS  = {'IOD', 'Triceps', 'Deltoid'};
n_muscles      = numel(TARGET_MUSCLES);

% Preprocessing
BP_LOW   = 20;  BP_HIGH = 180;  BP_ORDER = 4;
NOTCH_FREQS = [50, 100];  NOTCH_BW = 2;  NOTCH_ORDER = 2;
ENV_CUTOFF  = 5;  ENV_ORDER = 4;
NORM_PCTILE = 95;

% Baseline window (s) — pre-event rest period
BASELINE_WIN_ONSET = [-0.5, 0];
BASELINE_WIN_PULL  = [-0.3, 0];

% Display
MUSCLE_COLS = [0.15 0.65 0.25;   % IOD     — green
               0.80 0.30 0.10;   % Triceps — red-orange
               0.20 0.40 0.80];  % Deltoid — blue
COL_SHADE   = 0.25;

%% ===== BIOSIG STUB FIX =====
if ~isempty(which('eeglab'))
    biosig_stubs = fullfile(fileparts(which('eeglab')), ...
        'plugins', 'Biosig3.8.4', 'biosig', 'maybe-missing');
    if exist(biosig_stubs, 'dir'), rmpath(biosig_stubs); end
end

%% ===== INIT =====
eeglab nogui;
if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end

n_subj    = numel(SUBJECTS);
t_onset   = [];
t_pull    = [];

% Per-subject averaged activations: {n_muscles × n_subj}, each [1 × n_samp]
emg_onset_subj = cell(n_muscles, n_subj);
emg_pull_subj  = cell(n_muscles, n_subj);

%% ===== SUBJECT LOOP =====
for s = 1:n_subj
    subj     = SUBJECTS{s};
    emg_path = fullfile(BASE_PATH, subj, INPUT_SUB);
    set_files = dir(fullfile(emg_path, '*.set'));

    if isempty(set_files)
        warning('No .set files for %s — skipping.', subj);
        continue;
    end

    fprintf('\n===== Subject: %s  (%d block(s)) =====\n', subj, numel(set_files));
    use_dominant = ismember(subj, DOMINANT_SUBJECTS);

    % LFP files are the authoritative source of kinematic events —
    % EMG_KIN .set files contain signal only, no A_T*/D_T* markers.
    lfp_path  = fullfile(BASE_PATH, subj, LFP_SUB);
    lfp_files = dir(fullfile(lfp_path, '*_wEv.set'));
    if isempty(lfp_files)
        warning('No LFP *_wEv.set for %s — cannot source events, skipping.', subj);
        continue;
    end

    % Accumulators: each cell{m} = [n_samp × n_trials_total]
    ep_onset = cell(n_muscles, 1);
    ep_pull  = cell(n_muscles, 1);
    fs_subj  = NaN;

    for f = 1:numel(set_files)
        EEG = pop_loadset('filename', set_files(f).name, 'filepath', emg_path);
        fs  = EEG.srate;
        if isnan(fs_subj), fs_subj = fs; end

        % ---- Inject kinematic events from preprocessed LFP file ----
        % EMG_KIN .set files carry no A_T*/D_T* event markers; copy them
        % from the matching LFP block, converting latencies via fs ratio.
        if f <= numel(lfp_files)
            LFP_ev   = pop_loadset('filename', lfp_files(f).name, 'filepath', lfp_path);
            kin_mask = ~cellfun(@isempty, regexp({LFP_ev.event.type}, '^[A-F]_T\d+$'));
            kin_evs  = LFP_ev.event(kin_mask);
            if ~isempty(kin_evs)
                fs_lfp = LFP_ev.srate;
                for k = 1:numel(kin_evs)
                    kin_evs(k).latency = round(kin_evs(k).latency * fs / fs_lfp);
                end
                EEG.event = [EEG.event, kin_evs];
                EEG        = eeg_checkset(EEG);
            end
        else
            warning('Block %d: no matching LFP file for events — skipping block.', f);
            continue;
        end

        % ---- Select channels ----
        ch_labels = lower({EEG.chanlocs.labels});
        if use_dominant
            sel_mask = contains(ch_labels, '_dominant_emg');
        else
            sel_mask = contains(ch_labels, '_emg');
        end
        ch_all = find(sel_mask);
        if isempty(ch_all), continue; end

        % Map each target muscle to a channel index (first match)
        muscle_ch = zeros(1, n_muscles);
        for m = 1:n_muscles
            hits = find(contains(ch_labels(ch_all), TARGET_MUSCLES{m}));
            if ~isempty(hits), muscle_ch(m) = ch_all(hits(1)); end
        end

        % ---- Build filters once per block (fs may vary) ----
        [b_bp,  a_bp]  = butter(BP_ORDER,  [BP_LOW, BP_HIGH]/(fs/2), 'bandpass');
        [b_env, a_env]  = butter(ENV_ORDER, ENV_CUTOFF/(fs/2),        'low');
        notch_filt = cell(numel(NOTCH_FREQS), 2);
        for ni = 1:numel(NOTCH_FREQS)
            fn = NOTCH_FREQS(ni);
            if fn < fs/2
                bw = NOTCH_BW/(fs/2); fc = fn/(fs/2);
                [notch_filt{ni,1}, notch_filt{ni,2}] = ...
                    butter(NOTCH_ORDER, [fc-bw/2, fc+bw/2], 'stop');
            end
        end

        % ---- Preprocess each muscle ----
        env_data = zeros(n_muscles, EEG.pnts);
        for m = 1:n_muscles
            if muscle_ch(m) == 0, continue; end
            sig = double(EEG.data(muscle_ch(m), :));
            sig = sig - mean(sig);
            sig = filtfilt(b_bp, a_bp, sig);
            for ni = 1:numel(NOTCH_FREQS)
                if ~isempty(notch_filt{ni,1})
                    sig = filtfilt(notch_filt{ni,1}, notch_filt{ni,2}, sig);
                end
            end
            sig = abs(sig);
            sig = filtfilt(b_env, a_env, sig);
            p95 = prctile(sig, NORM_PCTILE);
            if p95 > 0, sig = sig / p95; end
            env_data(m, :) = sig;
        end

        % ---- Build minimal EEG struct with processed data for pop_epoch ----
        EEG_proc          = EEG;
        EEG_proc.data     = env_data;
        EEG_proc.nbchan   = n_muscles;
        for m = 1:n_muscles
            EEG_proc.chanlocs(m).labels = MUSCLE_LABELS{m};
        end
        EEG_proc.chanlocs = EEG_proc.chanlocs(1:n_muscles);

        % ---- Epoch around A (movement onset) ----
        ev_types = {EEG.event.type};
        types_A  = unique(ev_types(~cellfun(@isempty, regexp(ev_types, '^A_T\d+$'))));
        types_D  = unique(ev_types(~cellfun(@isempty, regexp(ev_types, '^D_T\d+$'))));

        if ~isempty(types_A)
            EEG_ep = pop_epoch(EEG_proc, types_A, [-PRE_ONSET_S, POST_ONSET_S]);
            if isempty(t_onset), t_onset = EEG_ep.times / 1000; end
            for m = 1:n_muscles
                if size(EEG_ep.data,3) > 0
                    trials = squeeze(EEG_ep.data(m, :, :));  % [samp × trials]
                    if isvector(trials), trials = trials(:); end
                    ep_onset{m} = [ep_onset{m}, trials];
                end
            end
        end

        if ~isempty(types_D)
            EEG_ep = pop_epoch(EEG_proc, types_D, [-PRE_PULL_S, POST_PULL_S]);
            if isempty(t_pull), t_pull = EEG_ep.times / 1000; end
            for m = 1:n_muscles
                if size(EEG_ep.data,3) > 0
                    trials = squeeze(EEG_ep.data(m, :, :));
                    if isvector(trials), trials = trials(:); end
                    ep_pull{m} = [ep_pull{m}, trials];
                end
            end
        end

    end  % block loop

    % ---- Baseline correction (per trial) and trial average ----
    if ~isempty(t_onset)
        bl_on = t_onset >= BASELINE_WIN_ONSET(1) & t_onset < BASELINE_WIN_ONSET(2);
        bl_pu = t_pull  >= BASELINE_WIN_PULL(1)  & t_pull  < BASELINE_WIN_PULL(2);

        for m = 1:n_muscles
            if ~isempty(ep_onset{m})
                bl = mean(ep_onset{m}(bl_on, :), 1);
                ep_onset{m} = ep_onset{m} - bl;
                emg_onset_subj{m,s} = mean(ep_onset{m}, 2)';   % [1 × n_samp]
                fprintf('  %s onset: %d trials\n', MUSCLE_LABELS{m}, size(ep_onset{m},2));
            end
            if ~isempty(ep_pull{m})
                bl = mean(ep_pull{m}(bl_pu, :), 1);
                ep_pull{m} = ep_pull{m} - bl;
                emg_pull_subj{m,s} = mean(ep_pull{m}, 2)';
            end
        end
    end

end  % subject loop

%% ===== GRAND AVERAGE =====
n_t_on = numel(t_onset);
n_t_pu = numel(t_pull);
ga_onset  = zeros(n_muscles, n_t_on);
ga_pull   = zeros(n_muscles, n_t_pu);
sem_onset = zeros(n_muscles, n_t_on);
sem_pull  = zeros(n_muscles, n_t_pu);

for m = 1:n_muscles
    valid_on = ~cellfun(@isempty, emg_onset_subj(m,:));
    if any(valid_on)
        mat = cell2mat(emg_onset_subj(m, valid_on)');
        ga_onset(m,:)  = mean(mat, 1);
        sem_onset(m,:) = std(mat, 0, 1) / sqrt(sum(valid_on));
    end
    valid_pu = ~cellfun(@isempty, emg_pull_subj(m,:));
    if any(valid_pu)
        mat = cell2mat(emg_pull_subj(m, valid_pu)');
        ga_pull(m,:)  = mean(mat, 1);
        sem_pull(m,:) = std(mat, 0, 1) / sqrt(sum(valid_pu));
    end
end

%% ===== FIGURE 1: Grand Average EMG Activation =====
shade = @(ax,t,mn,se,col) fill(ax, [t, fliplr(t)], [mn+se, fliplr(mn-se)], ...
    col, 'EdgeColor','none','FaceAlpha',COL_SHADE,'HandleVisibility','off');

figure('Name','Grand Average EMG Activation','Position',[50 50 1300 500],'Visible','off');

ax1 = subplot(1,2,1); hold on; grid on;
for m = 1:n_muscles
    shade(ax1, t_onset, ga_onset(m,:), sem_onset(m,:), MUSCLE_COLS(m,:));
    plot(ax1, t_onset, ga_onset(m,:), 'Color', MUSCLE_COLS(m,:), ...
        'LineWidth', 2, 'DisplayName', MUSCLE_LABELS{m});
end
xline(0, '--k', 'Mov onset (A)', 'LabelVerticalAlignment','bottom','LineWidth',1.5);
xlabel('Time (s)','FontSize',12);
ylabel('\DeltaNorm. amplitude (re baseline)','FontSize',12);
title(sprintf('Movement Onset (n=%d)', n_subj),'FontSize',13);
legend('Location','best','FontSize',10);
set(ax1,'FontSize',11); xlim([t_onset(1) t_onset(end)]);

ax2 = subplot(1,2,2); hold on; grid on;
for m = 1:n_muscles
    shade(ax2, t_pull, ga_pull(m,:), sem_pull(m,:), MUSCLE_COLS(m,:));
    plot(ax2, t_pull, ga_pull(m,:), 'Color', MUSCLE_COLS(m,:), ...
        'LineWidth', 2, 'DisplayName', MUSCLE_LABELS{m});
end
xline(0, '--k', 'Pull onset (D)', 'LabelVerticalAlignment','bottom','LineWidth',1.5);
xlabel('Time (s)','FontSize',12);
ylabel('\DeltaNorm. amplitude (re baseline)','FontSize',12);
title(sprintf('Pull Onset (n=%d)', n_subj),'FontSize',13);
legend('Location','best','FontSize',10);
set(ax2,'FontSize',11); xlim([t_pull(1) t_pull(end)]);

sgtitle('Grand Average EMG Activation — Contralateral Muscles  (IOD / Triceps / Deltoid)', ...
    'FontSize', 14);
set(gcf,'Units','pixels','Position',[50 50 1300 500]);
print(fullfile(RESULTS_DIR,'EMG_grandAvg_activation'), '-dpng', '-r300');
savefig(fullfile(RESULTS_DIR,'EMG_grandAvg_activation.fig'));

%% ===== FIGURES 2–3: Per-Subject grids per muscle =====
for m = 1:n_muscles
    % --- Movement Onset ---
    figure('Name', sprintf('Per-Subject EMG — %s — MovOnset', MUSCLE_LABELS{m}), ...
        'Position',[50 50 1400 700],'Visible','off');
    for s = 1:n_subj
        subplot(2,4,s); hold on; grid on;
        if isempty(emg_onset_subj{m,s})
            title([SUBJECTS{s} ' — MISSING'],'FontSize',9); continue;
        end
        plot(t_onset, emg_onset_subj{m,s}, 'Color', MUSCLE_COLS(m,:), 'LineWidth', 1.5);
        xline(0,'--k','LineWidth',1.2);
        title(SUBJECTS{s},'FontSize',11);
        xlabel('Time (s)','FontSize',9); ylabel('\DeltaAmp','FontSize',9);
        set(gca,'FontSize',9); xlim([t_onset(1) t_onset(end)]);
    end
    sgtitle(sprintf('Per-Subject EMG — %s — Movement Onset', MUSCLE_LABELS{m}),'FontSize',13);
    set(gcf,'Units','pixels','Position',[50 50 1400 700]);
    print(fullfile(RESULTS_DIR, sprintf('EMG_perSubject_MovOnset_%s', MUSCLE_LABELS{m})), ...
        '-dpng', '-r300');
    close(gcf);

    % --- Pull Onset ---
    figure('Name', sprintf('Per-Subject EMG — %s — Pull', MUSCLE_LABELS{m}), ...
        'Position',[50 50 1400 700],'Visible','off');
    for s = 1:n_subj
        subplot(2,4,s); hold on; grid on;
        if isempty(emg_pull_subj{m,s})
            title([SUBJECTS{s} ' — MISSING'],'FontSize',9); continue;
        end
        plot(t_pull, emg_pull_subj{m,s}, 'Color', MUSCLE_COLS(m,:), 'LineWidth', 1.5);
        xline(0,'--k','LineWidth',1.2);
        title(SUBJECTS{s},'FontSize',11);
        xlabel('Time (s)','FontSize',9); ylabel('\DeltaAmp','FontSize',9);
        set(gca,'FontSize',9); xlim([t_pull(1) t_pull(end)]);
    end
    sgtitle(sprintf('Per-Subject EMG — %s — Pull Onset', MUSCLE_LABELS{m}),'FontSize',13);
    set(gcf,'Units','pixels','Position',[50 50 1400 700]);
    print(fullfile(RESULTS_DIR, sprintf('EMG_perSubject_Pull_%s', MUSCLE_LABELS{m})), ...
        '-dpng', '-r300');
    close(gcf);
end

%% ===== SAVE =====
EMG_ACT.subjects        = SUBJECTS;
EMG_ACT.muscles         = MUSCLE_LABELS;
EMG_ACT.t_onset         = t_onset;
EMG_ACT.t_pull          = t_pull;
EMG_ACT.ga_onset        = ga_onset;   % [n_muscles × n_samp]
EMG_ACT.ga_pull         = ga_pull;
EMG_ACT.sem_onset       = sem_onset;
EMG_ACT.sem_pull        = sem_pull;
EMG_ACT.per_subj_onset  = emg_onset_subj;   % {n_muscles × n_subj} cell
EMG_ACT.per_subj_pull   = emg_pull_subj;
EMG_ACT.params.baseline_win_onset = BASELINE_WIN_ONSET;
EMG_ACT.params.baseline_win_pull  = BASELINE_WIN_PULL;
EMG_ACT.params.bp_band    = [BP_LOW, BP_HIGH];
EMG_ACT.params.env_cutoff = ENV_CUTOFF;
EMG_ACT.params.norm_pctile = NORM_PCTILE;

save(fullfile(RESULTS_DIR, 'EMG_activation.mat'), 'EMG_ACT', '-v7.3');
fprintf('\nSaved EMG activation results to %s\n', RESULTS_DIR);
