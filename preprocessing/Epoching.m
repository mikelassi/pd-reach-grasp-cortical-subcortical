%% Epoching.m
% =========================================================================
% Fixed-length epoching of unsegmented LFP data around two events of
% interest: movement onset (A) and pull onset (D).
%
% Epoch windows (from PreEpoching.m — 5th percentile criterion):
%   Movement onset  (locked to A):  [-0.5, +1.0] s
%     [-0.5, 0)  = pre-movement rest  → baseline for TF normalization
%     [0,  +1.0] = early reach        → covers 95% of reach durations
%
%   Pull onset  (locked to D):  [-0.3, +0.875] s
%     [-0.3, 0)  = late grasp context
%     [0, +0.875]= pull phase          → covers 95% of pull durations
%
% Outputs (saved to RESULTS_DIR/Epochs/):
%   MovOnset/  — per-subject EEGLAB .set files, all channels
%   Pull/      — per-subject EEGLAB .set files, all channels
%   Epochs_allSubjects.mat — consolidated struct with contralateral-channel
%                             data matrices + time vectors + metadata
%
% Author: Michael Lassi
% =========================================================================

close all; clear; clc;

%% ===== PARAMETERS =====

SUBJECTS    = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
BASE_PATH   = 'H:\Parkinson_ReachGrasp\Reprocessing';
PREPROC_SUB = fullfile('Preprocessed', 'LFP');
RESULTS_DIR = fullfile(BASE_PATH, 'RESULTS_final', 'Epochs');

% Epoch windows [s] — based on PreEpoching.m p5 results
PRE_ONSET_S  = 0.5;     % pre-A  rest baseline
POST_ONSET_S = 1.0;     % post-A reach
PRE_PULL_S   = 0.3;     % pre-D  grasp context
POST_PULL_S  = 0.875;   % post-D pull phase

% Contralateral hemisphere channel index (same convention as TF_analysis_LFP_v2.m)
CONTRA_CH = containers.Map({'wue02','wue03'}, {2, 2});
CH_DEFAULT = 1;

%% ===== INIT =====
eeglab nogui;

out_dir_onset = fullfile(RESULTS_DIR, 'MovOnset');
out_dir_pull  = fullfile(RESULTS_DIR, 'Pull');
if ~exist(out_dir_onset, 'dir'), mkdir(out_dir_onset); end
if ~exist(out_dir_pull,  'dir'), mkdir(out_dir_pull);  end

% Master struct for consolidated .mat output
EPOCHS = struct();

%% ===== SUBJECT LOOP =====
for s = 1:numel(SUBJECTS)
    subj         = SUBJECTS{s};
    preproc_path = fullfile(BASE_PATH, subj, PREPROC_SUB);
    set_files    = dir(fullfile(preproc_path, '*_wEv.set'));

    if isempty(set_files)
        warning('No *_wEv.set files found for %s — skipping.', subj);
        continue;
    end

    fprintf('\n===== Subject: %s  (%d block(s)) =====\n', subj, numel(set_files));

    % Contralateral channel for this subject
    if isKey(CONTRA_CH, subj)
        ch = CONTRA_CH(subj);
    else
        ch = CH_DEFAULT;
    end

    % Accumulate epochs across blocks before saving
    EEG_onset_all = [];
    EEG_pull_all  = [];
    trial_meta_onset = struct('subject', {}, 'block', {}, 'trial_in_block', {});
    trial_meta_pull  = struct('subject', {}, 'block', {}, 'trial_in_block', {});

    for f = 1:numel(set_files)
        EEG = pop_loadset('filename', set_files(f).name, 'filepath', preproc_path);
        fs  = EEG.srate;

        ev_types = {EEG.event.type};

        % ---- Build event-type lists matching A_T* and D_T* ----
        is_A = ~cellfun(@isempty, regexp(ev_types, '^A_T\d+$'));
        is_D = ~cellfun(@isempty, regexp(ev_types, '^D_T\d+$'));
        types_A = unique(ev_types(is_A));
        types_D = unique(ev_types(is_D));

        n_A = sum(is_A);
        n_D = sum(is_D);

        if isempty(types_A)
            warning('  Block %d: no A events found — skipping.', f);
            continue;
        end

        % ---- Epoch around A (movement onset) ----
        EEG_onset = pop_epoch(EEG, types_A, [-PRE_ONSET_S, POST_ONSET_S], ...
            'epochinfo', 'yes');
        EEG_onset = pop_rmbase(EEG_onset, []);  % do NOT apply baseline correction —
                                                 % leave raw power for TF normalisation
        n_onset_kept = EEG_onset.trials;
        fprintf('  Block %d | MovOnset: %d / %d epochs kept\n', f, n_onset_kept, n_A);

        % ---- Epoch around D (pull onset) ----
        if ~isempty(types_D)
            EEG_pull = pop_epoch(EEG, types_D, [-PRE_PULL_S, POST_PULL_S], ...
                'epochinfo', 'yes');
            EEG_pull = pop_rmbase(EEG_pull, []);
            n_pull_kept = EEG_pull.trials;
            fprintf('  Block %d | Pull    : %d / %d epochs kept\n', f, n_pull_kept, n_D);
        else
            warning('  Block %d: no D events found — pull epochs skipped.', f);
            EEG_pull = [];
            n_pull_kept = 0;
        end

        % ---- Concatenate across blocks ----
        if isempty(EEG_onset_all)
            EEG_onset_all = EEG_onset;
        else
            EEG_onset_all = pop_mergeset(EEG_onset_all, EEG_onset);
        end

        if ~isempty(EEG_pull)
            if isempty(EEG_pull_all)
                EEG_pull_all = EEG_pull;
            else
                EEG_pull_all = pop_mergeset(EEG_pull_all, EEG_pull);
            end
        end

        % ---- Store metadata ----
        for k = 1:n_onset_kept
            trial_meta_onset(end+1) = struct('subject', subj, ...
                'block', f, 'trial_in_block', k); %#ok<AGROW>
        end
        for k = 1:n_pull_kept
            trial_meta_pull(end+1) = struct('subject', subj, ...
                'block', f, 'trial_in_block', k); %#ok<AGROW>
        end
    end  % block loop

    if isempty(EEG_onset_all)
        warning('Subject %s: no onset epochs collected — skipping save.', subj);
        continue;
    end

    % ---- Save per-subject epoched .set files (all channels) ----
    fname_onset = sprintf('epochs_movOnset_%s.set', subj);
    fname_pull  = sprintf('epochs_pull_%s.set',     subj);

    pop_saveset(EEG_onset_all, 'filename', fname_onset, 'filepath', out_dir_onset);

    if ~isempty(EEG_pull_all)
        pop_saveset(EEG_pull_all, 'filename', fname_pull, 'filepath', out_dir_pull);
    end

    fprintf('  Saved: %s  (%d total epochs)\n', fname_onset, EEG_onset_all.trials);
    if ~isempty(EEG_pull_all)
        fprintf('  Saved: %s  (%d total epochs)\n', fname_pull, EEG_pull_all.trials);
    end

    % ---- Store contralateral-channel matrices in EPOCHS struct ----
    % Shape: [n_samples × n_trials]  (transposed for readability; TF scripts expect
    % [n_trials × n_samples] or use squeeze as needed)
    EPOCHS.(subj).fs           = fs;
    EPOCHS.(subj).ch_contra    = ch;

    EPOCHS.(subj).onset.data   = squeeze(EEG_onset_all.data(ch, :, :));  % [samp × trials]
    EPOCHS.(subj).onset.times  = EEG_onset_all.times / 1000;             % ms → s
    EPOCHS.(subj).onset.n      = EEG_onset_all.trials;
    EPOCHS.(subj).onset.meta   = trial_meta_onset;

    if ~isempty(EEG_pull_all)
        EPOCHS.(subj).pull.data  = squeeze(EEG_pull_all.data(ch, :, :));
        EPOCHS.(subj).pull.times = EEG_pull_all.times / 1000;
        EPOCHS.(subj).pull.n     = EEG_pull_all.trials;
        EPOCHS.(subj).pull.meta  = trial_meta_pull;
    end

end  % subject loop

%% ===== SAVE CONSOLIDATED .mat =====
EPOCHS.params.subjects     = SUBJECTS;
EPOCHS.params.pre_onset_s  = PRE_ONSET_S;
EPOCHS.params.post_onset_s = POST_ONSET_S;
EPOCHS.params.pre_pull_s   = PRE_PULL_S;
EPOCHS.params.post_pull_s  = POST_PULL_S;
EPOCHS.params.description  = ...
    'Contralateral-channel LFP epochs. onset: [-0.5,+1.0]s around A. pull: [-0.3,+0.875]s around D.';

if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end
save(fullfile(RESULTS_DIR, 'Epochs_allSubjects.mat'), 'EPOCHS', '-v7.3');
fprintf('\nSaved consolidated mat: %s\n', fullfile(RESULTS_DIR, 'Epochs_allSubjects.mat'));

%% ===== SUMMARY TABLE =====
fprintf('\n');
fprintf('===========================================================\n');
fprintf('  EPOCHING SUMMARY\n');
fprintf('===========================================================\n');
fprintf('%-8s  %10s  %10s\n', 'Subject', 'MovOnset N', 'Pull N');
fprintf('%s\n', repmat('-', 1, 34));
total_onset = 0; total_pull = 0;
for s = 1:numel(SUBJECTS)
    subj = SUBJECTS{s};
    if ~isfield(EPOCHS, subj), fprintf('%-8s  %10s  %10s\n', subj, 'SKIPPED', 'SKIPPED'); continue; end
    n_on = EPOCHS.(subj).onset.n;
    n_pu = 0;
    if isfield(EPOCHS.(subj), 'pull'), n_pu = EPOCHS.(subj).pull.n; end
    fprintf('%-8s  %10d  %10d\n', subj, n_on, n_pu);
    total_onset = total_onset + n_on;
    total_pull  = total_pull  + n_pu;
end
fprintf('%s\n', repmat('-', 1, 34));
fprintf('%-8s  %10d  %10d\n', 'TOTAL', total_onset, total_pull);
fprintf('\nEpoch windows:\n');
fprintf('  MovOnset : [%.1f, +%.1f] s  (%d samples at %d Hz)\n', ...
    -PRE_ONSET_S, POST_ONSET_S, round((PRE_ONSET_S + POST_ONSET_S) * 400), 400);
fprintf('  Pull     : [%.1f, +%.3f] s  (%d samples at %d Hz)\n', ...
    -PRE_PULL_S, POST_PULL_S, round((PRE_PULL_S + POST_PULL_S) * 400), 400);
