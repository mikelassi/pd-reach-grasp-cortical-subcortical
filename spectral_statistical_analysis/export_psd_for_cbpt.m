%% export_psd_for_cbpt.m
% =========================================================================
% Flatten PSD_STRUCT_allSubjects.mat (saved -v7.3 by TF_analysis_LFP_v2.m)
% into a tidy, scipy-readable (-v7) file for the subject-level CBPT.
%
% Aggregation: for each subject and kinematic phase, take the MEDIAN over
% all of that subject's trials (pooled across blocks) of the per-trial,
% rest-normalised phase PSD (dB).  This collapses the trial dimension so the
% statistical unit becomes the SUBJECT (n=8) — the correct level of
% inference for a group claim (fixes the pseudoreplication in the old CBPT).
%
% Output (cbpt_input.mat, -v7):
%   freqs            [n_fr x 1]      shared frequency axis (Hz)
%   subjects         {1 x n_subj}    subject ids in order
%   phase_names      {1 x n_phase}   {'Rest_pre','Reach','Grasp','Pull','Rest_post'}
%   psd_subj_phase   [n_subj x n_phase x n_fr]  median-over-trials, dB re rest
%   psd_subj_phase_mean [n_subj x n_phase x n_fr]  mean-over-trials (sensitivity)
%   n_trials_subj    [n_subj x 1]    pooled trial count per subject
%
% Author: Michael Lassi
% =========================================================================

clear; clc;

IN_MAT  = 'F:\Projects\Parkinson_ReachGrasp\Reprocessing\RESULTS_compare\v2\Baseline_01_median\PSD_STRUCT_allSubjects.mat';
OUT_MAT = fullfile(fileparts(mfilename('fullpath')), 'cbpt_input.mat');

phase_names = {'Rest_pre','Reach','Grasp','Pull','Rest_post'};
n_phase     = numel(phase_names);

fprintf('Loading %s ...\n', IN_MAT);
S = load(IN_MAT, 'PSD_STRUCT');
PSD_STRUCT = S.PSD_STRUCT;

subjects = fieldnames(PSD_STRUCT)';   % row cell of subject ids
n_subj   = numel(subjects);

% Shared frequency axis (v2 uses a fixed cwtfilterbank → identical for all)
freqs  = PSD_STRUCT.(subjects{1}).block(1).frequencies(:);
n_fr   = numel(freqs);

psd_subj_phase      = nan(n_subj, n_phase, n_fr);
psd_subj_phase_mean = nan(n_subj, n_phase, n_fr);
n_trials_subj       = zeros(n_subj, 1);

% Trial-level (long format) accumulators for the precision-weighted and
% mixed-effects CBPT variants — every trial kept, tagged with its subject.
trial_psd_cell      = {};   % each: [n_trials_s x n_phase x n_fr]
trial_subject_cell  = {};   % each: [n_trials_s x 1] integer subject index

for s = 1:n_subj
    subj   = subjects{s};
    blocks = PSD_STRUCT.(subj).block;

    % Pool per-phase trial PSDs across all blocks → [n_fr x n_trials_total]
    pool = cell(1, n_phase);
    for p = 1:n_phase, pool{p} = zeros(n_fr, 0); end

    for b = 1:numel(blocks)
        if isempty(blocks(b).trial), continue; end
        trials = blocks(b).trial;
        for t = 1:numel(trials)
            % Guard against any frequency-length drift (shouldn't happen in v2)
            for p = 1:n_phase
                v = trials(t).phase(p).psd_norm(:);
                if numel(v) ~= n_fr
                    v = interp1(linspace(0,1,numel(v))', v, linspace(0,1,n_fr)', 'linear');
                end
                pool{p}(:, end+1) = v; %#ok<AGROW>
            end
        end
    end

    n_trials_subj(s) = size(pool{1}, 2);
    for p = 1:n_phase
        psd_subj_phase(s, p, :)      = median(pool{p}, 2);
        psd_subj_phase_mean(s, p, :) = mean(pool{p},  2);
    end

    % Pack this subject's trials into [n_trials_s x n_phase x n_fr]
    nt = n_trials_subj(s);
    subj_trials = zeros(nt, n_phase, n_fr);
    for p = 1:n_phase
        subj_trials(:, p, :) = pool{p}';   % [n_trials_s x n_fr]
    end
    trial_psd_cell{end+1}     = subj_trials;        %#ok<AGROW>
    trial_subject_cell{end+1} = repmat(s, nt, 1);   %#ok<AGROW>

    fprintf('  %-7s : %d trials pooled\n', subj, n_trials_subj(s));
end

% Concatenate trial-level long-format arrays across subjects
trial_psd     = cat(1, trial_psd_cell{:});       % [n_total x n_phase x n_fr]
trial_subject = cat(1, trial_subject_cell{:});   % [n_total x 1]

save(OUT_MAT, 'freqs', 'subjects', 'phase_names', ...
    'psd_subj_phase', 'psd_subj_phase_mean', 'n_trials_subj', ...
    'trial_psd', 'trial_subject', '-v7');
fprintf('\nSaved %s\n', OUT_MAT);
fprintf('  psd_subj_phase: [%d subj x %d phase x %d freq]\n', n_subj, n_phase, n_fr);
fprintf('  trial_psd:      [%d trials x %d phase x %d freq]\n', size(trial_psd,1), n_phase, n_fr);
fprintf('  freq range: %.2f - %.2f Hz\n', freqs(end), freqs(1));
