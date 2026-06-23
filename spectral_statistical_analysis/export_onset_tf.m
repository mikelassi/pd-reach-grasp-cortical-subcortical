%% export_onset_tf.m
% Flatten the onset-locked TF (from TF_analysis_epoched.m, -v7.3) into a
% scipy-readable -v7 file for the event-locked ERD cluster test.
clear; clc;

IN  = 'F:\Projects\Parkinson_ReachGrasp\Reprocessing\RESULTS_final\TF_epoched\TF_allSubjects.mat';
OUT = fullfile(fileparts(mfilename('fullpath')), 'onset_tf.mat');

S = load(IN, 'TF'); TF = S.TF;

per_subj_onset = TF.per_subj_onset;     % [n_fr x n_samp x n_subj], dB re pre-move rest
t_onset        = TF.t_onset(:);         % seconds, 0 = movement onset
fr             = TF.fr(:);              % Hz (descending)
subjects       = TF.subjects;
if isfield(TF, 'coi_mask_onset')
    coi_mask_onset = TF.coi_mask_onset; % [n_fr x n_samp] logical (true = valid)
else
    coi_mask_onset = true(numel(fr), numel(t_onset));
end

fprintf('per_subj_onset: %s\n', mat2str(size(per_subj_onset)));
fprintf('t_onset: %d samples  [%.3f .. %.3f] s\n', numel(t_onset), t_onset(1), t_onset(end));
fprintf('fr: %d  [%.2f .. %.2f] Hz\n', numel(fr), fr(end), fr(1));

save(OUT, 'per_subj_onset', 't_onset', 'fr', 'subjects', 'coi_mask_onset', '-v7');
fprintf('Saved %s\n', OUT);
