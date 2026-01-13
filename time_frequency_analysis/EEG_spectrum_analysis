%% =======================================================
% EEG Time–Frequency Analysis and Trial-wise Normalization
%
% This script performs EEG power spectral analysis for multiple subjects
% across behavioral phases (Rest, Reach, Grasp, Pull) and channels (126).
% The workflow includes:
%   1. Loading preprocessed EEG data and trial-specific segments.
%   2. Computing time-frequency decomposition using Continuous Wavelet Transform (Morse wavelets).
%   3. Calculating trial-wise power spectral density (PSD) and baseline-normalized PSD
%      using pre- and post-Rest periods.
%   4. Aggregating PSDs across trials and blocks to obtain subject-level median PSD.
%   5. Organizing data into structures suitable for topographic visualization
%      and connectivity analyses in MATLAB and Python.
%
% Outputs:
%   - EEG_PSD_ALL: Subject-level PSD (raw and normalized), median across trials.
%   - EEG_PSD_ALL_CONNECTIVITY: Trial-wise PSD arrays for connectivity analyses,
%     median PSD across trials, compatible with HDF5 format.
%
% This pipeline enables detailed examination of band-specific EEG dynamics
% (High/Low Beta, Motor-Related Beta) during motor tasks and supports
% subsequent topoplot and connectivity analyses.
%% =======================================================

clear; clc; close all;

subject_list = {'wue02', 'wue03', 'wue05', 'wue06', 'wue07', 'wue09', 'wue10', 'wue11'}; 
base_path    = 'C:\Users\tomma\OneDrive - University of Pisa\Desktop\TESI\Dataset_tesi';
preproc_subfolder = 'Preprocessed_01\EEG_manual';

n_subjects = length(subject_list);
n_channels = 126;
eeglab;

phase_names = {'Rest', 'Reach', 'Grasp', 'Pull'};
n_phases = length(phase_names);

% Store PSD norm and not normalized
all_subj_psd = cell(n_phases,n_channels, n_subjects);
all_subj_psd_norm = cell(n_phases,n_channels, n_subjects);

% Store Kin and freq for plot across subjects
all_subj_fr = cell(1, n_subjects);

% Struttura finale 
EEG_PSD_ALL = struct();

EEG_PSD_STRUCT = {};

for s = 1:length(subject_list)
    
    subject_id = subject_list{s};
    preproc_path = fullfile(base_path, subject_id, preproc_subfolder);
    
    % Load EEG data (.set)
    set_files = dir(fullfile(preproc_path, '*.set'));
    mat_files = dir(fullfile(preproc_path, '*_EEG_trialsByRegionAndPhase.mat'));
    if isempty(set_files) || isempty(mat_files)
        warning('No EEG files found for %s', subject_id);
        continue;
    end

    n_blocks = length(mat_files);
    
    % Store for global frequency axis 
    all_subj_freq = cell(1, n_blocks);
    all_blocks_freq_length = zeros(1,n_blocks);
    
    all_blocks_phase_psd = cell(n_phases, n_channels, n_blocks);
    all_blocks_phase_psd_norm = cell(n_phases, n_channels, n_blocks);

    % --- Loop sui blocchi --- %
    for f = 1:length(mat_files)
    
        set_file = set_files(f).name;
        mat_file = mat_files(f).name;

        EEG = pop_loadset('filename', set_file, 'filepath', preproc_path);
        E = load(fullfile(preproc_path, mat_file));
        fs = EEG.srate;
        
        n_channels = length(EEG.chanlocs);
        all_ch_idxs = 1:length(EEG.chanlocs);
        all_labels.(subject_id) = {EEG.chanlocs.labels};
    
        % === Define behavioral phases ===
        EEG_phases = E.EEG_phases;
        EEG_trials = E.EEG_trials;
        baseline_EEG = E.baseline_EEG;
        
        [num_trials, num_phases] = size(EEG_phases);

        % --- Loop sui canali ---
        for ch = 1:n_channels
            
            ch_idx = ch;
            
            % Baseline CWT 
            baseline_data = baseline_EEG(ch_idx, :);
            baseline_length = length(baseline_data);
            baseline_idxs = 1:baseline_length;

            % Global CWT computation 
            EEG_channel_data = EEG.data(ch_idx,:);
            [wt_glob, fr_glob] = cwt(EEG_channel_data, fs, 'morse');
            power_wt_glob = abs(wt_glob).^2;
            
            if ch == 1 
                all_subj_freq{f} = fr_glob;
                n_fr = size(fr_glob,1);
                all_blocks_freq_length(1,f) = n_fr;
            end

            phase_psd_all_trials = zeros(4, n_fr, num_trials);
            phase_psd_norm_all_trials = zeros(4, n_fr, num_trials);

            for t = 1:num_trials
                trial_data = EEG_trials{t}(ch_idx,:);

                trial_length = length(trial_data);
    
                n_eeg_phases = size(EEG_phases,2);

                for p = 1:n_eeg_phases

                    phase_data = EEG_phases{t,p}(ch_idx,:);
                    phase_length = length(phase_data);
    
                    % Estrazione regione trial (dall'intero segnale preprocessato
                    % l'inizio del primo trial si trova 250 campioni prima
                    % dell'evento A)
                    if t == 1 && p == 1
                        start_idx =  baseline_length - 250;
                    else
                        start_idx = end_idx + 1;
                    end
                    end_idx = start_idx + length(phase_data) - 1;
                    
                    % Salvo anche regione della fase corrente -->
                    % potenza non normalizzata
                    pow_wt_phase_no_norm = power_wt_glob(:, start_idx:end_idx);

                    % --- COSTRUZIONE BASELINE (REST PRE + POST)
                    if p == 1
                        pow_rest_bef = pow_wt_phase_no_norm;
                        length_rest_aft = length(EEG_phases{t,5}(ch_idx,:));
                        idx_rest_aft = (start_idx + trial_length - length_rest_aft): (start_idx + trial_length - 1);
                        pow_rest_aft = power_wt_glob(:,idx_rest_aft);
                        pow_rest = cat(2, pow_rest_bef, pow_rest_aft);
                
                        % PSD media del rest (baseline del trial)
                        mean_baseline_trial = mean(pow_rest, 2);

                        % Salvataggio PSD non normalizzata del rest
                        phase_psd_all_trials(1,:,t) = mean(pow_rest, 2);
                        
                        % Salvataggio PSD normalizzata del rest 
                        pow_rest_norm = (pow_rest - mean_baseline_trial)./mean_baseline_trial;
                        phase_psd_norm_all_trials(p,:,t) = mean(pow_rest_norm,2);
    
                    elseif p > 1 && p < 5
                        phase_psd_all_trials(p,:,t) = mean(pow_wt_phase_no_norm,2);

                        % Salvataggio PSD normalizzata nelle fasi
                        % cinematiche
                        pow_wt_phase = (pow_wt_phase_no_norm - mean_baseline_trial)./mean_baseline_trial;
                        phase_psd_norm_all_trials(p,:,t) = mean(pow_wt_phase,2);
                    end

                    % Normalization by trial baseline for all phases
                    pow_wt_phase = (pow_wt_phase_no_norm - mean_baseline_trial)./mean_baseline_trial;
                end
            end

            % === Salvataggio PSD per-trial ===
            subj_field = subject_id;
            for t = 1:num_trials
                for p = 1:n_phases
                    EEG_PSD_STRUCT.(subj_field).block(f).channel(ch).trial(t).phase(p).psd_norm = squeeze(phase_psd_norm_all_trials(p,:,t));
                    EEG_PSD_STRUCT.(subj_field).block(f).channel(ch).trial(t).phase(p).psd = squeeze(phase_psd_all_trials(p,:,t));
                    EEG_PSD_STRUCT.(subject_id).block(f).frequencies = fr_glob;
                end
            end
            
            % Average psd and baseline 
            for p = 1:n_phases
                all_blocks_phase_psd{p,ch,f} = median(phase_psd_all_trials(p,:,:),3); % Median across trials
                all_blocks_phase_psd_norm{p,ch,f} = median(phase_psd_norm_all_trials(p,:,:),3);
            end
        end
    end

   
    % Trovo Lunghezza minima fr_glob per soggetto
    min_fr_glob = min(all_blocks_freq_length);
    min_fr_glob_mask = 1:min_fr_glob;

    % Fr glob tagliato a min_length
    subj_fr = all_subj_freq{1}(1:min_fr_glob);
    
    % Matrices for unfolding PSD cell arrays
    all_psd_subj = zeros(n_phases,min_fr_glob, n_blocks, n_channels);
    all_psd_norm_subj = zeros(n_phases, min_fr_glob, n_blocks, n_channels);

    for ch = 1:n_channels
        ch_name = all_labels.(subject_id){ch};
    
        for p = 1:n_phases
            all_trials_psd = zeros(min_fr_glob, 0);
            all_trials_psd_norm = zeros(min_fr_glob, 0);
    
            % Concateno tutti i trials di tutti i blocchi per la fase p
            for k = 1:n_blocks
                n_trials_block = numel(EEG_PSD_STRUCT.(subject_id).block(k).channel(ch).trial);
                for t = 1:n_trials_block
                    psd_vec = EEG_PSD_STRUCT.(subject_id).block(k).channel(ch).trial(t).phase(p).psd(1:min_fr_glob);
                    psd_norm_vec = EEG_PSD_STRUCT.(subject_id).block(k).channel(ch).trial(t).phase(p).psd_norm(1:min_fr_glob);
    
                    % force column and concatenate as new column
                    all_trials_psd(:, end+1) = psd_vec(:);
                    all_trials_psd_norm(:, end+1) = psd_norm_vec(:);
                end
            end

            EEG_PSD_ALL(s).PSD_all{p} = all_trials_psd;
            EEG_PSD_ALL(s).PSD_norm_all{p} = all_trials_psd_norm;
    
            % Mediana tra tutti i trials della fase p
            all_subj_psd{p,ch,s} = median(all_trials_psd, 2);
            all_subj_psd_norm{p,ch,s} = median(all_trials_psd_norm, 2);
        end
    end

    EEG_PSD_ALL(s).subject_id = subject_id;
    EEG_PSD_ALL(s).freq_axis = subj_fr;
    
    EEG_PSD_ALL(s).PSD = all_subj_psd(:,:,s);           % [phases x freq x channels]
    EEG_PSD_ALL(s).PSD_norm = all_subj_psd_norm(:,:,s); % [phases x freq x channels]
end

data_file = fullfile(base_path, 'RESULTS', 'Topoplot_NORM_median_01', 'Mat_data', 'All_Data_TF.mat');
save(data_file, 'EEG_PSD_ALL', '-v7.3');


%% Costruzione struttura finale per import in Python
EEG_PSD_ALL_CONNECTIVITY = struct();

for s = 1:length(subject_list)
    
    subject_id = subject_list{s};
    
    % Trova numero totale di trials (dal primo canale)
    n_blocks = length(EEG_PSD_STRUCT.(subject_id).block);
    total_trials = 0;
    for f = 1:n_blocks
        total_trials = total_trials + numel(EEG_PSD_STRUCT.(subject_id).block(f).channel(1).trial);
    end
    
    n_phases = length(phase_names);
    n_channels = 126;
    
    % Trova numero minimo di frequenze tra tutti i blocchi
    min_n_freqs = inf;
    for f = 1:n_blocks
        n_freq_curr = length(EEG_PSD_STRUCT.(subject_id).block(f).frequencies);
        if n_freq_curr < min_n_freqs
            min_n_freqs = n_freq_curr;
        end
    end
    
    PSD_all_array = zeros(n_phases, min_n_freqs, n_channels, total_trials);
    PSD_norm_all_array = zeros(n_phases, min_n_freqs, n_channels, total_trials);
    
    trial_idx = 1;
    for f = 1:n_blocks
        n_trials_block = numel(EEG_PSD_STRUCT.(subject_id).block(f).channel(1).trial);
        for t = 1:n_trials_block
            for ch = 1:n_channels
                for p = 1:n_phases
                    PSD_all_array(p, :, ch, trial_idx) = ...
                        EEG_PSD_STRUCT.(subject_id).block(f).channel(ch).trial(t).phase(p).psd(1:min_n_freqs);
                    PSD_norm_all_array(p, :, ch, trial_idx) = ...
                        EEG_PSD_STRUCT.(subject_id).block(f).channel(ch).trial(t).phase(p).psd_norm(1:min_n_freqs);
                end
            end
            trial_idx = trial_idx + 1;   % Incrementa l'indice globale
        end
    end
    
    % Salva nella struttura finale
    EEG_PSD_ALL_CONNECTIVITY(s).subject_id = subject_id;
    EEG_PSD_ALL_CONNECTIVITY(s).freq_axis = EEG_PSD_STRUCT.(subject_id).block(1).frequencies(1:min_n_freqs);
    
    % Mediana su tutti i trials
    EEG_PSD_ALL_CONNECTIVITY(s).PSD = median(PSD_all_array, 4);        % [phases x freq x channels]
    EEG_PSD_ALL_CONNECTIVITY(s).PSD_norm = median(PSD_norm_all_array, 4);
    
    EEG_PSD_ALL_CONNECTIVITY(s).PSD_all = PSD_all_array;               % [phases x freq x channels x trials]
    EEG_PSD_ALL_CONNECTIVITY(s).PSD_norm_all = PSD_norm_all_array;
    
end

% Salva come file v7.3 compatibile HDF5
save(fullfile(base_path, 'RESULTS', 'Topoplot_NORM_median_01', 'Mat_data', 'All_Data_TF_CONNECTIVITY.mat'), ...
     'EEG_PSD_ALL_CONNECTIVITY', '-v7.3');

disp('Struttura EEG_PSD_ALL_CONNECTIVITY salvata correttamente.');

