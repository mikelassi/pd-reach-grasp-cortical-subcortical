%% EEG ANALYSIS -- Several ROI Power TF
% Time-frequency analysis with trial-wise interpolation and baseline normalization.
dbstop if error

clear; clc; close all;

subject_list = {'wue02', 'wue03', 'wue05', 'wue06', 'wue07', 'wue09', 'wue10', 'wue11'}; 
base_path    = 'C:\Users\tomma\OneDrive - University of Pisa\Desktop\TESI\Dataset_tesi';
%preproc_subfolder = 'Preprocessed_manual_EEG';
preproc_subfolder = 'Preprocessed_01\EEG_manual';

n_subjects = length(subject_list);


phase_names = {'Rest', 'Reach', 'Grasp', 'Pull'};
n_phases = length(phase_names);

% Definizioni regioni d'interesse (ROI) 

% Supplementary Motor Area (SMA, BA 6 medial)
ROI_SMA = {'FC1','FC2','FC3','FC4', 'FCz'}; 

% Primary Motor Cortex (Precentral gyrus, BA 4)
ROI_M1_left  = {'C3','C1','FCC3h','FCC1h', 'CCP1h', 'CCP3h'};
ROI_M1_right = {'C4','C2','FCC4h','FCC2h', 'CCP2h', 'CPP4h'};

ROI_M1 = [ROI_M1_left, ROI_M1_right, 'Cz'];

% Ventrolateral Prefrontal Cortex (BA 44/45/47)
ROI_VLPFC_left  = {'F7','F5','FT7','FC5','FFC5h','AF7'};
ROI_VLPFC_right = {'F8','F6','FT8','FC6','FFC6h','AF8'};

ROI_VLPFC = [ROI_VLPFC_left, ROI_VLPFC_right];

% Temporal Cortex (Superior Temporal Gyrus, BA 21/22)
ROI_Temporal_left  = {'T7','TP7','TTP7h'};
ROI_Temporal_right = {'T8','TP8','TTP8h'};

% Parietal Cortex (BA 7 - superior/inferior parietal lobule)
ROI_Parietal = {'P3','P1','P5','P7','P9' 'P4','P2','P6','P8','P10'};

% Occipital Cortex (BA 17-19)
ROI_Occipital = {'Oz','O1','O2','POz','POO1','POO2', 'PO1h','PO2h','PO3','PO4'};

% 26 channels
ROI_CENTRAL_SCALP = {'FC5','FC3','FC1','FC2','FC4','FC6','C5', 'C3','C1','C2','C4','C6', 'CP5','CP3','CP1','CP2','CP4','CP6','CCP3h','CCP1h','CCP2h','CCP4h','FCC3h','FCC1h','FCC2h','FCC4h'};

ROI_CENTRAL_left = {'FC5','FC3','FC1','C5', 'C3','C1', 'CP5','CP3','CP1','CCP3h','CCP1h','FCC3h','FCC1h'};

ROI_CENTRAL_right = {'FC2','FC4','FC6','C2','C4','C6', 'CP2','CP4','CP6','CCP2h','CCP4h','FCC2h','FCC4h'};

% Clusters identified in Topoplot Analysis --> 18 channels
% --- Left cluster
LEFT_CLUSTER_CHANNELS = {'C1','C3','CCP1h','CCP3h','CCP5h','CP1','CP3','CPP3h','FCC3h'};

% --- Right cluster
RIGHT_CLUSTER_CHANNELS = {'C4','CCP4h','CCP6h','CP2','CP4','CP6','CPP4h','CPP6h','P4'};

CLUSTER_CHANNELS = [LEFT_CLUSTER_CHANNELS, RIGHT_CLUSTER_CHANNELS];

%ROI_names = {'PSCContro', 'DLPFCContro', 'M1Contro','M1Ipsi','M1','TemporalContro','SMA','Occipital', 'VLPFC', 'Parietal'};
ROI_names = {'SMA','M1Contro','M1Ipsi','M1', 'VLPFC', 'Parietal', 'Central', 'CentralIpsi','CentralContro','Occipital', 'Cluster', 'ClusterIpsi', 'ClusterContro'};

n_ROI = length(ROI_names);

eeglab;

percentile_caxis = 0.95;

% Store Avg Power TF of Baseline for plot across subjects
all_ROI_med_subj_baseline = cell(n_ROI, n_subjects);
all_subj_wt = cell(n_ROI, n_subjects);
all_subj_wt_no_norm = cell(n_ROI, n_subjects);

% Store PSD phases
all_subj_psd = cell(n_phases, n_ROI, n_subjects);
all_subj_psd_norm = cell(n_phases, n_ROI, n_subjects);

all_subj_psd_trial = cell(n_ROI, n_subjects);
all_subj_psd_norm_trial = cell(n_ROI, n_subjects);

% Store Kin and freq for plot across subjects
all_subj_kin = cell(1, n_subjects);
all_subj_fr = cell(1, n_subjects);


% Store for kinematic latencies
all_kin_latencies = zeros(4,n_subjects);

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

    % % Handle corrupted block 
    % if any(strcmp(subject_id, {'wue06'})) 
    %     mat_files = mat_files(1:2);      % Third block is corrupted for wue06
    % end

    n_blocks = length(mat_files);

    % Store for global frequency axis 
    all_subj_freq = cell(1, n_blocks);
    all_blocks_freq_length = zeros(1,n_blocks);

    % Store for Baseline Power average across blocks in a ROI
    all_blocks_ROI_baseline = cell(n_ROI, n_blocks);
    
    % Store for ROI Power TF
    all_blocks_ROI_power_wt = cell(n_ROI, n_blocks);
    all_blocks_ROI_power_wt_no_norm = cell(n_ROI, n_blocks);
    all_blocks_kin = cell(n_blocks, 1);

    % Store for PSD mean across blocks
    all_blocks_ROI_psd = cell(n_phases, n_ROI, n_blocks);
    all_blocks_ROI_psd_norm = cell(n_phases, n_ROI, n_blocks);

    all_blocks_ROI_psd_trial = cell(n_ROI, n_blocks);
    all_blocks_ROI_psd_norm_trial = cell(n_ROI, n_blocks);

    % Store for kinematic event latencies across blocks
    onsets = [];
    pulls = [];
    grasps = [];
    offsets = [];
    
    % --- Loop sui blocchi --- %
    for f = 1:length(mat_files)
    
        set_file = set_files(f).name;
        mat_file = mat_files(f).name;

        EEG = pop_loadset('filename', set_file, 'filepath', preproc_path);
        E = load(fullfile(preproc_path, mat_file));
        fs = EEG.srate;

        % Select electrodes of ROI controlateral to hand performing movement

        if any(strcmp(subject_id, {'wue05', 'wue09'}))   % Left handed
            % [~, ch_idxs_PSC] = ismember(ROI_PSC_right, {EEG.chanlocs.labels});
            [~, ch_idxs_M1_contro] = ismember(ROI_M1_right, {EEG.chanlocs.labels});
            [~, ch_idxs_M1_ipsi] = ismember(ROI_M1_left, {EEG.chanlocs.labels});
            % [~, ch_idxs_DLPFC] = ismember(ROI_DLPFC_right, {EEG.chanlocs.labels});
            [~, ch_idxs_Central_contro] = ismember(ROI_CENTRAL_right, {EEG.chanlocs.labels});
            [~, ch_idxs_Central_ipsi] = ismember(ROI_CENTRAL_left, {EEG.chanlocs.labels});
            [~, ch_idxs_CLUSTER_IPSI] = ismember(LEFT_CLUSTER_CHANNELS, {EEG.chanlocs.labels});
            [~, ch_idxs_CLUSTER_CONTRO] = ismember(RIGHT_CLUSTER_CHANNELS, {EEG.chanlocs.labels});
        else
            % [~, ch_idxs_PSC] = ismember(ROI_PSC_left, {EEG.chanlocs.labels}); % Right handed
            [~, ch_idxs_M1_contro] = ismember(ROI_M1_left, {EEG.chanlocs.labels});
            [~, ch_idxs_M1_ipsi] = ismember(ROI_M1_right, {EEG.chanlocs.labels});
            % [~, ch_idxs_DLPFC] = ismember(ROI_DLPFC_left, {EEG.chanlocs.labels});
            [~, ch_idxs_Central_contro] = ismember(ROI_CENTRAL_left, {EEG.chanlocs.labels});
            [~, ch_idxs_Central_ipsi] = ismember(ROI_CENTRAL_right, {EEG.chanlocs.labels});
            [~, ch_idxs_CLUSTER_IPSI] = ismember(RIGHT_CLUSTER_CHANNELS, {EEG.chanlocs.labels});
            [~, ch_idxs_CLUSTER_CONTRO] = ismember(LEFT_CLUSTER_CHANNELS, {EEG.chanlocs.labels});
        end
        
        % Ch_idxs for entire ROIs (not only controlateral)
        [~, ch_idxs_M1] = ismember(ROI_M1, {EEG.chanlocs.labels});
        ch_idxs_M1 = ch_idxs_M1(ch_idxs_M1>0);
        [~, ch_idxs_SMA] = ismember(ROI_SMA, {EEG.chanlocs.labels});
        ch_idxs_SMA = ch_idxs_SMA(ch_idxs_SMA>0);
        [~, ch_idxs_Occipital] = ismember(ROI_Occipital, {EEG.chanlocs.labels});
        ch_idxs_Occipital = ch_idxs_Occipital(ch_idxs_Occipital>0);
        [~, ch_idxs_VLPFC] = ismember(ROI_VLPFC, {EEG.chanlocs.labels});
        ch_idxs_VLPFC = ch_idxs_VLPFC(ch_idxs_VLPFC>0);
        [~, ch_idxs_Parietal] = ismember(ROI_Parietal, {EEG.chanlocs.labels});
        ch_idxs_Parietal = ch_idxs_Parietal(ch_idxs_Parietal>0);
        [~, ch_idxs_Central] = ismember(ROI_CENTRAL_SCALP, {EEG.chanlocs.labels});
        ch_idxs_Central = ch_idxs_Central(ch_idxs_Central>0);
        [~, ch_idxs_CLUSTER] = ismember(CLUSTER_CHANNELS, {EEG.chanlocs.labels});
        ch_idxs_CLUSTER = ch_idxs_CLUSTER(ch_idxs_CLUSTER>0);

        % --- Loop sulle ROI ---%
        %ROIs_ch_idxs = {ch_idxs_PSC, ch_idxs_DLPFC, ch_idxs_M1_contro, ch_idxs_M1_ipsi, ch_idxs_M1, ch_idxs_Temporal_contro, ch_idxs_SMA, ch_idxs_Occipital, ch_idxs_VLPFC, ch_idxs_Parietal};
        ROIs_ch_idxs = {ch_idxs_SMA, ch_idxs_M1_contro, ch_idxs_M1_ipsi, ch_idxs_M1, ch_idxs_VLPFC, ch_idxs_Parietal, ch_idxs_Central, ch_idxs_Central_ipsi, ch_idxs_Central_contro, ch_idxs_Occipital, ch_idxs_CLUSTER, ch_idxs_CLUSTER_IPSI, ch_idxs_CLUSTER_CONTRO};

        % % Handle EEG corrupted trials for wue03 ( block 3 trial 4) and wue05
        % % (block 1 trial 4) 
        % 
        % if (any(strcmp(subject_id, {'wue03'})) && f == 3) || (any(strcmp(subject_id, {'wue05'})) && f == 1 )
        %     E.EEG_phases_vis = E.EEG_phases_vis([1:3,5:end],:);        % No trial 4
        %     E.EEG_phases = E.EEG_phases([1:3,5:end],:);        
        %     E.EEG_trials = E.EEG_trials([1:3,5:end],:);
        %     E.kinematic_4EEG_trials = E.kinematic_4EEG_trials([1:3,5:end],:);  
        % end

        all_labels.(subject_id) = {EEG.chanlocs.labels};
    
        % === Define behavioral phases ===
        EEG_phases = E.EEG_phases;
        KIN_trials = E.kinematic_4EEG_trials; 
        EEG_trials = E.EEG_trials;
        baseline_EEG = E.baseline_EEG;
        
        [num_trials, num_phases] = size(EEG_phases);

        % Compute median length of each phase across trials --> for
        % 1 electrode only
        len_phase1 = zeros(num_trials,1);
        len_phase2 = zeros(num_trials,1);
        len_phase3 = zeros(num_trials,1);
        len_phase4 = zeros(num_trials,1);
        len_phase5 = zeros(num_trials,1);
        
        for t = 1:num_trials
            len_phase1(t) = size(EEG_phases{t,1}, 2);
            len_phase2(t) = size(EEG_phases{t,2}, 2);
            len_phase3(t) = size(EEG_phases{t,3}, 2);
            len_phase4(t) = size(EEG_phases{t,4}, 2);
            len_phase5(t) = size(EEG_phases{t,5}, 2);
        end

        % Median values
        median_len_1 = round(median(len_phase1));
        median_len_2 = round(median(len_phase2));
        median_len_3 = round(median(len_phase3));
        median_len_4 = round(median(len_phase4));
        median_len_5 = round(median(len_phase5));
        median_lengths_phase = [median_len_1, median_len_2, median_len_3, median_len_4, median_len_5];
        tot_length = median_len_1 + median_len_2 + median_len_3 + median_len_4 + median_len_5;

        % Calcolo mediana latencies degli eventi cinematici
        onset = median_lengths_phase(1)/fs;
        grasp = onset + median_lengths_phase(2)/fs;
        pull = grasp + median_lengths_phase(3)/fs;
        offset = pull + median_lengths_phase(4)/fs;
        
        % One for block
        all_kin_norm_block = zeros(tot_length, num_trials);
        
        % --- Loop sulle ROI ---
        for r = 1:n_ROI
            ch_idxs = ROIs_ch_idxs{r};
            n_ch_idxs = length(ch_idxs);
            
            % Store for average across channels
            all_channels_ROI_power_wt = cell(n_ch_idxs, 1);
            all_channels_ROI_power_wt_no_norm = cell(n_ch_idxs, 1);
            all_channels_ROI_baseline = cell(n_ch_idxs, 1);

            psd_all_channels_ROI = cell(n_ch_idxs);
            psd_norm_all_channels_ROI = cell(n_ch_idxs);

            % Store for average psd across channels
            phase_psd_all_channels_ROI = cell(4, n_ch_idxs);
            phase_psd_norm_all_channels_ROI = cell(4, n_ch_idxs);

            % --- Loop sui canali della ROI ---%
            for ch = 1:n_ch_idxs
                
                ch_idx = ch_idxs(ch);
        
                % Baseline CWT 
                baseline_data = baseline_EEG(ch_idx, :);
                baseline_length = length(baseline_data);
                baseline_idxs = 1:baseline_length;

                % Global CWT computation 
                EEG_channel_data = EEG.data(ch_idx,:);
                [wt_glob, fr_glob] = cwt(EEG_channel_data, fs, 'morse');
                power_wt_glob = abs(wt_glob).^2;
        
                n_fr = size(fr_glob,1);
                
                % Initialize store for powers and times across trials
                all_wt_mat = zeros(n_fr, tot_length, num_trials);
                all_wt_mat_no_norm = zeros(n_fr, tot_length, num_trials);

                % Il segnale cinematico può essere processato una sola
                % volta per il primo segnale del blocco-->
                % Calcolo media segnale cinematico tra trials per
                % sovrapposizione

                if r == 1 && ch == 1

                    % Store frequency axis information (can be done one
                    % time only since signal length is the same across
                    % channels in a block)
                    all_subj_freq{f} = fr_glob;
                    all_blocks_freq_length(1,f) = n_fr;

                    for t = 1:num_trials
                        kin_signal = KIN_trials{t};
                        
                        % Store for interpolated kinematic signal
                        kin_trial_interp = zeros(tot_length, 1);

                        for p = 1:size(EEG_phases,2)
                            
                            phase_data = EEG_phases{t,p}(ch_idx,:);
                            phase_length = length(phase_data);
                            phase_median_length = median_lengths_phase(p);

                            % Segnale cinematico da sovraimporre 
                            if p == 1
                                start_kin = 1;
                            else 
                                start_kin = end_kin + 1;
                            end
                            end_kin = start_kin + phase_length - 1;
                            kin_phase = kin_signal(start_kin:end_kin);
                            
                            x_new_kin = linspace(1, length(kin_phase), phase_median_length);
                            kin_phase_interp = interp1((1:length(kin_phase)), kin_phase(:), x_new_kin, 'linear', 'extrap');
            
                            if p == 1
                                start_phase_interp = 1;
                            else 
                                start_phase_interp = end_phase_interp + 1;
                            end
                            end_phase_interp = start_phase_interp + phase_median_length - 1;

                            kin_trial_interp(start_phase_interp:end_phase_interp) = kin_phase_interp;
                        end
                        % Calcolo KIN_trial normalizzato
                        kin_norm = (kin_trial_interp / max(kin_trial_interp)) * (max(fr_glob)-min(fr_glob)) + min(fr_glob);

                        % Store current trial kinematic signal
                        all_kin_norm_block(:,t) = kin_norm;
                    end

                    % Mean of kin_velocity across trials
                    kin_mean_block = mean(all_kin_norm_block, 2);
            
                    onsets = [onsets, onset];
                    pulls = [pulls, pull];
                    grasps = [grasps, grasp];
                    offsets = [offsets, offset];

                    all_blocks_kin{f} = kin_mean_block;
                end

                all_trials_baseline = zeros(num_trials, n_fr);
                
                psd_all_trials = zeros(n_fr,num_trials);
                psd_norm_all_trials = zeros(n_fr,num_trials);

                phase_psd_all_trials = zeros(4, n_fr, num_trials);
                phase_psd_norm_all_trials = zeros(4, n_fr, num_trials);

                for t = 1:num_trials
                    trial_data = EEG_trials{t}(ch_idx,:);
                    trial_length = length(trial_data);
        
                    % Store for current trial power and kinematic signal
                    pow_wt_trial = zeros(n_fr,tot_length);
                    pow_wt_trial_no_norm = zeros(n_fr,tot_length);

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
                            all_trials_baseline(t,:) = mean_baseline_trial;
                    
                            % Salvataggio PSD non normalizzata del rest
                            phase_psd_all_trials(1,:,t) = mean(pow_rest, 2);

                            pow_rest_norm = (pow_rest - mean_baseline_trial)./mean_baseline_trial;
                            phase_psd_norm_all_trials(p,:,t) = mean(pow_rest_norm,2);

                            pow_wt_trial_original = power_wt_glob(:,start_idx:trial_length+start_idx-1);
                            psd_all_trials(:,t) = mean(pow_wt_trial_original,2);

                            pow_wt_trial_original_norm = (pow_wt_trial_original - mean_baseline_trial)./mean_baseline_trial;
                            psd_norm_all_trials(:,t) = mean(pow_wt_trial_original_norm,2);
        
                        elseif p > 1 && p < 5
                            % Normalization by trial baseline
                            pow_wt_phase = (pow_wt_phase_no_norm - mean_baseline_trial)./mean_baseline_trial;

                            phase_psd_all_trials(p,:,t) = mean(pow_wt_phase_no_norm,2);
                            phase_psd_norm_all_trials(p,:,t) = mean(pow_wt_phase,2);
                        end

                        % Normalization by trial baseline (All phases -->
                        % also rest_bef and rest aft)
                        pow_wt_phase = (pow_wt_phase_no_norm - mean_baseline_trial)./mean_baseline_trial;

                        % Interpolo la potenza alla lunghezza mediana
                        phase_median_length = median_lengths_phase(p);
                        x_new_wt_phase = linspace(1, phase_length, phase_median_length);
                        pow_wt_phase_interp = interp1((1:phase_length)', pow_wt_phase(:,:)', x_new_wt_phase', 'linear', 'extrap')';
                        
                        % Interpolo potenza non normalizzata
                        pow_wt_phase_no_norm_interp = interp1((1:phase_length)', pow_wt_phase_no_norm(:,:)', x_new_wt_phase', 'linear', 'extrap')';
        
                        if p == 1
                            start_phase_interp = 1;
                        else 
                            start_phase_interp = end_phase_interp + 1;
                        end
                        end_phase_interp = start_phase_interp + phase_median_length - 1;
                        pow_wt_trial(:, start_phase_interp:end_phase_interp) = pow_wt_phase_interp;

                        pow_wt_trial_no_norm(:, start_phase_interp:end_phase_interp) = pow_wt_phase_no_norm_interp;
                    end 

                    % Store current trial power
                    all_wt_mat(:,:,t) = pow_wt_trial;     
                    all_wt_mat_no_norm(:,:,t) = pow_wt_trial_no_norm;  

                end

                % Median power across trials
                power_wt_median = median(all_wt_mat, 3);
                power_wt_median_no_norm = median(all_wt_mat_no_norm, 3);
        
                % Store average normalized CWT for block and current
                % channel
                all_channels_ROI_power_wt{ch} = power_wt_median;
                all_channels_ROI_power_wt_no_norm{ch} = power_wt_median_no_norm;
                all_channels_ROI_baseline{ch} = mean(all_trials_baseline,1);

                for p = 1:n_phases
                    phase_psd_all_channels_ROI{p,ch} = median(phase_psd_all_trials(p,:,:),3); % Median across trials
                    phase_psd_norm_all_channels_ROI{p,ch} = median(phase_psd_norm_all_trials(p,:,:),3);
                end

                psd_all_channels_ROI{ch} = median(psd_all_trials,2);
                psd_norm_all_channels_ROI{ch} = median(psd_norm_all_trials,2);

            end
            
            % Matrice Power TF per tutti i canali della ROI corrente
            channel_power_wt_mat  = zeros(n_fr, tot_length, n_ch_idxs);
            channel_power_wt_mat_no_norm  = zeros(n_fr, tot_length, n_ch_idxs);
            channel_baseline_mat  = zeros(n_fr, n_ch_idxs);

            % Matrice Phases PSD per tutti i canali della ROI corrente
            channel_psd_mat = zeros(4, n_fr, n_ch_idxs);
            channel_psd_norm_mat = zeros(4, n_fr, n_ch_idxs);

            channel_psd_mat_trial = zeros(n_fr, n_ch_idxs);
            channel_psd_norm_mat_trial = zeros(n_fr, n_ch_idxs);

            % Average Power TF and phases over channels of the current ROI
            for c = 1:n_ch_idxs
                channel_power_wt_mat(:,:,c) = all_channels_ROI_power_wt{c};
                channel_power_wt_mat_no_norm(:,:,c) = all_channels_ROI_power_wt_no_norm{c};
                channel_baseline_mat(:,c) = all_channels_ROI_baseline{c};

                for p = 1:n_phases
                    channel_psd_mat(p,:,c) = phase_psd_all_channels_ROI{p,c};
                    channel_psd_norm_mat(p,:,c) = phase_psd_norm_all_channels_ROI{p,c};
                end

                channel_psd_mat_trial(:,c) = psd_all_channels_ROI{c};
                channel_psd_norm_mat_trial(:,c) = psd_norm_all_channels_ROI{c};
            end

            mean_power_wt_ROI = mean(channel_power_wt_mat, 3);
            mean_power_wt_ROI_no_norm = mean(channel_power_wt_mat_no_norm, 3);
            mean_baseline_ROI = mean(channel_baseline_mat, 2);
            
            for p = 1:n_phases
                all_blocks_ROI_psd{p,r,f} = mean(channel_psd_mat(p,:,:),3);
                all_blocks_ROI_psd_norm{p,r,f} = mean(channel_psd_norm_mat(p,:,:),3);
            end

            all_blocks_ROI_psd_trial{r,f} = mean(channel_psd_mat_trial,2);
            all_blocks_ROI_psd_norm_trial{r,f} = mean(channel_psd_norm_mat_trial,2);

            all_blocks_ROI_power_wt{r,f} = mean_power_wt_ROI;
            all_blocks_ROI_power_wt_no_norm{r,f} = mean_power_wt_ROI_no_norm;
            all_blocks_ROI_baseline{r,f} = mean_baseline_ROI;
        end
    end

    % Media latencies eventi kinematici per soggetto (across blocks)
    mean_subj_onset  = mean(onsets);
    mean_subj_grasp = mean(grasps);
    mean_subj_pull   = mean(pulls);
    mean_subj_offset = mean(offsets);

    % Trovo Lunghezza minima fr_glob per soggetto
    min_fr_glob = min(all_blocks_freq_length);
    min_fr_glob_mask = 1:min_fr_glob;

    % Fr glob tagliato a min_length
    subj_fr = all_subj_freq{1}(1:min_fr_glob);

    % Store for time length over all blocks
    block_wt_length = zeros(1, n_blocks);

    % Store Matrice avg baseline per average su soggetto
    subj_ROI_baseline_mat = zeros(min_fr_glob, n_blocks, n_ROI);
    % Average of baseline across blocks for each ROI
    mean_ROI_subj_baseline = zeros(min_fr_glob, n_ROI);
    
    for r = 1:n_ROI
        for b = 1:n_blocks
            % Lunghezze temporali potenza e costruzione matrice baseline
            block_wt_length(b) = size(all_blocks_ROI_power_wt{1,b}, 2);
            subj_ROI_baseline_mat(:,b,r) = all_blocks_ROI_baseline{r,b}(1:min_fr_glob);
        end
        mean_ROI_subj_baseline(:,r) = mean(subj_ROI_baseline_mat(:,:,r), 2);
    end

    % Lunghezza temporale media tra blocchi
    mean_length_blocks = round(mean(block_wt_length));

    % Interpolazione segnale cinematico alla lunghezza media tra blocchi del soggetto
    all_kin_subj = zeros(mean_length_blocks, n_blocks);
    
    for k = 1:n_blocks
        % Interpolazione KIN alla lunghezza media tra blocchi
        x_new_kin = linspace(1, length(all_blocks_kin{k}), mean_length_blocks);
        kin_block_interp = interp1((1:length(all_blocks_kin{k}))', all_blocks_kin{k}', ...
                                   x_new_kin', 'linear', 'extrap');
        all_kin_subj(:,k) = kin_block_interp;
    end

    % Media segnale cinematico tra blocchi
    mean_kin_subj = mean(all_kin_subj, 2);

    % Interpolazione Power TF per ROI per blocco alla lunghezza media tra blocchi del soggetto
    all_ROI_power_wt_subj = zeros(min_fr_glob, mean_length_blocks, n_blocks, n_ROI);
    all_ROI_mean_power_wt_subj = zeros(min_fr_glob, mean_length_blocks, n_ROI);
    
    all_ROI_power_wt_subj_no_norm = zeros(min_fr_glob, mean_length_blocks, n_blocks, n_ROI);
    all_ROI_mean_power_wt_subj_no_norm = zeros(min_fr_glob, mean_length_blocks, n_ROI);

    all_ROI_psd_subj = zeros(n_phases,min_fr_glob, n_ROI, n_blocks);
    all_ROI_mean_psd_subj = zeros(n_phases,min_fr_glob, n_ROI);
    all_ROI_psd_norm_subj = zeros(n_phases, min_fr_glob, n_ROI, n_blocks);
    all_ROI_mean_psd_norm_subj = zeros(n_phases, min_fr_glob, n_ROI);

    all_ROI_psd_subj_trial = zeros(min_fr_glob, n_ROI, n_blocks);
    all_ROI_psd_norm_subj_trial = zeros(min_fr_glob, n_ROI, n_blocks);

    for r=1:n_ROI

        ROI_name = ROI_names{r};
        for k = 1:n_blocks
            x_new_wt = linspace(1, block_wt_length(k), mean_length_blocks);
            power_wt_block_interp = interp1((1:block_wt_length(k))', all_blocks_ROI_power_wt{r,k}', x_new_wt', 'linear', 'extrap')';
            all_ROI_power_wt_subj(:,:,k,r) =  power_wt_block_interp(min_fr_glob_mask,:);

            power_wt_block_interp_no_norm = interp1((1:block_wt_length(k))', all_blocks_ROI_power_wt_no_norm{r,k}', x_new_wt', 'linear', 'extrap')';
            all_ROI_power_wt_subj_no_norm(:,:,k,r) =  power_wt_block_interp_no_norm(min_fr_glob_mask,:);

            for p = 1:n_phases
                all_ROI_psd_subj(p,:,r,k) = all_blocks_ROI_psd{p,r,k}(1:min_fr_glob);
                all_ROI_psd_norm_subj(p,:,r,k) = all_blocks_ROI_psd_norm{p,r,k}(1:min_fr_glob);
            end

            all_ROI_psd_subj_trial(:,r,k) = all_blocks_ROI_psd_trial{r,k}(1:min_fr_glob);
            all_ROI_psd_norm_subj_trial(:,r,k) = all_blocks_ROI_psd_norm_trial{r,k}(1:min_fr_glob);
            
        end
        % Average across blocks of the interpolated Power TF for ROI
        all_ROI_mean_power_wt_subj(:,:,r) = mean(all_ROI_power_wt_subj(:,:,:,r), 3);
        all_ROI_mean_power_wt_subj_no_norm(:,:,r) = mean(all_ROI_power_wt_subj_no_norm(:,:,:,r), 3);

        % Average across blocks of the phases PSD for ROI
        for p = 1:n_phases
            all_ROI_mean_psd_subj(p,:,r) = mean(all_ROI_psd_subj(p,:,r,:),4);
            all_ROI_mean_psd_norm_subj(p,:,r) = mean(all_ROI_psd_norm_subj(p,:,r,:),4);
            all_subj_psd{p,r,s} = squeeze(all_ROI_mean_psd_subj(p,:,r));
            all_subj_psd_norm{p,r,s} = squeeze(all_ROI_mean_psd_norm_subj(p,:,r));
        end

        all_subj_psd_trial{r,s} = squeeze(mean(all_ROI_psd_subj_trial(:,r,:),3));
        all_subj_psd_norm_trial{r,s} = squeeze(mean(all_ROI_psd_norm_subj_trial(:,r,:),3));

        % % === SAVE PSD DATA FOR CURRENT SUBJECT ===
        % save_folder = fullfile(base_path, 'RESULTS','Baseline_EEG_03_median','Single_subject_PSD');
        % save_filename = sprintf('PSD_Subj_%s.mat', subject_id);
        % save_path = fullfile(save_folder, save_filename);
        % 
        % % Struct for subject PSD saving
        % PSD_data = struct();
        % PSD_data.subject_id = subject_id;
        % PSD_data.freq = subj_fr;
        % PSD_data.phase_names = phase_names;
        % PSD_data.ROI_names = ROI_names;
        % PSD_data.PSD_trial = all_subj_psd_trial(:,s);
        % PSD_data.PSD_norm_trial = all_subj_psd_norm_trial(:,s);
        % PSD_data.PSD = all_ROI_mean_psd_subj;            % [n_phases x n_freq x n_ROI]
        % PSD_data.PSD_norm = all_ROI_mean_psd_norm_subj;  % [n_phases x n_freq x n_ROI]
        % PSD_data.baseline = mean_ROI_subj_baseline;      % [n_freq x n_ROI]
        % PSD_data.fs = fs;
        % 
        % save(save_path, 'PSD_data', '-v7.3');

        
        % Salvo baseline media per soggetto e ROI e fr_glob minima per soggetto
        all_ROI_med_subj_baseline{r,s} = mean_ROI_subj_baseline(:,r);
        all_subj_wt{r,s} =  all_ROI_mean_power_wt_subj(:,:,r);
        all_subj_wt_no_norm{r,s} =  all_ROI_mean_power_wt_subj_no_norm(:,:,r);
    
        % Salvo freq e segnale cinematico medio per soggetto e ROI
        all_subj_fr{s} = subj_fr;
        all_subj_kin{s} = mean_kin_subj;
    end

    % Salvo latencies medie eventi kin
    all_kin_latencies(1, s) = mean_subj_onset;
    all_kin_latencies(2, s) = mean_subj_grasp;
    all_kin_latencies(3, s) = mean_subj_pull;
    all_kin_latencies(4, s) = mean_subj_offset;
end


%% --- PLOT BASELINE MEDIA ACROSS SUBJECTS --- %

% Lunghezza minima in frequenza
min_freqs = min(cellfun(@(x) size(x,1), all_subj_wt(1,:)));
minimum_fr_length = min_freqs;

for r = 1:n_ROI

    % Costruisco matrice baseline e psd fasi per average tra soggetti
    all_baseline_mat = zeros(minimum_fr_length, n_subjects);
    all_phases_psd_mat = zeros(4,minimum_fr_length, n_subjects);
    all_phases_psd_norm_mat = zeros(4,minimum_fr_length, n_subjects);
    for s = 1:n_subjects
        all_baseline_mat(:,s) = all_ROI_med_subj_baseline{r,s}(1:minimum_fr_length);
        for p = 1:n_phases
            all_phases_psd_mat(p,:,s) = all_subj_psd{p,r,s}(1:minimum_fr_length);
            all_phases_psd_norm_mat(p,:,s) = all_subj_psd_norm{p,r,s}(1:minimum_fr_length);
        end
    end

    % Frequenze comuni (prendo quelle del primo soggetto tagliate alla min len)
    common_fr = all_subj_fr{1}(1:minimum_fr_length);

    % Media + SEM tra soggetti
    mean_baseline = mean(all_baseline_mat, 2);
    sem_baseline  = std(all_baseline_mat, 0, 2) ./ sqrt(n_subjects);

    % Plot
    figure; hold on;
    fill([common_fr' fliplr(common_fr')], ...
         [mean_baseline'+sem_baseline' fliplr(mean_baseline'-sem_baseline')], ...
         [0.8 0.8 1], 'EdgeColor','none', 'FaceAlpha',0.4);
    plot(common_fr, mean_baseline, 'b', 'LineWidth',2);

    xlabel('Frequency (Hz)', 'FontSize',12);
    ylabel('Normalized Power (a.u.)', 'FontSize',12);
    title([' Average across subjects - Baseline Spectrum ± SEM - ROI : ',ROI_names{r}], 'FontSize',14);
    grid on; set(gca,'FontSize',12);

    % Salvataggio in PNG e FIG
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 600]); 
    set(gca, 'FontSize', 12); % dimensione testo
    fig_name = sprintf(['Temporal_Mean_Baseline_all_subj_', ROI_names{r}]);
    out_dir = fullfile(base_path,'RESULTS_final', 'Baseline_EEG_03_median');
    print(fullfile(out_dir, fig_name), '-dpng', '-r300'); 
    savefig(fullfile(out_dir, [fig_name, '.fig']));

    % === Colori ===
    col_baseline = [0.7 0.7 0.7];    
    col_rest     = [0.8 0.2 0.2];      
    col_reach    = [0.3 0.8 0.3];      
    col_grasp    = [0.2 0.6 1];   
    col_pull     = [0.6 0.1 0.9]; 
    phase_colors = {col_rest, col_reach, col_grasp, col_pull};


    % =========================
    % 1) PLOT PSD NON NORMALIZZATA
    % ============================

    mean_phases_psd = zeros(4,minimum_fr_length);
    mean_phases_psd_norm = zeros(4,minimum_fr_length);
    sem_phases_psd = zeros(4,minimum_fr_length);
    sem_phases_psd_norm = zeros(4,minimum_fr_length);
    for p = 1:n_phases
        mean_phases_psd(p,:) = mean(all_phases_psd_mat(p,:,:), 3);
        mean_phases_psd_norm(p,:) = mean(all_phases_psd_norm_mat(p,:,:), 3);
        sem_phases_psd(p,:) = std(all_phases_psd_mat(p,:,:),0, 3) ./ sqrt(n_subjects);
        sem_phases_psd_norm(p,:) = std(all_phases_psd_norm_mat(p,:,:), 0, 3) ./ sqrt(n_subjects);
    end

    figure('Name',['PSD non normalizzata ', ROI_names{r}],'Position',[100 100 1200 600]); hold on; grid on;

    common_fr = common_fr(:)';  
    grand_mean_baseline = mean_baseline(:)';  
    sem_baseline = sem_baseline(:)';

    % % Baseline
    % h_base_fill = fill([common_fr fliplr(common_fr)], ...
    %      [grand_mean_baseline+sem_baseline fliplr(grand_mean_baseline-sem_baseline)], ...
    %      col_baseline, 'EdgeColor','none', 'FaceAlpha',0.3, 'DisplayName', 'SEM Baseline');
    % plot(common_fr, grand_mean_baseline, 'Color', col_baseline, 'LineWidth', 2.5, 'DisplayName', 'Baseline');

    % Fasi cinematiche
    for p = 1:n_phases
        mean_ph = squeeze(mean_phases_psd(p,:)*100);
        sem_ph  = squeeze(sem_phases_psd(p,:)*100);
        mean_ph = mean_ph(:)'; sem_ph = sem_ph(:)';

        h_fill = fill([common_fr fliplr(common_fr)], ...
             [mean_ph+sem_ph fliplr(mean_ph-sem_ph)], ...
             phase_colors{p}, 'EdgeColor','none', 'FaceAlpha',0.25,'HandleVisibility','off');
        plot(common_fr, mean_ph, 'Color', phase_colors{p}, 'LineWidth', 2.5, 'DisplayName', [phase_names{p}]);
    end

    xlabel('Frequency (Hz)', 'FontSize', 12);
    ylabel('Power (a.u.)', 'FontSize', 12);
    xlim([0 80]);
    title(['Grand-Average EEG Marginal Power Spectral Density - ', ROI_names{r}], 'FontSize', 14);
    legend('Location','best');
    set(gca, 'FontSize', 12);

    % Salvataggio
    fig_name = ['Mean_PSD_all_subjects',ROI_names{r}];
    print(fullfile(out_dir, fig_name), '-dpng', '-r300');
    saveas(gcf, fullfile(out_dir, [fig_name '.fig']));


    %% =========================
    % 2) PLOT PSD NORMALIZZATA
    % ============================
    figure('Name',['PSD normalizzata', ROI_names{r}],'Position',[100 100 1200 600]); hold on; grid on;

    % Fasi normalizzate
    for p = 1:length(phase_names)
        mean_ph_n = squeeze(mean_phases_psd_norm(p,:)*100);
        sem_ph_n  = squeeze(sem_phases_psd_norm(p,:)*100);
        mean_ph_n = mean_ph_n(:)'; sem_ph_n = sem_ph_n(:)';

        fill([common_fr fliplr(common_fr)], ...
         [mean_ph_n+sem_ph_n fliplr(mean_ph_n-sem_ph_n)], ...
         phase_colors{p}, 'EdgeColor','none', 'FaceAlpha',0.25, 'HandleVisibility','off');

        plot(common_fr, mean_ph_n, 'Color', phase_colors{p}, 'LineWidth', 2.5, 'DisplayName', [phase_names{p}]);
    end

    xlabel('Frequency (Hz)', 'FontSize', 12);
    ylabel('Relative Power (w.r.t. Rest)(%)', 'FontSize', 12);
    xlim([0 80]);
    ylim([-50 40]);
    title(['Grand-Average EEG Marginal Normalized Power (w.r.t. Rest) - ' ROI_names{r}], 'FontSize', 14);
    legend('Location','best');
    set(gca, 'FontSize', 12);

    % Salvataggio
    fig_name = ['Mean_PSD_norm_all_subjects',ROI_names{r}];
    print(fullfile(out_dir, fig_name), '-dpng', '-r300');
    saveas(gcf, fullfile(out_dir, [fig_name '.fig']));


end

% --- PLOT AVERAGE su soggetti POTENZA NORMALIZZATA TF per ROI --- %
% Chiedere riguardo il SEM per plot Vissani : fatto media temporale per
% fase di matrice SEM TF.

% Lunghezza mediana temporale per interpolazione
median_length_wt = round(median(cellfun(@(x) size(x,2), all_subj_wt(1,:))));

% Matrici 4D 
all_wt_mat  = zeros(min_freqs, median_length_wt, n_ROI, n_subjects);
all_wt_mat_no_norm  = zeros(min_freqs, median_length_wt, n_ROI, n_subjects);
all_kin_mat = zeros(median_length_wt, n_subjects);

% Frequency axis 
length_subj_fr = zeros(1,n_subjects);
for s = 1:n_subjects
    length_subj_fr = all_subj_fr{s};
end
[~, idx_max] = max(length_subj_fr);
all_fr_glob =  all_subj_fr{idx_max};

for r = 1:n_ROI
    for s = 1:n_subjects
        wt_subj  = all_subj_wt{r,s};
        wt_subj_no_norm  = all_subj_wt_no_norm{r,s};
        kin_subj = all_subj_kin{s};

        % Lunghezza originale e asse temporale
        or_length  = size(wt_subj,2);
        t_orig = 1:or_length;

        % Asse temporale target
        t_target = linspace(1,or_length,median_length_wt);

        % interpolazione lungo la dimensione temporale
        wt_interp = interp1(t_orig, wt_subj(1:min_freqs,:)', t_target, 'linear', 'extrap')';
        wt_interp_no_norm = interp1(t_orig, wt_subj_no_norm(1:min_freqs,:)', t_target, 'linear', 'extrap')';
        kin_interp = interp1(t_orig, kin_subj(:), t_target, 'linear', 'extrap');

        % salvo nelle matrici 3D
        all_wt_mat(:,:,r,s)  = wt_interp;
        all_wt_mat_no_norm(:,:,r,s)  = wt_interp_no_norm;
        all_kin_mat(:,s)   = kin_interp;
    end
    % Media Potenza e segnale kinematico tra soggetti
    all_subj_avg_wt  = mean(all_wt_mat(:,:,r,:), 4);
    all_subj_avg_wt_no_norm  = mean(all_wt_mat_no_norm(:,:,r,:), 4);
    all_subj_avg_kin = mean(all_kin_mat, 2);

    % SEM across subjects 
    sem_wt = std(all_wt_mat, 0, 3) ./ sqrt(n_subjects);
    sem_wt_no_norm = std(all_wt_mat_no_norm, 0, 3) ./ sqrt(n_subjects);

    % Calcolo asse tempi e frequenze
    time_axis = (1:median_length_wt) / fs;
    freq_axis = all_fr_glob(1:min_freqs);

    % Average latencies kinematics events
    mean_onset_all = mean(all_kin_latencies(1,:), 2);
    mean_grasp_all = mean(all_kin_latencies(2,:), 2);
    mean_pull_all = mean(all_kin_latencies(3,:), 2);
    mean_offset_all = mean(all_kin_latencies(4,:), 2);

    % --- Limit frequency axis to 0-90 Hz ---
    freq_limit = 90;
    freq_mask = freq_axis <= freq_limit;
    
    freq_axis_plot = freq_axis(freq_mask);
    all_subj_avg_wt_plot = all_subj_avg_wt(freq_mask, :);


    % --- Normalize kinematic signal to the 0-90 Hz frequency range ---
    kin_min = min(freq_axis_plot);
    kin_max = max(freq_axis_plot);
    all_subj_avg_kin_norm = (all_subj_avg_kin - min(all_subj_avg_kin)) ...
                            / (max(all_subj_avg_kin) - min(all_subj_avg_kin)); % 0-1
    all_subj_avg_kin_norm = all_subj_avg_kin_norm * (kin_max - kin_min) + kin_min;


    % Plot Grand Average WT
    figure('Name',['All Subjects Average Normalized Power CWT - ROI - ', ROI_names{r}],'NumberTitle','off');
    pcolor(time_axis, freq_axis_plot, all_subj_avg_wt_plot*100);
    shading interp; set(gca,'YScale','linear','FontSize',12);
    xlabel('Time (s)','FontSize',12);
    ylabel('Frequency (Hz)','FontSize',12);
    xline(mean_onset_all,  'Label','Mov-ONSET',   'Color','w','LineWidth',2.5);
    xline(mean_grasp_all,  'Label','Grasp Start', 'Color','w','LineWidth',2.5);
    xline(mean_pull_all,   'Label','Grasp End',   'Color','w','LineWidth',2.5);
    xline(mean_offset_all, 'Label','Mov-OFFSET',  'Color','w','LineWidth',2.5);
    title(['Grand-Average EEG Time–Frequency Power (CWT, Rest-Normalized) - ROI - ', ROI_names{r}],'FontSize',14);
    c = colorbar; ylabel(c,'Relative Power (w.r.t. Rest) (%)');
    colormap jet;
    alpha_caxis_grand = quantile(all_subj_avg_wt_plot(:)*100, 0.70);
    caxis([min(all_subj_avg_wt_plot(:)*100) max(all_subj_avg_wt_plot(:)*100)]);
    %clim('auto');

    % Sovrappongo segnale cinematico medio
    yyaxis right
    plot(time_axis, all_subj_avg_kin_norm, 'k','LineWidth',2);
    ax = gca;
    ax.YAxis(2).Visible = 'off';

    %ylabel('Grand average kinematic signal');

    % Salvataggio
    set(gcf,'Units','pixels','Position',[100,100,1200,600]);
    fig_name = ['All_Subjects_Average_Power_TF_',ROI_names{r}];
    out_dir_tf = fullfile(base_path,'RESULTS_final','EEG_TF_pow_norm_phase_results_03_median','TF_Plot');
    if ~exist(out_dir_tf, 'dir')
        mkdir(out_dir_tf);
    end

    print(fullfile(out_dir_tf,fig_name),'-dpng','-r300');
    savefig(fullfile(out_dir_tf,[fig_name '.fig']));
    
    % % Salva potenza media e segnale cinematico per plot Vissani
    % save(fullfile(base_path, 'RESULTS/EEG_TF_pow_norm_phase_results_03_median/Mat_data', ['Final_Results_All_Subj_' ROI_names{r} '.mat']), ...
    %      'all_subj_avg_wt', ...          % Potenza media TF ROI sui soggetti
    %      'all_subj_avg_kin', ...         % Segnale cinematico medio
    %      'sem_wt', ...                   % SEM della potenza
    %      'freq_axis', ...                % Asse delle frequenze
    %      'time_axis', ...                % Asse temporale
    %      'mean_onset_all', 'mean_grasp_all', 'mean_pull_all', 'mean_offset_all', ... % Latencies medie
    %      'n_subjects', ...               % Numero di soggetti
    %      'ROI_names' );                  % Nome delle ROI
    % 
    % % Salva potenza media NON NORMALIZZATA e segnale cinematico per plot Vissani
    % save(fullfile(base_path, 'RESULTS/EEG_TF_pow_norm_phase_results_03_median/Mat_data', ['Final_Results_All_Subj_' ROI_names{r} '_NO_NORM.mat']), ...
    %      'all_subj_avg_wt_no_norm', ...          % Potenza media TF ROI sui soggetti
    %      'all_subj_avg_kin', ...         % Segnale cinematico medio
    %      'sem_wt_no_norm', ...                   % SEM della potenza
    %      'freq_axis', ...                % Asse delle frequenze
    %      'time_axis', ...                % Asse temporale
    %      'mean_onset_all', 'mean_grasp_all', 'mean_pull_all', 'mean_offset_all', ... % Latencies medie
    %      'n_subjects', ...               % Numero di soggetti
    %      'ROI_names' );                  % Nome delle ROI

    %%  NOT NORMALIZED TF POWER

    all_subj_avg_wt_no_norm_plot = all_subj_avg_wt_no_norm(freq_mask, :);


    % Plot Grand Average WT
    figure('Name',['All Subjects Average Not Normalized Power CWT - ROI - ', ROI_names{r}],'NumberTitle','off');
    pcolor(time_axis, freq_axis_plot, all_subj_avg_wt_no_norm_plot);
    shading interp; set(gca,'YScale','linear','FontSize',12);
    xlabel('Time (s)','FontSize',12);
    ylabel('Frequency (Hz)','FontSize',12);
    xline(mean_onset_all,  'Label','Mov-ONSET',   'Color','w','LineWidth',2.5);
    xline(mean_grasp_all,  'Label','Grasp Start', 'Color','w','LineWidth',2.5);
    xline(mean_pull_all,   'Label','Grasp End',   'Color','w','LineWidth',2.5);
    xline(mean_offset_all, 'Label','Mov-OFFSET',  'Color','w','LineWidth',2.5);
    title(['Grand-Average EEG Time–Frequency Power (CWT, Not Normalized) - ROI - ', ROI_names{r}],'FontSize',14);
    c = colorbar; ylabel(c,'Power (a.u.)');
    colormap jet;
    alpha_caxis_grand = quantile(all_subj_avg_wt_no_norm_plot(:), 0.85);
    caxis([min(all_subj_avg_wt_no_norm_plot(:)) alpha_caxis_grand]);
    %clim('auto');

    % Sovrappongo segnale cinematico medio
    yyaxis right
    plot(time_axis, all_subj_avg_kin_norm, 'k','LineWidth',2);
    ax = gca;
    ax.YAxis(2).Visible = 'off';

    % Salvataggio
    set(gcf,'Units','pixels','Position',[100,100,1200,600]);
    fig_name = ['All_Subjects_Average_Not_Norm_Power_TF_',ROI_names{r}];
    out_dir_tf_no_norm = fullfile(base_path,'RESULTS_final','EEG_TF_pow_norm_phase_results_03_median','TF_Plot_No_Norm');
    if ~exist(out_dir_tf_no_norm, 'dir')
        mkdir(out_dir_tf_no_norm);
    end

    print(fullfile(out_dir_tf_no_norm,fig_name),'-dpng','-r300');
    savefig(fullfile(out_dir_tf_no_norm,[fig_name '.fig']));
end


