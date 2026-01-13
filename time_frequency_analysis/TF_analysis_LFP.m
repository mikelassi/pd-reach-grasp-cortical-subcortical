%% ==========================================
%  LFP Time-Frequency Analysis
%%  Interpolation for each PHASE - TRIAL-WISE NORMALIZATION
%   1) Compute the CWT (Morse wavelet) on the block
%   2) Normalize power with respect to the REST phase
%   3) Identify behavioral regions (rest_bef, reach,
%      grasp, pull, rest_aft) and extract median length for each phase
%   4) Interpolate to the median length for each phase
%   5) Apply the same processing to the kinematic signal
%
%  Author: Tommaso Marcantoni
%% ==========================================

close all;
subject_list = {'wue02', 'wue03', 'wue05', 'wue06', 'wue07', 'wue09', 'wue10', 'wue11'}; % List of subjects
n_subjects = length(subject_list);
phase_names = {'Rest', 'Reach', 'Grasp', 'Pull'};
n_phases = length(phase_names);

base_path    = 'C:\Users\tomma\OneDrive - University of Pisa\Desktop\TESI\Dataset_tesi';
preproc_subfolder = 'Preprocessed_01\LFP_01';

eeglab;

all_med_subj_baseline = cell(1, n_subjects);
all_med_phases_psd_subj  = cell(4, n_subjects);
all_med_phases_psd_norm_subj  = cell(4, n_subjects);
all_subj_fr = cell(1, n_subjects);

all_med_psd_subj  = cell(1, n_subjects);
all_med_psd_norm_subj  = cell(1, n_subjects);

% Store WT and KIN signal for plot across subjects
all_subj_kin = cell(1, n_subjects);
all_subj_wt = cell(1, n_subjects);
all_subj_wt_no_norm = cell(1, n_subjects);

% Store for kinematic latencies
all_kin_latencies = zeros(4,n_subjects);

for s = 1:n_subjects
    
    subject_id = subject_list{s};
    preproc_path = fullfile(base_path, subject_id, preproc_subfolder);
    mat_files = dir(fullfile(preproc_path, '*_LFP_trialsByRegionAndPhase.mat'));
    set_files = dir(fullfile(preproc_path, '*.set'));
    n_blocks = length(mat_files);

    % Define channel of interest ( CONTROLATERAL to hand performing
    % movement)
    if any(strcmp(subject_id, {'wue02', 'wue03'})) 
        channel_of_interest = 2;
    else
        channel_of_interest = 1;
    end

    % Store for kinematic event latencies across blocks
    onsets = [];
    pulls = [];
    grasps = [];
    offsets = [];

    % Store for Power and Kinematic signal average across blocks
    all_blocks_power_wt = cell(n_blocks, 1);
    all_blocks_power_wt_no_norm = cell(n_blocks, 1);
    all_blocks_kin = cell(n_blocks, 1);

    all_subj_baseline = cell(1, n_blocks);
    phases_psd_subj  = cell(4, n_blocks);
    phases_psd_norm_subj  = cell(4, n_blocks);

    all_blocks_psd_trial = cell(1, n_blocks);
    all_blocks_psd_norm_trial = cell(1, n_blocks);

    % Store for global frequency axis 
    all_subj_freq = cell(1, n_blocks);
    all_blocks_freq_length = zeros(1,n_blocks);

    tot_trials = 0;

    for f = 1:length(mat_files)
        mat_file = mat_files(f).name;
        fprintf('\nProcessing LFP segments: %s\n', mat_file);

        set_file = set_files(f).name;
        LFP = pop_loadset('filename', set_file, 'filepath', preproc_path);
        fs = LFP.srate;

        % Load the LFP cell array 
        M = load(fullfile(preproc_path, mat_file));
        
        % Load Baseline (From 50 samples after TenStrigger to A_T1)
        baseline_LFP = M.baseline_LFP;
        all_baseline_length = length(baseline_LFP(channel_of_interest, :));
       
        % Load Trials, Phases, Kinematic signal
        LFP_trials = M.LFP_trials;
        LFP_phases = M.LFP_phases;
        KIN_trials = M.kinematic_4LFP_trials;  

        % Number of trials of current block
        num_trials = size(LFP_trials,1);

        tot_trials = tot_trials + num_trials;

        % Compute median length of each phase across trials
        len_phase1 = zeros(num_trials,1);
        len_phase2 = zeros(num_trials,1);
        len_phase3 = zeros(num_trials,1);
        len_phase4 = zeros(num_trials,1);
        len_phase5 = zeros(num_trials,1);
        
        for t = 1:num_trials
            len_phase1(t) = size(LFP_phases{t,1}, 2);
            len_phase2(t) = size(LFP_phases{t,2}, 2);
            len_phase3(t) = size(LFP_phases{t,3}, 2);
            len_phase4(t) = size(LFP_phases{t,4}, 2);
            len_phase5(t) = size(LFP_phases{t,5}, 2);
        end

        % Median values
        median_len_1 = round(median(len_phase1));
        median_len_2 = round(median(len_phase2));
        median_len_3 = round(median(len_phase3));
        median_len_4 = round(median(len_phase4));
        median_len_5 = round(median(len_phase5));
        median_lengths_phase = [median_len_1, median_len_2, median_len_3, median_len_4, median_len_5];
        tot_length = median_len_1 + median_len_2 + median_len_3 + median_len_4 + median_len_5;

        % Initial baseline CWT 
        baseline_data = baseline_LFP(channel_of_interest, :);
        baseline_idxs = 1:all_baseline_length;

        all_data = LFP.data(channel_of_interest,:);
       
        % CWT of entire block and power computation
        [wt_glob, fr_glob] = cwt(all_data, fs, 'morse');
        power_wt_glob = abs(wt_glob).^2;

        n_fr = size(fr_glob,1);
        
        % Store frequency axis information
        all_subj_freq{f} = fr_glob;
        all_blocks_freq_length(1,f) = n_fr;

        % Initialize store for powers and times
        all_wt_mat = zeros(n_fr, tot_length, num_trials);
        all_wt_mat_no_norm = zeros(n_fr, tot_length, num_trials);
        all_kin_norm_block = zeros(tot_length, num_trials);

        phase_psd_all_trials = zeros(4, n_fr, num_trials);
        phase_psd_norm_all_trials = zeros(4, n_fr, num_trials);

        psd_all_trials = zeros(n_fr,num_trials);
        psd_norm_all_trials = zeros(n_fr,num_trials);

        all_trials_baseline = zeros(num_trials, n_fr);

        start_idx = 0;
        start_idx_2 = 0;
        for t = 1:num_trials
            trial_data = LFP_trials{t}(channel_of_interest,:);
            kin_signal = KIN_trials{t};

            trial_length = length(trial_data);

            % Store for current trial power and kinematic signal
            pow_wt_trial = zeros(n_fr,tot_length);
            pow_wt_trial_no_norm = zeros(n_fr,tot_length);
            kin_trial_interp = zeros(tot_length);

            n_lfp_phases = size(LFP_phases,2);

            for p = 1:n_lfp_phases

                % Data extraction for current phase
                phase_data = LFP_phases{t,p}(channel_of_interest,:);
                phase_length = length(phase_data);
            
                % Index relative to entire trials
                rest_offset_baseline = 250;
                if t == 1 && p == 1
                    start_idx = all_baseline_length - rest_offset_baseline;
                else
                    start_idx = end_idx + 1;
                end
                end_idx = start_idx + phase_length - 1;
            
                % Non normalized phase power
                pow_wt_phase_no_norm = power_wt_glob(:, start_idx:end_idx);

         
                % --- EXTRACTION REST PHASE for normalization
                if p == 1
                    pow_rest_bef = pow_wt_phase_no_norm;
                    length_rest_aft = length(LFP_phases{t,5}(channel_of_interest,:));
                    idx_rest_aft = (start_idx + trial_length - length_rest_aft): (start_idx + trial_length - 1);
                    pow_rest_aft = power_wt_glob(:,idx_rest_aft);
                    pow_rest = cat(2, pow_rest_bef, pow_rest_aft);
            
                    % PSD MEAN IN REST
                    mean_baseline_trial = mean(pow_rest, 2);
                    all_trials_baseline(t,:) = mean_baseline_trial;
            
                    % Save not normalized power in rest
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
                
                % Normalization for all phases
                pow_wt_phase = (pow_wt_phase_no_norm - mean_baseline_trial)./mean_baseline_trial;

                % --- INTERPOLATION OF POWER AND KINEMATIC SIGNAL ---
                phase_median_length = median_lengths_phase(p);
                x_new_wt_phase = linspace(1, phase_length, phase_median_length);
            
                pow_wt_phase_interp = interp1((1:phase_length)', pow_wt_phase(:,:)', x_new_wt_phase', 'linear', 'extrap')';
                pow_wt_phase_interp_no_norm = interp1((1:phase_length)', pow_wt_phase_no_norm(:,:)', x_new_wt_phase', 'linear', 'extrap')';
            
                % --- Kinematic signal associated with the phase ---
                if p == 1
                    start_kin = 1;
                else 
                    start_kin = end_kin + 1;
                end
                end_kin = start_kin + phase_length - 1;
                kin_phase = kin_signal(start_kin:end_kin);
            
                x_new_kin = linspace(1, length(kin_phase), phase_median_length);
                kin_phase_interp = interp1((1:length(kin_phase)), kin_phase(:), x_new_kin, 'linear', 'extrap');
            
                % --- Assemble the complete trial ---
                if p == 1
                    start_phase_interp = 1;
                else 
                    start_phase_interp = end_phase_interp + 1;
                end
                end_phase_interp = start_phase_interp + phase_median_length - 1;
            
                pow_wt_trial(:, start_phase_interp:end_phase_interp) = pow_wt_phase_interp;
                pow_wt_trial_no_norm(:, start_phase_interp:end_phase_interp) = pow_wt_phase_interp_no_norm;
                kin_trial_interp(start_phase_interp:end_phase_interp) = kin_phase_interp;
            end
            
            % Median length across phases
            onset = median_lengths_phase(1)/fs;
            grasp = onset + median_lengths_phase(2)/fs;
            pull = grasp + median_lengths_phase(3)/fs;
            offset = pull + median_lengths_phase(4)/fs;

            % KIN_trial normalized
            kin_norm = (kin_trial_interp / max(kin_trial_interp)) * (max(fr_glob)-min(fr_glob)) + min(fr_glob); 

            % Store current trial power
            all_wt_mat(:,:,t) = pow_wt_trial;
            all_wt_mat_no_norm(:,:,t) = pow_wt_trial_no_norm;

            % Store current trial kinematic signal
            all_kin_norm_block(:,t) = kin_norm;
        end

         % === Save PSD normalized per trial ===
        for t = 1:num_trials
            PSD_STRUCT.(subject_id).block(f).trial(t).psd_norm = psd_norm_all_trials(:,t);
            PSD_STRUCT.(subject_id).block(f).trial(t).psd = psd_all_trials(:,t);

            for p = 1:n_phases
                PSD_STRUCT.(subject_id).block(f).trial(t).phase(p).psd_norm = squeeze(phase_psd_norm_all_trials(p,:,t))';
                PSD_STRUCT.(subject_id).block(f).trial(t).phase(p).psd = squeeze(phase_psd_all_trials(p,:,t))';
                PSD_STRUCT.(subject_id).block(f).frequencies = fr_glob;
            end
        end

        
        % Store psd phases per block
        for p=1:n_phases
            phases_psd_subj{p,f} = median(phase_psd_all_trials(p,:,:),3);
            phases_psd_norm_subj{p,f} = median(phase_psd_norm_all_trials(p,:,:),3);
        end
        
        % Store psd baseline per block
        all_subj_baseline{f} = mean(all_trials_baseline,1);
        
        % Mean power across trials
        power_wt_mean = median(all_wt_mat, 3);
        power_wt_mean_no_norm = median(all_wt_mat_no_norm, 3);

        % Mean of kin_velocity across trials
        kin_mean_block = mean(all_kin_norm_block, 2);

        onsets = [onsets, onset];
        pulls = [pulls, pull];
        grasps = [grasps, grasp];
        offsets = [offsets, offset];

        % Store average normalized CWT for block
        all_blocks_power_wt{f} = power_wt_mean;
        all_blocks_power_wt_no_norm{f} = power_wt_mean_no_norm;
        all_blocks_kin{f} = kin_mean_block;

        all_blocks_psd_trial{f} = mean(psd_all_trials,2);
        all_blocks_psd_norm_trial{f} = mean(psd_norm_all_trials,2);
    end

    % Min frequency lenght per block
    min_fr_glob = min(all_blocks_freq_length);
    [max_fr_glob, idx_max_fr_glob] = max(all_blocks_freq_length);
    min_fr_glob_mask = 1:min_fr_glob;
    
    % Mean psd across blocks
    mat_psd_phases = zeros(4,min_fr_glob,n_blocks);
    mat_psd_norm_phases = zeros(4,min_fr_glob,n_blocks);

    mat_psd = zeros(min_fr_glob,n_blocks);
    mat_psd_norm = zeros(min_fr_glob,n_blocks);

    subj_baseline_mat = zeros(min_fr_glob, n_blocks);
    block_wt_length = zeros(1, n_blocks);

    for b = 1:n_blocks
        for p=1:n_phases
            mat_psd_phases(p,:,b) = phases_psd_subj{p,b}(1:min_fr_glob);
            mat_psd_norm_phases(p,:,b) = phases_psd_norm_subj{p,b}(1:min_fr_glob);
        end
        mat_psd(:,b) = all_blocks_psd_trial{b}(1:min_fr_glob);
        mat_psd_norm(:,b) = all_blocks_psd_norm_trial{b}(1:min_fr_glob);

        block_wt_length(b) = size(all_blocks_power_wt{b}, 2);
        subj_baseline_mat(:,b) = all_subj_baseline{b}(1:min_fr_glob);
    end

    % === Compute PSD median across all trials of all blocks ===
    all_psd_phases = cell(1, n_phases);
    all_psd_norm_phases = cell(1, n_phases);

    for f = 1:numel(PSD_STRUCT.(subject_id).block)
        for t = 1:numel(PSD_STRUCT.(subject_id).block(f).trial)
            for p = 1:n_phases
                psd_norm_phase = PSD_STRUCT.(subject_id).block(f).trial(t).phase(p).psd_norm(1:min_fr_glob);
                psd_phase = PSD_STRUCT.(subject_id).block(f).trial(t).phase(p).psd(1:min_fr_glob);

                all_psd_norm_phases{p} = [all_psd_norm_phases{p}, psd_norm_phase(:)];
                all_psd_phases{p} = [all_psd_phases{p}, psd_phase(:)];
            end
        end
    end


    % === Compute the median PSD across all trials (for each phase) ===
    for p = 1:n_phases
        all_med_phases_psd_norm_subj{p,s} = median(all_psd_norm_phases{p}, 2);
        all_med_phases_psd_subj{p,s} = median(all_psd_phases{p}, 2);
    end

    
    all_med_psd_subj{s} = mean(mat_psd,2);
    all_med_psd_norm_subj{s} = mean(mat_psd_norm,2);

    mean_subj_baseline = mean(subj_baseline_mat, 2);

    % Mean latencies kinematic events for subject
    mean_subj_onset  = mean(onsets);
    mean_subj_grasp = mean(grasps);
    mean_subj_pull   = mean(pulls);
    mean_subj_offset = mean(offsets);

    % Fr glob cut to min_length
    subj_fr = all_subj_freq{idx_max_fr_glob}(1:min_fr_glob);

    % Median temporal length across blocks
    mean_length_blocks = round(mean(block_wt_length));

    % Interpolation of kinematic signal and power to the mean length across blocks of the subject
    all_power_wt_block_interp = zeros(min_fr_glob, mean_length_blocks, n_blocks);
    all_power_wt_block_interp_no_norm = zeros(min_fr_glob, mean_length_blocks, n_blocks);
    all_kin_subj = zeros(mean_length_blocks, n_blocks);
    for k = 1:n_blocks
        x_new_wt_rest = linspace(1, block_wt_length(k), mean_length_blocks);
        power_wt_block_interp = interp1((1:block_wt_length(k))', all_blocks_power_wt{k}', x_new_wt_rest', 'linear', 'extrap')';
        power_wt_block_interp_no_norm = interp1((1:block_wt_length(k))', all_blocks_power_wt_no_norm{k}', x_new_wt_rest', 'linear', 'extrap')';

        all_power_wt_block_interp(:,:,k) =  power_wt_block_interp(min_fr_glob_mask,:);
        all_power_wt_block_interp_no_norm(:,:,k) =  power_wt_block_interp_no_norm(min_fr_glob_mask,:);

        % Interpolation KIN 
        x_new_kin = linspace(1, length(all_blocks_kin{k}), mean_length_blocks);
        kin_block_interp = interp1((1:length(all_blocks_kin{k}))', all_blocks_kin{k}', ...
                                   x_new_kin', 'linear', 'extrap');

        all_kin_subj(:,k) = kin_block_interp;
    end

    % Mean power interpolated to median length
    mean_wt_blocks_interp = mean(all_power_wt_block_interp, 3);
    mean_wt_blocks_interp_no_norm = mean(all_power_wt_block_interp_no_norm, 3);

    % Mean kinematic signal across blocks
    mean_kin_subj = mean(all_kin_subj, 2);

    % Plot TF
    figure('Name',['Avg Normalized Power wt - Subj - ' subject_id ' All Blocks '] ,'NumberTitle','off');
    pcolor((1:mean_length_blocks)/fs, subj_fr, mean_wt_blocks_interp(:,:)*100);
    shading interp; set(gca, 'YScale', 'linear', 'FontSize', 12);
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Frequency (Hz)', 'FontSize', 12);
    xline(mean_subj_onset,  'Label','Mov-ONSET',   'Color','w','LineWidth',2.5);
    xline(mean_subj_grasp,  'Label','Grasp-start', 'Color','w','LineWidth',2.5);
    xline(mean_subj_pull,   'Label','Grasp-end',   'Color','w','LineWidth',2.5);
    xline(mean_subj_offset, 'Label','Mov-OFFSET',  'Color','w','LineWidth',2.5);
    title(['Averaged normalized Power CWT - All Blocks - ' subject_id ' Ch ' num2str(channel_of_interest)]);
    c = colorbar; ylabel(c, 'Power wrt Baseline %');
    colormap jet;
    alpha_caxis_subj = quantile(mean_wt_blocks_interp(:)*100, 0.95);
    %alpha_caxis_subj_low  = quantile(mean_wt_blocks_interp(:)*100, 1-percentile_caxis);
    caxis([min(mean_wt_blocks_interp(:)*100), max(mean_wt_blocks_interp(:)*100)]);
    %clim('auto');

    yyaxis right
    plot((1:mean_length_blocks)/fs, mean_kin_subj, 'k', 'LineWidth', 2);
    ylabel('Median over blocks of median kinematic signal');

    % Save in PNG with high resolution
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 600]); 
    set(gca, 'FontSize', 12); 
    fig_name = sprintf('Averaged_TF_%s_all_blocks', subject_id);

    out_dir_TF = fullfile(base_path,'RESULTS_final', 'TF_pow_norm_phase_results_01');
    if ~exist(out_dir_TF, 'dir'), mkdir(out_dir_TF); end

    print(fullfile(out_dir_TF, fig_name), '-dpng', '-r300'); % -r300 = 300 dpi


    % % Plot baseline and other kinematic phases PSDs
    % figure('Name',['Baseline and Phases PSD - ' subject_id ' (All Blocks)'], 'NumberTitle','off');
    % hold on; grid on;
    % 
    % === Definition colours for phase ===
    col_rest     = [0.8 0.2 0.2];      
    col_reach    = [0.3 0.8 0.3];      
    col_grasp    = [0.2 0.6 1];   
    col_pull     = [0.6 0.1 0.9]; 


    % Plot all kinematic phases NORMALIZED PSDs 
    figure('Name',['Phases Normalized PSD - ' subject_id ' (All Blocks)'], 'NumberTitle','off');
    hold on; grid on;

    rest_psd_norm_subj     = all_med_phases_psd_norm_subj{1,s}(1:min_fr_glob);
    reach_psd_norm_subj    = all_med_phases_psd_norm_subj{2,s}(1:min_fr_glob);
    grasp_psd_norm_subj    = all_med_phases_psd_norm_subj{3,s}(1:min_fr_glob);
    pull_psd_norm_subj     = all_med_phases_psd_norm_subj{4,s}(1:min_fr_glob);

    % === Plot ===
    p1 = plot(subj_fr, rest_psd_norm_subj, 'Color', col_rest, 'LineWidth', 2.5, 'DisplayName', 'Rest (before/after)');
    p2 = plot(subj_fr, reach_psd_norm_subj, 'Color', col_reach, 'LineWidth', 2.5, 'DisplayName', 'Reach');
    p3 = plot(subj_fr, grasp_psd_norm_subj, 'Color', col_grasp, 'LineWidth', 2.5, 'DisplayName', 'Grasp');
    p4 = plot(subj_fr, pull_psd_norm_subj, 'Color', col_pull, 'LineWidth', 2.5, 'DisplayName', 'Pull');

    xlabel('Frequency (Hz)', 'FontSize', 12);
    ylabel('Normalized Power wrt Rest', 'FontSize', 12);
    xlim([min(subj_fr) 50])
    ylim([-0.8 0.8]);
    title(['Median Normalized PSD across Kinematic Phases - ' subject_id], 'FontSize', 14);
    legend([p1 p2 p3 p4], 'Location', 'best');
    set(gca, 'FontSize', 12);

    out_dir_norm = fullfile(base_path,'RESULTS_final', 'Baseline_01_median','NORM');
    if ~exist(out_dir_norm, 'dir'), mkdir(out_dir_norm); end

    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 6]);
    saveas(gcf, fullfile(out_dir_norm, ['PSD_norm_phases_' subject_id '.fig']));
    print(fullfile(out_dir_norm, ['PSD_norm_phases_' subject_id '.png']), '-dpng', '-r300');

    all_med_subj_baseline{s} = mean_subj_baseline;
    all_subj_fr{s} = subj_fr;

    % Save mean power across blocks
    all_subj_wt{s} =  mean_wt_blocks_interp;
    all_subj_wt_no_norm{s} =  mean_wt_blocks_interp_no_norm;
    all_subj_kin{s} = mean_kin_subj;

    % Save latencies mean 
    all_kin_latencies(1, s) = mean_subj_onset;
    all_kin_latencies(2, s) = mean_subj_grasp;
    all_kin_latencies(3, s) = mean_subj_pull;
    all_kin_latencies(4, s) = mean_subj_offset;
end

out_dir_struct = fullfile(base_path,'RESULTS_final', 'Baseline_01_median');
if ~exist(out_dir_struct, 'dir'), mkdir(out_dir_struct); end

save(fullfile(out_dir_struct, 'PSD_STRUCT_allSubjects.mat'), 'PSD_STRUCT', '-v7.3');

% --- PLOT MEAN BASELINE ACROSS SUBJECTS --- %

% Find minimum baseline across subjects
minumum_length_baseline = min(cellfun(@length, all_med_subj_baseline));

%Build average baseline across subjects
all_baseline_mat = zeros(minumum_length_baseline, n_subjects);
all_phases_psd_mat = zeros(4,minumum_length_baseline, n_subjects);
all_phases_psd_norm_mat = zeros(4,minumum_length_baseline, n_subjects);
for s = 1:n_subjects
    all_baseline_mat(:,s) = all_med_subj_baseline{s}(1:minumum_length_baseline);
    for p = 1:n_phases
        all_phases_psd_mat(p,:,s) = all_med_phases_psd_subj{p,s}(1:minumum_length_baseline);
        all_phases_psd_norm_mat(p,:,s) = all_med_phases_psd_norm_subj{p,s}(1:minumum_length_baseline);
    end
end

% Common freqs across subjects
common_fr = all_subj_fr{1}(1:minumum_length_baseline);

% Mean + SEM across subjects
grand_mean_baseline= mean(all_baseline_mat, 2);
sem_baseline  = std(all_baseline_mat, 0, 2) ./ sqrt(n_subjects);

mean_phases_psd = zeros(4,minumum_length_baseline);
mean_phases_psd_norm = zeros(4,minumum_length_baseline);
sem_phases_psd = zeros(4,minumum_length_baseline);
sem_phases_psd_norm = zeros(4,minumum_length_baseline);
for p = 1:n_phases
    mean_phases_psd(p,:) = mean(all_phases_psd_mat(p,:,:), 3);
    mean_phases_psd_norm(p,:) = mean(all_phases_psd_norm_mat(p,:,:), 3);
    sem_phases_psd(p,:) = std(all_phases_psd_mat(p,:,:),0, 3) ./ sqrt(n_subjects);
    sem_phases_psd_norm(p,:) = std(all_phases_psd_norm_mat(p,:,:), 0, 3) ./ sqrt(n_subjects);
end

% === Colours ===
col_baseline = [0.7 0.7 0.7];    
col_rest     = [0.8 0.2 0.2];      
col_reach    = [0.3 0.8 0.3];      
col_grasp    = [0.2 0.6 1];   
col_pull     = [0.6 0.1 0.9]; 
phase_colors = {col_rest, col_reach, col_grasp, col_pull};

%% =========================
% 1) PLOT PSD NOT NORMALIZED
% ============================
figure('Name','PSD Non normalizzata','Position',[100 100 1200 600]); hold on; grid on;

% Normalized phases
for p = 1:length(phase_names)
    mean_ph_n = squeeze(mean_phases_psd(p,:));
    sem_ph_n  = squeeze(sem_phases_psd(p,:));
    mean_ph_n = mean_ph_n(:)'; sem_ph_n = sem_ph_n(:)';

    common_fr = common_fr(:)';  

    fill([common_fr fliplr(common_fr)], ...
     [mean_ph_n+sem_ph_n fliplr(mean_ph_n-sem_ph_n)], ...
     phase_colors{p}, ...
     'EdgeColor','none', ...
     'FaceAlpha',0.25, ...
     'HandleVisibility','off');
    plot(common_fr, mean_ph_n, 'Color', phase_colors{p}, 'LineWidth', 2.5, 'DisplayName', [phase_names{p}]);
end

xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Power (a.u.)', 'FontSize', 12);
xlim([0 80]);
title('Grand-Average LFP Marginal Power Spectral Density', 'FontSize', 14);
legend('Location','best');
set(gca, 'FontSize', 12);

out_dir = fullfile(base_path,'RESULTS/Baseline_01_median');

% Salvataggio
print(fullfile(out_dir, 'Mean_PSD_all_subjects'), '-dpng', '-r300');
saveas(gcf, fullfile(out_dir, 'Mean_PSD_all_subjects.fig'));


%% =========================
% 2) PLOT PSD NORMALIZED
% ============================
figure('Name','PSD normalizzata','Position',[100 100 1200 600]); hold on; grid on;

% Normalized phases
for p = 1:length(phase_names)
    mean_ph_n = squeeze(mean_phases_psd_norm(p,:)*100);
    sem_ph_n  = squeeze(sem_phases_psd_norm(p,:)*100);
    mean_ph_n = mean_ph_n(:)'; sem_ph_n = sem_ph_n(:)';

    common_fr = common_fr(:)';  

    fill([common_fr fliplr(common_fr)], ...
     [mean_ph_n+sem_ph_n fliplr(mean_ph_n-sem_ph_n)], ...
     phase_colors{p}, ...
     'EdgeColor','none', ...
     'FaceAlpha',0.25, ...
     'HandleVisibility','off');
    plot(common_fr, mean_ph_n, 'Color', phase_colors{p}, 'LineWidth', 2.5, 'DisplayName', [phase_names{p}]);
end

xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Relative Power (w.r.t. Rest)(%)', 'FontSize', 12);
xlim([0 80]);
title('Grand-Average LFP Marginal Normalized Power (w.r.t. Rest)', 'FontSize', 14);
legend('Location','best');
set(gca, 'FontSize', 12);
ax = gca;

out_dir = fullfile(base_path,'RESULTS/Baseline_01_median');

% Save
print(fullfile(out_dir, 'Mean_PSD_norm_all_subjects'), '-dpng', '-r300');
saveas(gcf, fullfile(out_dir, 'Mean_PSD_norm_all_subjects.fig'));



%% --- PLOT AVERAGE NORMALIZED POWER TF --- %

% Median temporal length for interpolation
median_length_wt = round(median(cellfun(@(x) size(x,2), all_subj_wt)));
% Min length in frequency
min_freqs        = min(cellfun(@(x) size(x,1), all_subj_wt));

% Matrix 3D 
all_wt_mat  = zeros(min_freqs, median_length_wt, n_subjects);
all_kin_mat = zeros(median_length_wt, n_subjects);

for s = 1:n_subjects
    wt_subj  = all_subj_wt{s};
    kin_subj = all_subj_kin{s};

    % temporal original axis
    or_length  = size(wt_subj,2);
    t_orig = 1:or_length;

    % Target temporal axis
    t_target = linspace(1,or_length,median_length_wt);

    % Time interoplation
    wt_interp = interp1(t_orig, wt_subj(1:min_freqs,:)', t_target, 'linear', 'extrap')';
    kin_interp = interp1(t_orig, kin_subj(:), t_target, 'linear', 'extrap');

    % Save 3D matrices
    all_wt_mat(:,:,s)  = wt_interp;
    all_kin_mat(:,s)   = kin_interp;
end

% Mean Power and kinematic signal across subjects
all_subj_avg_wt  = mean(all_wt_mat, 3);
all_subj_avg_kin = mean(all_kin_mat, 2);

% SEM across subjects of Power TF matrix 
sem_wt = std(all_wt_mat, 0, 3) ./ sqrt(n_subjects);


time_axis = (1:median_length_wt) / fs;


freq_axis = fr_glob(1:min_freqs);

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
sem_wt_plot = sem_wt(freq_mask, :);

% --- Normalize kinematic signal to the 0-90 Hz frequency range ---
kin_min = min(freq_axis_plot);
kin_max = max(freq_axis_plot);
all_subj_avg_kin_norm = (all_subj_avg_kin - min(all_subj_avg_kin)) ...
                        / (max(all_subj_avg_kin) - min(all_subj_avg_kin)); % 0-1
all_subj_avg_kin_norm = all_subj_avg_kin_norm * (kin_max - kin_min) + kin_min;

% Plot Grand Average WT
figure('Name','Grand-Average LFP Time–Frequency Power (CWT, Rest-Normalized)','NumberTitle','off');
%pcolor(time_axis, freq_axis, all_subj_avg_wt);
pcolor(time_axis, freq_axis_plot, all_subj_avg_wt_plot*100);
shading interp; set(gca,'YScale','linear','FontSize',12);
xlabel('Time (s)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
xline(mean_onset_all,  'Label','Mov-ONSET',   'Color','w','LineWidth',2.5);
xline(mean_grasp_all,  'Label','Grasp Start', 'Color','w','LineWidth',2.5);
xline(mean_pull_all,   'Label','Grasp End',   'Color','w','LineWidth',2.5);
xline(mean_offset_all, 'Label','Mov-OFFSET',  'Color','w','LineWidth',2.5);
title('Grand-Average LFP Time–Frequency Power (CWT, Rest-Normalized)','FontSize',14);
c = colorbar; ylabel(c,'Relative Power (w.r.t. Rest) (%)');
colormap jet;
alpha_caxis_grand = quantile(all_subj_avg_wt(:)*100, 0.99);
%alpha_caxis_grand_low = quantile(all_subj_avg_wt(:)*100, 0.10);
caxis([min(all_subj_avg_wt(:)*100), max(all_subj_avg_wt(:)*100)]);
% clim('auto');


yyaxis right
%plot(time_axis, all_subj_avg_kin, 'k','LineWidth',2);
plot(time_axis, all_subj_avg_kin_norm, 'k','LineWidth',2);
ax = gca;
ax.YAxis(2).Visible = 'off';
%ylabel('Grand average kinematic signal');

% Save
set(gcf,'Units','pixels','Position',[100,100,1200,600]);
fig_name = sprintf('All_Subjects_Average_Power_TF');

out_dir_struct_avg = fullfile(base_path,'RESULTS_final', 'TF_pow_norm_phase_results_01');
if ~exist(out_dir_struct_avg, 'dir'), mkdir(out_dir_struct_avg); end

print(fullfile(out_dir_struct_avg,fig_name),'-dpng','-r300');
savefig(fullfile(out_dir_struct_avg,[fig_name '.fig']));



%% --- PLOT AVERAGE NOT NORMALIZED POWER TF --- %

% Median temporal length
median_length_wt = round(median(cellfun(@(x) size(x,2), all_subj_wt)));
% Lunghezza minima in frequenza
min_freqs        = min(cellfun(@(x) size(x,1), all_subj_wt));

% 3D Matrices
all_wt_mat  = zeros(min_freqs, median_length_wt, n_subjects);
all_kin_mat = zeros(median_length_wt, n_subjects);

for s = 1:n_subjects
    wt_subj  = all_subj_wt_no_norm{s};
    kin_subj = all_subj_kin{s};

    % Original temporal axis
    or_length  = size(wt_subj,2);
    t_orig = 1:or_length;

    % Target temporal axis
    t_target = linspace(1,or_length,median_length_wt);

    % Time interpolation
    wt_interp = interp1(t_orig, wt_subj(1:min_freqs,:)', t_target, 'linear', 'extrap')';
    kin_interp = interp1(t_orig, kin_subj(:), t_target, 'linear', 'extrap');

    % Save 3D Matrices
    all_wt_mat(:,:,s)  = wt_interp;
    all_kin_mat(:,s)   = kin_interp;
end

% Mean Power and KIN across subjects
all_subj_avg_wt  = mean(all_wt_mat, 3);
all_subj_avg_kin = mean(all_kin_mat, 2);

% SEM across subjects of Power TF matrix 
sem_wt = std(all_wt_mat, 0, 3) ./ sqrt(n_subjects);


time_axis = (1:median_length_wt) / fs;


freq_axis = fr_glob(1:min_freqs);

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
sem_wt_plot = sem_wt(freq_mask, :);

% --- Normalize kinematic signal to the 0-90 Hz frequency range ---
kin_min = min(freq_axis_plot);
kin_max = max(freq_axis_plot);
all_subj_avg_kin_norm = (all_subj_avg_kin - min(all_subj_avg_kin)) ...
                        / (max(all_subj_avg_kin) - min(all_subj_avg_kin)); % 0-1
all_subj_avg_kin_norm = all_subj_avg_kin_norm * (kin_max - kin_min) + kin_min;

% Plot Grand Average WT
figure('Name','Grand-Average LFP Time–Frequency Power (CWT, Not Normalized)','NumberTitle','off');
%pcolor(time_axis, freq_axis, all_subj_avg_wt);
pcolor(time_axis, freq_axis_plot, all_subj_avg_wt_plot*100);
shading interp; set(gca,'YScale','linear','FontSize',12);
xlabel('Time (s)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
xline(mean_onset_all,  'Label','Mov-ONSET',   'Color','w','LineWidth',2.5);
xline(mean_grasp_all,  'Label','Grasp Start', 'Color','w','LineWidth',2.5);
xline(mean_pull_all,   'Label','Grasp End',   'Color','w','LineWidth',2.5);
xline(mean_offset_all, 'Label','Mov-OFFSET',  'Color','w','LineWidth',2.5);
title('Grand-Average LFP Time–Frequency Power (CWT, Not Normalized) ','FontSize',14);
c = colorbar; ylabel(c,'Power (a.u.)');
colormap jet;
alpha_caxis_grand = quantile(all_subj_avg_wt(:)*100, 0.95);
%alpha_caxis_grand_low = quantile(all_subj_avg_wt(:)*100, 0.10);
caxis([min(all_subj_avg_wt(:)*100), alpha_caxis_grand]);
% clim('auto');


yyaxis right
%plot(time_axis, all_subj_avg_kin, 'k','LineWidth',2);
plot(time_axis, all_subj_avg_kin_norm, 'k','LineWidth',2);
ax = gca;
ax.YAxis(2).Visible = 'off';
%ylabel('Grand average kinematic signal');

% Save
set(gcf,'Units','pixels','Position',[100,100,1200,600]);
fig_name = sprintf('All_Subjects_Average_Power_TF');

out_dir_struct_avg_NO_NORM = fullfile(base_path,'RESULTS_final', 'TF_pow_norm_phase_results_NO_NORM');
if ~exist(out_dir_struct_avg_NO_NORM, 'dir'), mkdir(out_dir_struct_avg_NO_NORM); end

print(fullfile(out_dir_struct_avg_NO_NORM,fig_name),'-dpng','-r300');
savefig(fullfile(out_dir_struct_avg_NO_NORM,[fig_name '.fig']));

