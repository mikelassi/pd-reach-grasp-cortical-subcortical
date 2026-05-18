
%% EEG PREPROCESSING PIPELINE (EEGLAB)
% ------------------------------------------------------------
% Project: Cortical–Subcortical Connectivity in Parkinson’s Disease
% Author: Tommaso Marcantoni
% Dependencies:
%   - EEGLAB (tested with v2025.0.0)
%   - clean_rawdata
%   - ZapLine-plus
%   - ICLabel
%   - DIPFIT
%
% Input:
%   - EEGLAB .set files with events:
%       TENStrigger, A_Tx, F_Tx
%
% Output:
%   - Preprocessed EEG .set files
%   - Subject-level preprocessing summary (.txt)
%
% Notes:
%   - Dataset paths are currently hard-coded
%   - Script provided as RAW reference implementation
% ------------------------------------------------------------

clear; clc;

subject_list = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'}; 

eeglab nogui;

%% PATHS (EDIT THIS SECTION)

base_path = 'C:\Users\tomma\OneDrive - University of Pisa\Desktop\TESI\Dataset_tesi';

eeglabroot = 'C:\Users\tomma\Downloads\eeglab_current\eeglab2025.0.0';
templateChanFile = fullfile(eeglabroot, 'plugins', 'dipfit', 'standard_BEM', 'elec', 'standard_1005.elc');
hdmFile = fullfile(eeglabroot, 'plugins', 'dipfit', 'standard_BEM', 'standard_vol.mat');
mriFile = fullfile(eeglabroot, 'plugins', 'dipfit', 'standard_BEM', 'standard_mri.mat');

for s = 1:length(subject_list)

    subject_id = subject_list{s};
    set_path = fullfile(base_path, subject_id, 'Extracted', 'EEG_wEv');
    
    % List of .set files in the folder
    set_files = dir(fullfile(set_path, '*.set'));

    % Handle corrupted block  --> subject wue06 
    if any(strcmp(subject_id, {'wue06'})) 
        set_files = set_files(1:2);      % Third block EEG is corrupted for wue06
    end
    
    % Struct for preprocessing informations for a subject
    prep_info = struct();
    
    for f = 1:length(set_files)
        
        % Store file name for naming preprocessed set
        set_file = set_files(f).name;
        fprintf('\nCarico file: %s\n', set_file);
    
        % Load .set EEG file
        EEG = pop_loadset('filename', set_file, 'filepath', set_path);
        [~, name_base] = fileparts(set_file);
        
        % Maintain the data in double precision
        pop_editoptions('option_single', false);
    
        % Manual division of signal between 50 sample after TENSTRIGGER  and FT10 + 300 samples
        % All events 
        event_types = {EEG.event.type};

        % Latenza di TENSTRIGGER (synchronization signal)
        lat_TEN_idx_all = find(strcmp(event_types, 'TENStrigger'));

        % Latency of A_T1
        AT1_latency_idx = find(strcmp(event_types, 'A_T1'));
        AT1_latency = EEG.event(AT1_latency_idx).latency;
        
        % First TENSTRIGGER
        lat_TEN_idx = lat_TEN_idx_all(1);
        lat_TEN = EEG.event(lat_TEN_idx).latency;

        % Find event F_Tx latency (first need to find x--> n° of trials)
        is_FT = startsWith(event_types, 'F_T');                                      % all events that starts with F_T
        nums = cellfun(@(x) str2double(extractAfter(x, 'F_T')), event_types(is_FT)); % extract num after F_T
        [max_num, idx_max] = max(nums);                                              % find maximum trial (eg. 10)
        last_event_label = ['F_T' num2str(max_num)];
        second_latency_idx = find(strcmp(event_types, ['F_T' num2str(max_num)]));
        second_latency = EEG.event(second_latency_idx).latency;                      % finally F_Tx latency
        
        start_point = lat_TEN + 50;  % TENStrigger latency + 50 samples
        end_point = second_latency + 300;
    
        % Mask for the plot --> 1 only in the EEG region that is kept
        manual_mask = false(1, EEG.pnts); 
        manual_mask((start_point):(end_point)) = 1;
    
        % === Select the region (pop_select) ===
        original_EEG = EEG;
        EEG = pop_select(EEG, 'point', [(start_point) (end_point)]);

        % Store proportion of removed data
        EEG.etc.clean_info.windowRejectionRate = (length(original_EEG.data) - length(EEG.data))/length(original_EEG.data);
    
        % Plot to visually check MANUALLY removed portions  
        figure;
        yyaxis left
        plot(original_EEG.data(1,:), 'b'); % signal
        ylabel('Amplitude (\muV)');
        xlabel('Sample index');
        title('EEG channel 1 with MANUAL clean mask and events');
        grid on;

        hold on;
        xline(AT1_latency, 'k--', 'A_T1', ...
              'LabelVerticalAlignment','bottom', ...
              'LabelHorizontalAlignment','right', ...
              'FontSize',6);

        xline(start_point, 'k--', 'Start_win', ...
              'LabelVerticalAlignment','bottom', ...
              'LabelHorizontalAlignment','right', ...
              'FontSize',6);

        xline(end_point, 'k--', 'End_win', ...
              'LabelVerticalAlignment','bottom', ...
              'LabelHorizontalAlignment','right', ...
              'FontSize',6);

        yyaxis right
        plot(manual_mask, 'r', 'LineWidth', 1.2); % mask
        ylabel('Clean sample mask');
        ylim([-0.1 1.1]);

        legend({'EEG signal', 'Clean mask'}, 'Location', 'best');

        event_types = {EEG.event.type};
        event_latencies = round([EEG.event.latency]); % sample indices

        % Handle corrupted trials for wue03 ( block 3 trial 4) and wue05
        % (block 1 trial 4)
        if (strcmp(subject_id, 'wue03') && f == 3)
            trial_to_remove = 4;
        elseif (strcmp(subject_id, 'wue05') && f == 1)
            trial_to_remove = 4;
        else
            trial_to_remove = [];
        end
        
        if ~isempty(trial_to_remove)
        
            % Identify trial region to remove
            F_idx_pre  = find(strcmp(event_types, sprintf('F_T%d', trial_to_remove - 1)));
            A_idx_curr = find(strcmp(event_types, sprintf('A_T%d', trial_to_remove)));
            F_idx_curr = find(strcmp(event_types, sprintf('F_T%d', trial_to_remove)));
            A_idx_next = find(strcmp(event_types, sprintf('A_T%d', trial_to_remove + 1)));
        
            % Mean point in the rest phase
            F_lat_pre  = event_latencies(F_idx_pre);
            A_lat_curr = event_latencies(A_idx_curr);
            F_lat_curr = event_latencies(F_idx_curr);
            A_lat_next = event_latencies(A_idx_next);
    
            mid_FA_pre  = round(F_lat_pre + (A_lat_curr - F_lat_pre)/2);
            mid_FA_post = round(F_lat_curr + (A_lat_next - F_lat_curr)/2);
    
            % Window removal
            EEG = pop_select(EEG, 'rmpoint', [mid_FA_pre mid_FA_post]);
        end
    
        % === HIGH PASS FILTER ===
        EEG_raw = EEG;
        EEG = pop_eegfiltnew(EEG, 'locutoff',1); % High pass filtering
    
        % % === REMOVAL OF LINE NOISE ===  Comparing cleanLine and ZapLine
        % the latter seems faster and more accurate in removing 50 Hz and
        % harmonics
    
        EEG  = clean_data_with_zapline_plus_eeglab_wrapper(EEG, struct('noisefreqs',50)); % By default even just with 50 it removes also the harmonics
        EEG_zap = EEG;
        
        % figure;
        % pop_spectopo(EEG_zap, 1, [], 'EEG');
        
    
        % === CLEAN RAW DATA === --> removed removal of windows. Channels are
        % identified as artifactual if they are less correlated than 0.6, flat
        % for more than 10 seconds or line noise power wrt signal power bigger
        % than 4 sd the mean among the other channels
        EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',10,'ChannelCriterion',0.6,'LineNoiseCriterion',4,'BurstCriterion','off','BurstRejection','off','WindowCriterion','off'); 
        EEG_clean = EEG;
        % figure;
        % pop_eegplot( EEG, 1, 1, 1);
    
        % Store channel rejection statistics [ EEG.etc.clean_info ]
        removed_channels = find(~EEG.etc.clean_channel_mask);
        EEG.etc.clean_info.num_interp_channels = sum(~EEG.etc.clean_channel_mask);          % Save number of channels to interpolate
        EEG.etc.clean_info.name_interp_channels = {EEG_zap.chanlocs(removed_channels).labels};
    
        % === INTERPOLATION OF SELECTED CHANNELS ===
        original_chanlocs = EEG_zap.chanlocs;
        EEG = pop_interp(EEG, original_chanlocs, 'spherical');
    
        % Print warning message if more than 30 channels are interpolated
        if EEG.etc.clean_info.num_interp_channels > 30
            fprintf('Number of interpolated channels larger than 30 : n° channels %d \n', EEG.etc.clean_info.num_interp_channels);
        end
    
        % == MAKOTO'S MEAN REFERENCING === --> add zero channel to avoid rank
        % deficiency in ICA
        EEG.nbchan = EEG.nbchan+1;
        EEG.data(end+1,:) = zeros(1, EEG.pnts);
        EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
        EEG = pop_reref(EEG, []);
        EEG = pop_select( EEG,'nochannel',{'initialReference'});
    
        % Compute data RANK --> this determine the number of ICs to find (
        % depends on the number of channels interpolated)
        dataCov = cov(double(EEG.data'));
        eigenVals = eig(dataCov);
        dataRank = sum(eigenVals > 1e-6); 
        fprintf('Data rank estimate: %d\n', dataRank);
    
        % === ICA with algorithm picard ===
        EEG = pop_runica(EEG, 'icatype', 'picard', 'maxiter', 1000, 'verbose', 'off', 'pca', dataRank); 
        %EEG = pop_runica(EEG, 'icatype','runica','extended',1,'pca',dataRank);

        % === DIPFIT ===
        eeglabroot = 'C:\Users\tomma\Downloads\eeglab_current\eeglab2025.0.0';

        % We are using digitized channel locations that have 10-20 system names (Fz, Cz, ...) 
        [~, coordinateTransformParameters] = coregister(EEG.chanlocs, 'C:\Users\tomma\Downloads\eeglab_current\eeglab2025.0.0/plugins/dipfit/standard_BEM/elec/standard_1005.elc', 'warp', 'auto', 'manual', 'off');

        EEG = pop_dipfit_settings(EEG, 'hdmfile', hdmFile, 'coordformat', 'MNI', ...
            'mrifile', mriFile, 'chanfile', templateChanFile, ...
            'coord_transform', coordinateTransformParameters, ...
            'chansel', 1:EEG.nbchan);

        % === DIPFIT === Compute dipole fitting
        EEG = pop_multifit(EEG, 1:EEG.nbchan, 'threshold', 100, 'dipplot', 'off');

        % Search for and estimate symmetrically constrained bilateral dipoles
        EEG = fitTwoDipoles(EEG, 'LRR', 35);

        % === ICLABEL CLASSIFICATION ===
        EEG = pop_iclabel(EEG, 'default');

        % % Visualize IC classification
        % disp('Classificazioni IC (ICLabel):');
        % disp(EEG.etc.ic_classification.ICLabel.classifications);

        % === Identification badICs based on ICLabel and dipfit_RV === 
        % EEG.etc.ic_classification.ICLabel.classifications(:,7) > 0.9 | ... % other
        badICs_ICLabel = find(EEG.etc.ic_classification.ICLabel.classifications(:,6) > 0.9 | ...                   % channel
            EEG.etc.ic_classification.ICLabel.classifications(:,5) > 0.9 | ...                   % line
            EEG.etc.ic_classification.ICLabel.classifications(:,4) > 0.9 | ...                   % heart
            EEG.etc.ic_classification.ICLabel.classifications(:,3) > 0.9 | ...                   % muscle
                      EEG.etc.ic_classification.ICLabel.classifications(:,2) > 0.9);             % eye

        badICs_dipfit = find([EEG.dipfit.model.rv] > 0.5);                                       % RV larger than 0.5

        badICs = union(badICs_ICLabel, badICs_dipfit);
        all_ICs = 1:size(EEG.icaweights,1);
        kept_ICs = setdiff(all_ICs, badICs);

        fprintf('Removal of %d artifact ICs \n', length(badICs));             % Print n° removed ICs

        % === Variance explained per IC ===
        ica_act = EEG.icaact;
        total_var = sum(var(EEG.data, 0, 2)); % Total variance
        ic_var   = var(ica_act, 0, 2);   % Single IC variance
        ic_explained = ic_var ./ total_var;   % Explained variance per IC
        EEG.ic_explained = ic_explained;                                      % Save explained variance by each IC in set [EEG.ic_explained]

        % % Topoplots of good and bad ICs
        % allICs = 1:size(EEG.icaweights,1);
        % goodICs = setdiff(allICs, badICs);  
        % 
        % % Visualize bad ICs 
        % figure('Name', sprintf('%s - Bad ICs', set_file), 'NumberTitle','off');
        % pop_topoplot(EEG, 0, badICs);
        % sgtitle('Bad ICs');  % Title for the figure
        % 
        % % Visualize good ICs
        % figure('Name', sprintf('%s - Good ICs', set_file), 'NumberTitle','off');
        % pop_topoplot(EEG, 0, goodICs(1:20));
        % sgtitle('Good ICs');
        % 
        % % Optional: also plot time course of first good IC
        % figure('Name', sprintf('%s - Time course IC %d', set_file, goodICs(1)), 'NumberTitle','off');
        % plot(EEG.icaact(goodICs(1), :));
        % xlabel('Sample');
        % ylabel('Amplitude (\muV)');
        % title(sprintf('Time course of IC %d', goodICs(1)));
        % grid on;

        % === REMOVAL badICs (pop_subcomp) ===
        EEG_preprocessed = pop_subcomp(EEG, badICs, 0);

        % Store in EEG_preprocessed
        EEG_preprocessed.kept_ICs = kept_ICs;
        EEG_preprocessed.removed_ICs = badICs;
        EEG_preprocessed.all_ICs = all_ICs;
        EEG_preprocessed.all_ICs_iclabel = EEG.etc.ic_classification.ICLabel.classifications(:, :);
       
        % === SAVE Preprocessed set ===
        out_path = fullfile(base_path, subject_id, 'Preprocessed_01', 'EEG');
        pop_saveset(EEG_preprocessed, 'filename', set_file, 'filepath', out_path);
        
        % === FILL STRUCT FOR PREPROCESSING INFORMATIONS ===
        prep_info(f).filename = set_file;
        prep_info(f).num_interp_channels = EEG_preprocessed.etc.clean_info.num_interp_channels;
        prep_info(f).interp_channel_names = EEG_preprocessed.etc.clean_info.name_interp_channels;
        prep_info(f).num_removed_ICs = length(badICs);
        prep_info(f).removed_ICs = badICs;
        prep_info(f).removed_ICs_iclabel = EEG.etc.ic_classification.ICLabel.classifications(badICs, :);
        prep_info(f).num_kept_ICs = length(kept_ICs);
        prep_info(f).kept_ICs = kept_ICs;
        prep_info(f).kept_ICs_iclabel = EEG.etc.ic_classification.ICLabel.classifications(kept_ICs, :);
        prep_info(f).ICs_dipfit_RV = [EEG.dipfit.model.rv];
        prep_info(f).ic_explained = ic_explained; % Explained variance of all ICs
       

        % Store in ALLEEG for convenience
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG_preprocessed, f);

    end
    
    % Print a SUMMARY of preprocessing statistics per subject
    fprintf('\n=== SUMMARY per subject %s ===\n', subject_id);
    for f = 1:length(prep_info)
        fprintf('\nFile: %s\n', prep_info(f).filename);
        fprintf('  Interpolated channels (%d): %s\n', prep_info(f).num_interp_channels, strjoin(prep_info(f).interp_channel_names, ', '));
        fprintf('  Removed ICs: %d (%s)\n', prep_info(f).num_removed_ICs, mat2str(prep_info(f).removed_ICs));
        fprintf('  Kept ICs: %d (%s)\n', prep_info(f).num_kept_ICs, mat2str(prep_info(f).kept_ICs));
        fprintf('  Mean dipole RV: %.3f\n', mean(prep_info(f).ICs_dipfit_RV));
        fprintf('  Total variance explained by removed ICs: %.2f%%\n', sum(prep_info(f).ic_explained(prep_info(f).removed_ICs)) * 100);
        fprintf('  Total variance explained by kept ICs: %.2f%%\n', sum(prep_info(f).ic_explained(setdiff(1:length(prep_info(f).ic_explained), prep_info(f).removed_ICs))) * 100);
    end

    % Save to a .TXT file 
    outfile = fullfile(out_path, sprintf('summary_%s.txt', subject_id));
    fid = fopen(outfile, 'w');   

    % === Intro note: where to find all preprocessing statistics ===
    fprintf(fid, '=============================\n');
    fprintf(fid, ' PREPROCESSING SUMMARY FILE\n');
    fprintf(fid, ' Subject: %s\n', subject_id);
    fprintf(fid, '=============================\n\n');
    fprintf(fid, 'Note: This file contains a short summary.\n');
    fprintf(fid, 'For complete details please consult:\n');
    fprintf(fid, ' - EEG_preprocessed.etc.clean_info  -> info on clean_rawdata (proportion of removed data, number and name of interpolated channels)\n');
    fprintf(fid, ' - EEG_preprocessed.etc.ic_classification.ICLabel.classifications ->  ICLabel probabilities of kept ICs\n');
    fprintf(fid, ' - EEG_preprocessed.removed_ICs_iclabel ->  ICLabel probabilities of removed ICs\n');
    fprintf(fid, ' - EEG_preprocessed.kept_ICs_iclabel ->  ICLabel probabilities of kept ICs\n');
    fprintf(fid, ' - EEG_preprocessed.dipfit.model  -> Parametes of dipoles for kept ICs \n');
    fprintf(fid, ' - EEG_preprocessed.bad_ICs_dipfit_model  -> Parametes of dipoles for removed ICs \n');
    fprintf(fid, ' - EEG_preprocessed.ic_explained -> Explained variance of all the ICs\n\n');

    fprintf(fid, '\n=== SUMMARY per subject %s ===\n', subject_id);
    for f = 1:length(prep_info)
        fprintf(fid, '\nFile: %s\n', prep_info(f).filename);
        fprintf(fid, '  Interpolated channels (%d): %s\n', prep_info(f).num_interp_channels, strjoin(prep_info(f).interp_channel_names, ', '));
        fprintf(fid, '  Removed ICs: %d (%s)\n', prep_info(f).num_removed_ICs, mat2str(prep_info(f).removed_ICs));
        fprintf(fid, '  Kept ICs: %d (%s)\n', prep_info(f).num_kept_ICs, mat2str(prep_info(f).kept_ICs));
        fprintf(fid, '  Mean dipole RV: %.3f\n', mean(prep_info(f).ICs_dipfit_RV));
        fprintf(fid, '  Total variance explained by removed ICs: %.2f%%\n', sum(prep_info(f).ic_explained(prep_info(f).removed_ICs)) * 100);
        fprintf(fid, '  Total variance explained by kept ICs: %.2f%%\n', sum(prep_info(f).ic_explained(setdiff(1:length(prep_info(f).ic_explained), prep_info(f).removed_ICs))) * 100);

        % === ICLabel classification of removed ICs ===
        if isfield(prep_info(f), 'removed_ICs_iclabel') && ~isempty(prep_info(f).removed_ICs_iclabel)
            fprintf(fid, '\n  ICLabel classification (Removed ICs):\n');

            % Add classes names
            iclabel_classes = {'Brain','Muscle','Eye','Heart','Line Noise','Channel Noise','Other'};
            fprintf(fid, '    %8s', 'IC#');
            for c = 1:length(iclabel_classes)
                fprintf(fid, '  %10s', iclabel_classes{c});
            end
            fprintf(fid, '\n');

            % Print ICLabel info to txt
            for r = 1:size(prep_info(f).removed_ICs_iclabel,1)
                ic_num = prep_info(f).removed_ICs(r);
                fprintf(fid, '    %8d', ic_num);
                fprintf(fid, '  %10.2f', prep_info(f).removed_ICs_iclabel(r,:) * 100); 
                fprintf(fid, '\n');
            end
        end
        if isfield(prep_info(f), 'kept_ICs_iclabel') && ~isempty(prep_info(f).kept_ICs_iclabel)
            fprintf(fid, '\n  ICLabel classification (Kept ICs):\n');
        
            % Add classes names
            iclabel_classes = {'Brain','Muscle','Eye','Heart','Line Noise','Channel Noise','Other'};
            fprintf(fid, '    %8s', 'IC#');
            for c = 1:length(iclabel_classes)
                fprintf(fid, '  %10s', iclabel_classes{c});
            end
            fprintf(fid, '\n');
        
            % Print ICLabel info to txt
            for r = 1:size(prep_info(f).kept_ICs_iclabel,1)
                ic_num = prep_info(f).kept_ICs(r);
                fprintf(fid, '    %8d', ic_num);
                fprintf(fid, '  %10.2f', prep_info(f).kept_ICs_iclabel(r,:) * 100); 
                fprintf(fid, '\n');
            end
        end
    end

    fclose(fid);
end

eeglab redraw;
