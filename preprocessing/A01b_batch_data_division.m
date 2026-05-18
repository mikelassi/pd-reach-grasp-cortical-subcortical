%% ==============================================================
%  BATCH SCRIPT: Split BrainVision files (Input Drive -> Output Drive)
% ===============================================================

clear; clc;

% Define your list of subjects
subject_list = {'wue02', 'wue03', 'wue05', 'wue06', 'wue07', 'wue09', 'wue10', 'wue11'}; 

% --- DUAL DRIVE CONFIGURATION ---
% Where the RAW data lives
input_base_path = 'E:\OneDrive - Scuola Superiore Sant''Anna\Data\Parkinson_ReachGrasp\Raw\';

% Where the PROCESSED data will go
output_base_path = 'H:\Parkinson_ReachGrasp\Reprocessing\';

% Initialize EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

for s = 1:length(subject_list)
    
    subject_id = subject_list{s};
    subj_input_dir = fullfile(input_base_path, subject_id);
    
    fprintf('\n=== Processing subject %s ===\n', subject_id);
    
    % 1. Dynamically find the date folder containing 'Sync_KIN' on the input drive
    date_folders = dir(subj_input_dir);
    date_folders = date_folders([date_folders.isdir] & ~ismember({date_folders.name}, {'.', '..'}));
    
    sync_kin_path = '';
    for d = 1:length(date_folders)
        potential_path = fullfile(subj_input_dir, date_folders(d).name, 'Sync_KIN');
        if exist(potential_path, 'dir')
            sync_kin_path = potential_path;
            break;
        end
    end
    
    % Check if Sync_KIN was found
    if isempty(sync_kin_path)
        warning('No Sync_KIN folder found for subject %s on input drive. Skipping...', subject_id);
        continue;
    end
    
    % 2. Define and create the output folders on the H: drive
    save_dir_base = fullfile(output_base_path, subject_id, '01_Extracted');
    
    % Specific modality folders
    save_path_eeg = fullfile(save_dir_base, 'EEG');
    save_path_lfp = fullfile(save_dir_base, 'LFP');
    save_path_kin = fullfile(save_dir_base, 'EMG_KIN');
    
    % Create directories if they don't exist
    if ~exist(save_path_eeg, 'dir'), mkdir(save_path_eeg); end
    if ~exist(save_path_lfp, 'dir'), mkdir(save_path_lfp); end
    if ~exist(save_path_kin, 'dir'), mkdir(save_path_kin); end
    
    % 3. Find all .vhdr files in the Sync_KIN folder
    vhdr_files = dir(fullfile(sync_kin_path, '*.vhdr'));
    
    if isempty(vhdr_files)
        warning('No .vhdr files found in %s', sync_kin_path);
        continue;
    end
    
    % 4. Loop through each .vhdr file (Block)
    for f = 1:length(vhdr_files)
        
        filename_vhdr = vhdr_files(f).name;
        [~, name_base, ~] = fileparts(filename_vhdr);
        fprintf('  > Splitting file: %s\n', filename_vhdr);
        
        % Load the vhdr file
        EEG = pop_loadbv(sync_kin_path, filename_vhdr);
        
        % --- EEG channel extraction ---
        EEG_eeg = pop_select(EEG, 'channel', 1:126);
        EEG_eeg.setname = ['eeg_' name_base];
        
        % --- LFP channel extraction ---
        EEG_lfp = pop_select(EEG, 'channel', 127:128);
        EEG_lfp.setname = ['lfp_' name_base];
        
        % --- EMG-KINEMATIC channel extraction ---
        EEG_emgkin = pop_select(EEG, 'channel', 129:EEG.nbchan);
        EEG_emgkin.setname = ['kinematic_' name_base];
        
        % --- Save datasets to their specific target folders ---
        % Uncomment the EEG and LFP lines if you want to extract and save them as well
        pop_saveset(EEG_eeg, 'filename', [EEG_eeg.setname '.set'], 'filepath', save_path_eeg);
        pop_saveset(EEG_lfp, 'filename', [EEG_lfp.setname '.set'], 'filepath', save_path_lfp);
        pop_saveset(EEG_emgkin, 'filename', [EEG_emgkin.setname '.set'], 'filepath', save_path_kin);
        
        % Clear EEGLAB memory to prevent RAM crash across loops
        ALLEEG = [];
        CURRENTSET = 0;
    end
end

disp('Dual-drive batch splitting completed successfully!');