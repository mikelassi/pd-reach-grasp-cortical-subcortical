%% ==============================================================
%  Script for splitting a BrainVision EEG file (.vhdr)
%  into three separate datasets: EEG, LFP, and EMG/Kinematic.
%
%  Workflow:
%   1. Loads the .vhdr file from the specified path.
%   2. Extracts and saves separate datasets:
%        - EEG  : channels 1–126
%        - LFP  : channels 127–128
%        - EMG/KINEMATIC : channels 129–EEG.nbchan
%   3. Saves the datasets in .set format in the destination folder.
%
%  Requirements:
%   - EEGLAB installed and added to the MATLAB path
%   - BrainVision files (.vhdr) and corresponding .eeg/.vmrk
%     located in the specified folder
%
%  Output:
%   - One .set file for each channel group, ready for event
%     addition and preprocessing.
% ===============================================================

folder_path = 'C:\Users\tomma\OneDrive - University of Pisa\Desktop\TESI\Dataset_tesi\wue11\Sync_KIN\';
filename_vhdr = 'reach2grasp_pcs-normal-off_off-wue11-20141028-EegPcsEmgKin-3.vhdr';

% File name --> used to rename datasets
[~, name_base, ~] = fileparts(filename_vhdr);

% File loading (vhdr file)
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG = pop_loadbv(folder_path, filename_vhdr);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 1);

% EEG channel extraction
EEG_eeg = pop_select(EEG, 'channel', 1:126);
EEG_eeg.setname = ['eeg_' name_base];
[ALLEEG, EEG_eeg, ~] = eeg_store(ALLEEG, EEG_eeg, length(ALLEEG)+1);

% LFP channel extraction
EEG_lfp = pop_select(EEG, 'channel', 127:128);
EEG_lfp.setname = ['lfp_' name_base];
[ALLEEG, EEG_lfp, ~] = eeg_store(ALLEEG, EEG_lfp, length(ALLEEG)+1);

% EMG-KINEMATIC channel extraction
%EEG_emgkin = pop_select(EEG, 'channel', 129:154); 
EEG_emgkin = pop_select(EEG, 'channel', 129:EEG.nbchan);
EEG_emgkin.setname = ['kinematic_' name_base];
[ALLEEG, EEG_emgkin, ~] = eeg_store(ALLEEG, EEG_emgkin, length(ALLEEG)+1);

% Save datasets with appropriate names
save_path = 'C:\Users\tomma\OneDrive - University of Pisa\Desktop\TESI\Dataset_tesi\wue11\Extracted\'; 
% pop_saveset(EEG_eeg, 'filename', [EEG_eeg.setname '.set'], 'filepath', save_path);
% pop_saveset(EEG_lfp, 'filename', [EEG_lfp.setname '.set'], 'filepath', save_path);
pop_saveset(EEG_emgkin, 'filename', [EEG_emgkin.setname '.set'], 'filepath', save_path);

disp('Splitting completed successfully!');
