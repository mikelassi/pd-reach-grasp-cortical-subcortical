%% =======================================================
% EEG TOPOPLOT VISUALIZATION FROM PRECOMPUTED PSD STRUCTURES
%
% This script generates topographic maps (topoplots) of EEG normalized
% power across multiple frequency bands (Theta, Alpha, Beta, Gamma, Whole Beta)
% and behavioral phases (Rest, Reach, Grasp, Pull) using precomputed
% trial-wise PSD data stored in EEG_PSD_ALL_CONNECTIVITY.
%
% Key functionalities:
%   1. Load preprocessed EEG PSD structures for all subjects.
%   2. Compute mean and median normalized power per frequency band,
%      optionally mirroring left-handed subjects.
%   3. Plot individual subject topoplots and grand-average topoplots
%      across subjects for each frequency band and phase.
%   4. Generate animated topoplots illustrating power changes across phases.
%   5. Identify spatial clusters of desynchronization in the Reach+Grasp
%      phases for Whole Beta using percentile-based thresholds.
%   6. Save all figures and cluster information for further analysis.
%
% Outputs:
%   - Subject-specific topoplots (static and animated) per band and phase.
%   - Grand-average topoplots across subjects.
%   - Cluster masks for desynchronized channels (Reach+Grasp, Whole Beta).
%   - Saved .mat files and images compatible with subsequent analyses.
%
% Dependencies:
%   - EEGLAB (for EEG data structures and topoplot function)
%   - Precomputed PSD structures (EEG_PSD_ALL_CONNECTIVITY)
%
% Author: Tommaso Marcantoni
%% =======================================================

close all;

% === Percorso file ===
base_path    = 'C:\Users\tomma\OneDrive - University of Pisa\Desktop\TESI\Dataset_tesi';
data_file = fullfile(base_path, 'RESULTS/Topoplot_NORM_median_01/Mat_data', 'All_Data_TF.mat');
data_file_conn = fullfile(base_path, 'RESULTS/Topoplot_NORM_median_01/Mat_data', 'All_Data_TF_CONNECTIVITY.mat');
base_save_path = fullfile(base_path, 'RESULTS_final/Topoplot_NORM_median_post_avg_01');

save_dir_conv  = fullfile(base_save_path, 'High_Low_Beta');
save_dir_atg= fullfile(base_save_path, 'ATG');
save_dir_anim= fullfile(base_save_path, 'Animations');
save_dir_conv_whole  = fullfile(base_save_path, 'Whole_Beta');


% All subjects
base_save_path_all = fullfile(base_save_path, 'All_subj');
save_dir_conv_all  = fullfile(base_save_path_all, 'All_subj\High_Low_Beta');
save_dir_atg_all= fullfile(base_save_path_all, 'ATG');
save_dir_anim_all= fullfile(base_save_path_all, 'Animations');
save_dir_conv_whole_all = fullfile(base_save_path_all, 'Whole_Beta');

preproc_path = 'C:\Users\tomma\OneDrive - University of Pisa\Desktop\TESI\Dataset_tesi\wue02\Preprocessed_New';
file_EEG_ref = fullfile(preproc_path, 'EEG/', 'eeg_reach2grasp_pcs-normal-off_off-wue02-20150822-EegPcsEmgKin-1_wEv.set');

eeglab;
EEG = pop_loadset(file_EEG_ref);  % File needed only for chanlocs and sampling rate

fprintf('Loading EEG_PSD_ALL_CONNECTIVITY ...\n');

load(data_file_conn, 'EEG_PSD_ALL_CONNECTIVITY', '-mat');


n_subjects = numel(EEG_PSD_ALL_CONNECTIVITY);
fprintf('Subjects found: %d\n', n_subjects);

% Construct saving directories
all_save_dirs = { save_dir_conv,save_dir_conv_whole, save_dir_atg, save_dir_anim, ...
                  save_dir_conv_all, save_dir_conv_whole_all, save_dir_atg_all, save_dir_anim_all };

all_save_dirs{end+1} = fullfile(base_path, 'RESULTS_final', 'Topoplot_NORM_median_post_avg_01');

for d = 1:numel(all_save_dirs)
    dirpath = all_save_dirs{d};
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
        fprintf('Created folder: %s\n', dirpath);
    end
end

%% Function to assign hemisphere and controlateral label
function [hemi, contra_label] = assign_hemisphere(labels)
    n = numel(labels);
    hemi = cell(1,n);
    contra_label = cell(1,n);  % inizializziamo array per canali controlaterali
    
    % Costruisco un mapping base L <-> R
    for i = 1:n
        label = labels{i};
        
        % Check centrale
        if contains(label,'z','IgnoreCase',true)
            hemi{i} = 'C';
            contra_label{i} = '';  % nessun controlaterale
        else
            % Trova il primo numero nel nome
            nums = regexp(label,'\d+','match');
            if ~isempty(nums)
                num_val = str2double(nums{1});
                if mod(num_val,2) == 0
                    hemi{i} = 'R';
                else
                    hemi{i} = 'L';
                end
                
                % Trova il controlaterale
                % cambio il numero da pari a dispari o viceversa
                if mod(num_val,2) == 0
                    contra_num = num_val - 1;
                else
                    contra_num = num_val + 1;
                end
                
                % Sostituisco il numero nella stringa
                contra_label{i} = regexprep(label, '\d+', num2str(contra_num), 'once');
                
                % Se il controlaterale non è tra le labels esistenti, metto stringa vuota
                if ~ismember(contra_label{i}, labels)
                    contra_label{i} = '';
                end
            else
                hemi{i} = '';
                contra_label{i} = '';
            end
        end
    end
end
% %% Test funzione
% [hemi,contra_label] = assign_hemisphere({EEG.chanlocs.labels});


%% === Parametri ===
band_names = {'LowBeta','HighBeta','Theta','Alpha','Gamma','WholeBeta'};
band_regions = {[13 20],[20 30],[4 8],[8 13],[30 80], [13 30]};
phase_names = {'Rest','Reach','Grasp','Pull'};

% === Lunghezze delle frequenze ===
min_len = min(arrayfun(@(x) length(x.freq_axis), EEG_PSD_ALL_CONNECTIVITY));
[max_len, idx_max] = max(arrayfun(@(x) length(x.freq_axis), EEG_PSD_ALL_CONNECTIVITY));
freq_axis = EEG_PSD_ALL_CONNECTIVITY(idx_max).freq_axis(1:min_len);

n_channels = size(EEG_PSD_ALL_CONNECTIVITY(1).PSD, 3);
n_bands = numel(band_names);
n_phases = numel(phase_names);
fs = EEG.srate; 

% === Preallocazione ===
perc_band_power = zeros(n_phases, n_channels, n_bands, n_subjects);
mean_band_power_norm = zeros(n_phases, n_channels, n_bands, n_subjects);

%% === Ottieni labels e hemisferi + controlaterali ===
chan_labels = {EEG.chanlocs.labels};
[hemi_labels, contra_labels] = assign_hemisphere(chan_labels);

%% === Loop per soggetto con specchiamento per mancini ===
idx_left = [3,6];
for s = 1:n_subjects
    subj_id = EEG_PSD_ALL_CONNECTIVITY(s).subject_id;
    fprintf('\nProcessing subject: %s\n', subj_id);
    
    PSD_norm_all = EEG_PSD_ALL_CONNECTIVITY(s).PSD_norm_all*100; 
    freqs = EEG_PSD_ALL_CONNECTIVITY(s).freq_axis(:)';

    for p = 1:n_phases
       
        psd_p = squeeze(PSD_norm_all(p, :, :, :));

        for b = 1:n_bands
            band_range = band_regions{b};
            freq_idx = freqs >= band_range(1) & freqs <= band_range(2);

            band_mean_trials = squeeze(mean(psd_p(freq_idx, :, :), 1));
            
            band_median_norm = median(band_mean_trials, 2)';

            % Se soggetto mancino, specchio i valori L <-> R
            if ismember(s, idx_left)  % soggetto mancino
                % Trova coppie L-R uniche
                for ch = 1:n_channels
                    if strcmp(hemi_labels{ch}, 'L') && ~isempty(contra_labels{ch})
                        contra_idx = find(strcmp(contra_labels{ch}, chan_labels));
                        if ~isempty(contra_idx)
                            % Swap solo se il canale di controllo è R
                            if strcmp(hemi_labels{contra_idx}, 'R')
                                temp = band_median_norm(ch);
                                band_median_norm(ch) = band_median_norm(contra_idx);
                                band_median_norm(contra_idx) = temp;
                            end
                        end
                    end
                end
            end
            
            mean_band_power_norm(p, :, b, s) = band_median_norm;

        end
    end
end

% === Media tra soggetti ===

avg_mean_power_norm_band = squeeze(mean(mean_band_power_norm, 4));


%% ================================================================
%  🔹 1) TOPOPLOT PER OGNI SOGGETTO 
% ================================================================
plot_individuals = true; 

if plot_individuals
    for s = 1:n_subjects
        subj_id = EEG_PSD_ALL_CONNECTIVITY(s).subject_id;
        fprintf('Plotting subject: %s\n', subj_id);
        
        for b = 1:n_bands
            % Calcolo min e max su tutte le fasi per lo stesso soggetto e banda
            band_data = squeeze(mean_band_power_norm(:, :, b, s)); % [phase x channels]
            min_val = min(band_data, [], 'all');
            max_val = max(band_data, [], 'all');
            
            f = figure('Name', sprintf('Subject %s - %s', subj_id, band_names{b}), 'Color', 'w');
            
            % ciclo sulle fasi
            for p = 1:n_phases
                subplot(2, ceil(n_phases/2), p); 
                topoplot(squeeze(mean_band_power_norm(p, :, b, s)), EEG.chanlocs, ...
                         'electrodes', 'off', ...
                         'maplimits', [min_val max_val], ... % stessa scala per tutte le fasi
                         'style', 'fill');
                colorbar;

                % 🔻 Nome fase in basso
                xlim([-0.6 0.6]); ylim([-0.6 0.6]);
                text(0, -0.7, phase_names{p}, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'top', ...
                    'FontSize', 11, 'FontWeight', 'bold');
            end
            
            range_str = sprintf('%.1f–%.1f', band_regions{b}(1), band_regions{b}(2));  
            sgtitle(sprintf('EEG Mean Power in band - %s - %s [%s]Hz', subj_id, band_names{b}, range_str), ...
                'FontSize', 11, 'FontWeight', 'bold');
            
            % === Salvataggio figura ===
            if ismember(band_names{b}, {'LowBeta','HighBeta'})
                save_name = sprintf('Topoplot_%s_%s.png', subj_id, band_names{b});
                exportgraphics(f, fullfile(save_dir_conv, save_name), 'Resolution', 300);
                close(f);
            elseif ismember(band_names{b}, {'LowMI','HighMI'})
                save_name = sprintf('Topoplot_%s_%s.png', subj_id, band_names{b});
                exportgraphics(f, fullfile(save_dir_MI, save_name), 'Resolution', 300);
                close(f);
            elseif ismember(band_names{b}, {'Theta','Alpha','Gamma'})
                save_name = sprintf('Topoplot_%s_%s.png', subj_id, band_names{b});
                exportgraphics(f, fullfile(save_dir_atg, save_name), 'Resolution', 300);
                close(f);
            elseif ismember(band_names{b}, {'LowRange','HighRange'})
                save_name = sprintf('Topoplot_%s_%s.png', subj_id, band_names{b});
                exportgraphics(f, fullfile(save_dir_Range, save_name), 'Resolution', 300);
                close(f);
            elseif ismember(band_names{b}, {'WholeBeta'})
                save_name = sprintf('Topoplot_%s_%s.png', subj_id, band_names{b});
                exportgraphics(f, fullfile(save_dir_conv_whole, save_name), 'Resolution', 300);
                close(f);
            end

        end
    end
end
%% ================================================================
%  🔹 1) TOPOPLOT MEDIO TRA SOGGETTI
% ================================================================

bands_to_plot = {'LowBeta','HighBeta','Theta', 'Gamma', 'Alpha', 'WholeBeta'};

for b = 1:n_bands
    if ~ismember(band_names{b}, bands_to_plot)
        continue;
    end

    band_data_all = squeeze(avg_mean_power_norm_band(:, :, b)); % [phase x channels]
    clim = [min(band_data_all(:)) max(band_data_all(:))];

    f = figure('Name', sprintf('Average Topoplot - %s', band_names{b}), 'Color', 'w');

    for p = 1:n_phases
        subplot(2, ceil(n_phases/2), p);
        topoplot(squeeze(band_data_all(p, :)), EEG.chanlocs, ...
                 'electrodes', 'off', 'maplimits', clim, 'style', 'fill');
        colorbar;

        xlim([-0.6 0.6]); ylim([-0.6 0.6]);
        text(0, -0.75, phase_names{p}, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', ...
             'FontSize', 12, 'FontWeight', 'bold');
  
    end
    
    title('');   
    sgtitle(sprintf(' Grand-Average EEG Normalized Power (w.r.t. Rest) - %s [%d-%d] Hz ', band_names{b}, band_regions{b}(1), band_regions{b}(2)), ...
        'FontSize', 11, 'FontWeight', 'bold');

    % Salvataggio figura
    save_name = sprintf('Average_Topoplot_%s.png', band_names{b});
    exportgraphics(f, fullfile(save_dir_conv_all, save_name), 'Resolution', 300);
    close(f);
end


%% ================================================================
% 🔹 ANIMAZIONE TOPOPLOT PER TUTTI I SOGGETTI
% ================================================================

bands_to_animate = {'LowBeta','HighBeta', 'WholeBeta'}; % bande di interesse

for s = 1:n_subjects
    subj_id = EEG_PSD_ALL(s).subject_id;
    
    for b = 1:length(bands_to_animate)
        if ~ismember(band_names{b}, bands_to_animate)
            continue;
        end
        
        fprintf('Creating animated topoplot for subject %s - %s...\n', subj_id, band_names{b});
        f_anim = figure('Name', sprintf('%s - Animated Topoplot - %s', subj_id, band_names{b}), 'Color', 'w');
        colormap jet;
        
        data_all = squeeze(mean_band_power_norm(:, :, b, s));  % [phase x channels]
        clim = [min(data_all(:)) max(data_all(:))];  % limiti colore comuni
        
        im = cell(1, n_phases);
        
        for p = 1:n_phases
            cla;
            topoplot(squeeze(mean_band_power_norm(p, :, b, s)), EEG.chanlocs, ...
                     'electrodes', 'on', 'maplimits', clim, 'style', 'fill');
            colorbar;
           
            title(phase_names{p}, 'FontSize', 14, 'FontWeight', 'bold');
            

            annotation('textbox', [0, 0.02, 1, 0.05], ...
                'String', sprintf('Band: %s - Subject: %s', bands_to_animate{b}, subj_id), ...
                'FontSize', 12, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'EdgeColor', 'none');
            drawnow;
            pause(0.8);
            
            % Salva frame per GIF
            frame = getframe(f_anim);
            im{p} = frame2im(frame);
        end
        
        % Salvataggio GIF
        gif_name = fullfile(save_dir_anim, sprintf('Animated_Topoplot_%s_%s.gif', band_names{b}, subj_id));
        
        for p = 1:length(im)
            [A,map] = rgb2ind(im{p},256);
            if p == 1
                imwrite(A,map,gif_name,'gif','LoopCount',Inf,'DelayTime',0.8);
            else
                imwrite(A,map,gif_name,'gif','WriteMode','append','DelayTime',0.8);
            end
        end
        
        close(f_anim);
    end
end

%% ================================================================
% 🔹 ANIMAZIONE TOPOPLOT MEDIA TRA SOGGETTI
% ================================================================

bands_to_animate_all = {'LowBeta','HighBeta', 'Alpha', 'Theta', 'Gamma', 'WholeBeta'};  % bande di interesse

for b = 1:n_bands
    if ~ismember(band_names{b}, bands_to_animate_all)
        continue;
    end

    fprintf('Creating animated average topoplot for band %s...\n', band_names{b});

    f_anim_all = figure('Name', sprintf('Animated Average Topoplot - %s', band_names{b}), 'Color', 'w');
    colormap jet;

    data_all = squeeze(avg_mean_power_norm_band(:, :, b));  % [phase x channels]
    clim = [min(data_all(:)) max(data_all(:))];             % limiti colore comuni tra fasi

    im = cell(1, n_phases);

    for p = 1:n_phases
        cla;
        topoplot(squeeze(avg_mean_power_norm_band(p, :, b)), EEG.chanlocs, ...
                 'electrodes', 'on', 'maplimits', clim, 'style', 'fill', 'whitebk', 'on');
        colorbar;
        title(phase_names{p}, 'FontSize', 14, 'FontWeight', 'bold');

        annotation('textbox', [0, 0.02, 1, 0.05], ...
            'String', sprintf('Band: %s - Group Average', band_names{b}), ...
            'FontSize', 12, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'EdgeColor', 'none');
        drawnow;
        pause(0.8);

        % Salva frame per GIF
        frame = getframe(f_anim_all);
        im{p} = frame2im(frame);
    end

    % === Salvataggio GIF ===
    gif_name = fullfile(save_dir_anim_all, sprintf('Animated_Average_Topoplot_%s.gif', band_names{b}));

    for p = 1:length(im)
        [A,map] = rgb2ind(im{p},256);
        if p == 1
            imwrite(A,map,gif_name,'gif','LoopCount',Inf,'DelayTime',0.8);
        else
            imwrite(A,map,gif_name,'gif','WriteMode','append','DelayTime',0.8);
        end
    end

    close(f_anim_all);
end

%% ================================================================
% 🔹 CLUSTER IDENTIFICATION IN REACH+GRASP PHASE (Whole Beta only)
% Identify spatial clusters of desynchronization based on the mean
% normalized power across subjects, averaging the Reach and Grasp phases.
% ================================================================

fprintf('\n=== Identifying desynchronization clusters (Reach+Grasp, Whole Beta) ===\n');

% --- Parameters
phase_pair = {'Reach','Grasp'};
phase_idx = find(ismember(phase_names, phase_pair));

band_wholebeta_idx = find(strcmp(band_names, 'WholeBeta'));
radius_thresh = 0.40;  % radius (in normalized 2D scalp coordinates)

% --- Extract group mean data for Reach and Grasp, then average
data_reach = squeeze(avg_mean_power_norm_band(phase_idx(1), :, band_wholebeta_idx));
data_grasp = squeeze(avg_mean_power_norm_band(phase_idx(2), :, band_wholebeta_idx));
data_mean_RG = mean([data_reach; data_grasp], 1);  % average across phases

chanlocs = EEG.chanlocs;
xyz = [[chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]'];

% === Assign hemispheres using helper function ===
chan_labels = {EEG.chanlocs.labels};
[hemi_labels, ~] = assign_hemisphere(chan_labels);

% --- Separate indices per hemisphere
left_idx  = find(strcmp(hemi_labels, 'L'));
right_idx = find(strcmp(hemi_labels, 'R'));



%% === Cluster based on percentile threshold instead of spatial radius ===

fprintf('\n=== Identifying desynchronization clusters using Lower 5th percentile ===\n');

% Compute 5th percentile
perc5 = prctile(data_mean_RG, 5);

% Logical mask: channels below threshold
mask_desync = data_mean_RG <= perc5;

% Split into hemispheres
cluster_left_wholebeta  = left_idx(mask_desync(left_idx));
cluster_right_wholebeta = right_idx(mask_desync(right_idx));

fprintf('Whole Beta clusters identified by percentile: Left (%d ch), Right (%d ch)\n', ...
    numel(cluster_left_wholebeta), numel(cluster_right_wholebeta));

fprintf('\nLeft cluster channels:\n');
fprintf('%s ', chan_labels{cluster_left_wholebeta});
fprintf('\n\nRight cluster channels:\n');
fprintf('%s ', chan_labels{cluster_right_wholebeta});
fprintf('\n');

% Save also the threshold used
cluster_info = struct();
cluster_info.WholeBeta.Left  = cluster_left_wholebeta;
cluster_info.WholeBeta.Right = cluster_right_wholebeta;
cluster_info.Percentile_Threshold = perc5;
cluster_info.Percentile = 5;
cluster_info.Phases = phase_pair;

save(fullfile(base_path, 'RESULTS', 'Topoplot_NORM_median_01', ...
    'Clusters_ReachGrasp_WholeBeta_PERCENTILE.mat'), 'cluster_info');

fprintf('Cluster indices saved successfully (Reach+Grasp, Whole Beta).\n');

%% ================================================================
% 🔹 VISUALIZATION
% 1) Full scalp map of mean power (Reach+Grasp)
% 2) Highlighted clusters only
% ================================================================

% --- 1️⃣ Full topomap of mean desynchronization
fig1 = figure('Color','w');
topoplot(data_mean_RG, chanlocs, ...
         'electrodes','off', ...
         'style','fill', ...
         'shading','interp', ...
         'whitebk','on');
cmin = min(data_mean_RG);
cmax = max(data_mean_RG);

caxis([cmin cmax]);
colorbar;
title('');   

sgtitle('Whole Beta – Mean Normalized Power (Reach & Grasp)', ...
        'FontSize', 13, 'FontWeight', 'bold');
drawnow;

print(gcf, fullfile(save_dir_conv_whole_all, ...
    'Topoplot_WholeBeta_MeanPower_ReachGrasp'), '-dpng', '-r300');

savefig(gcf, fullfile(save_dir_conv_whole_all, ...
    'Topoplot_WholeBeta_MeanPower_ReachGrasp.fig'));



