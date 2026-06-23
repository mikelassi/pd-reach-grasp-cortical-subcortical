%% cortex_stn_muscle_connectivity.m
% =========================================================================
% MULTIVARIATE CORTEX - STN - MUSCLE INTERACTION
% Parkinson's disease reach-to-grasp (8 STN-DBS patients, OFF medication).
%
% Computes, on the SAME trials and from a SINGLE multitaper cross-spectral
% density (CSD) matrix per condition, three complementary connectivity
% measures over a network of nodes:
%
%   NODES
%     * Cortical ROIs (groups of scalp EEG channels), laterality-corrected
%       so "contralateral" always refers to the hemisphere opposite the
%       moving hand:  M1_contra, PM_contra (premotor), S1par_contra
%       (post-central / parietal), SMA (mesial midline), M1_ipsi (control).
%     * STN  — contralateral bipolar LFP channel (index from Epoching.m).
%     * Muscles — contralateral IOD, Triceps, Deltoid (EMG 20-180 Hz).
%
%   MEASURES (all from the same CSD, hence mutually consistent)
%     1. Coherence            — undirected linear coupling.
%     2. Partial coherence    — coupling with all OTHER network nodes
%                               regressed out; separates DIRECT from
%                               INDIRECT links (e.g. is cortico-muscular
%                               coherence routed through the STN?).
%     3. Granger causality    — directed, frequency-resolved drive, via
%                               non-parametric Wilson factorisation.
%
%   CONDITIONS (windows reused from the lab's CMC / TF pipeline)
%     Movement : [0, +1.0]  s around movement onset (event A)
%     Pull     : [0, +0.875]s around pull onset     (event D)
%     Baseline : [-0.5, 0]  s pre-movement rest      (event A)
%
%   FOCUSED TRIAD  [M1_contra, STN, Muscle]  is analysed in detail
%   (3-node partial coherence conditions exactly on the third node), with
%   trial-shuffle surrogate significance for the Granger spectra.
%
% Estimators validated on synthetic ground truth in validate_connectivity.m
% (run that first). Outputs -> RESULTS_DIR/CortexSTNMuscle/.
%
% Author: Michael Lassi  (additional_analyses)
% =========================================================================

clear; clc;

%% ===== PATHS (auto-detect drive: data moved H: -> F:) =====
BASE_CANDIDATES = { 'F:\Projects\Parkinson_ReachGrasp\Reprocessing', ...
                    'H:\Parkinson_ReachGrasp\Reprocessing' };
BASE_PATH = '';
for i = 1:numel(BASE_CANDIDATES)
    if exist(BASE_CANDIDATES{i}, 'dir'), BASE_PATH = BASE_CANDIDATES{i}; break; end
end
assert(~isempty(BASE_PATH), 'Could not locate the Reprocessing data folder.');
fprintf('Data base: %s\n', BASE_PATH);

EEG_SUB     = fullfile('Preprocessed', 'EEG');     % *_manual.set (events embedded)
EMG_SUB     = fullfile('01_Extracted', 'EMG_KIN'); % *.set (no kin events)
LFP_SUB     = fullfile('Preprocessed', 'LFP');     % *_wEv.set (STN + events)
EVENTS_SUB  = fullfile('02_Kinematics', 'Events'); % *_kinematic_block.mat (native EMG-timeline events)
EPOCHS_MAT  = fullfile(BASE_PATH, 'RESULTS_final', 'Epochs', 'Epochs_allSubjects.mat');
RESULTS_DIR = fullfile(BASE_PATH, 'RESULTS_final', 'CortexSTNMuscle');

%% ===== RUN OPTIONS =====
QUICK_TEST = false;     % true: single subject, no surrogates (pipeline smoke test)
N_SURR     = 100;       % trial-shuffle surrogates for triad Granger (0 = skip;
                        % computed for the Movement condition, which is plotted)
SAVE_FIGS  = true;

%% ===== SUBJECTS & LATERALITY =====
SUBJECTS = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
% Hemisphere CONTRALATERAL to the moving hand. GROUND TRUTH is the wrist used
% for event segmentation in A02_events_segmentation.m (lines 91-96): wue05 and
% wue09 use the LEFT wrist; all others the RIGHT wrist.
%   left-hand movers  -> contralateral = RIGHT hemisphere (even electrodes)
%   right-hand movers -> contralateral = LEFT  hemisphere (odd  electrodes)
% NOTE (fixed): a previous version listed {wue02,wue03,wue10,wue11} here by
% conflating the EMG '_dominant_emg' subject set with the moving hand. That is
% WRONG and inverted contra/ipsi for 5/8 subjects.
LEFT_HAND_MOVERS  = {'wue05','wue09'};                    % contra = RIGHT
% (all others are right-hand movers -> contra = LEFT)
DOMINANT_SUBJECTS = {'wue02','wue03','wue10','wue11'};   % EMG '_dominant_emg' set (NOT the moving hand)

% wue06 block 3 EEG corrupted (per CMC_EEG_analysis.m)
SKIP_BLOCKS = containers.Map({'wue06'}, {3});

if QUICK_TEST, SUBJECTS = {'wue02'}; N_SURR = 0; end
n_subj = numel(SUBJECTS);

%% ===== CORTICAL ROIs (channel groups; hemisphere chosen per subject) =====
% Left/right channel lists; the script picks contra/ipsi per subject.
ROI_L = struct( ...
    'M1',   {{'C3','C1','C5','FCC3h','FCC1h','FCC5h','CCP3h','CCP1h','CCP5h'}}, ...
    'PM',   {{'FC3','FC1','FC5','FFC3h','FFC5h'}}, ...
    'S1par',{{'CP3','CP1','CP5','P3','CPP3h','CPP5h'}});
ROI_R = struct( ...
    'M1',   {{'C4','C2','C6','FCC4h','FCC2h','FCC6h','CCP4h','CCP2h','CCP6h'}}, ...
    'PM',   {{'FC4','FC2','FC6','FFC4h','FFC6h'}}, ...
    'S1par',{{'CP4','CP2','CP6','P4','CPP4h','CPP6h'}});
ROI_MID = {'FCz','Cz','Fz','FCC1h','FCC2h','FFC1h','FFC2h'};   % SMA / mesial (midline)

ROI_NAMES = {'M1_contra','PM_contra','S1par_contra','SMA','M1_ipsi'};
n_roi = numel(ROI_NAMES);

%% ===== MUSCLES =====
TARGET_MUSCLES = {'iod','triceps','deltoid'};
MUSCLE_LABELS  = {'IOD','Triceps','Deltoid'};
n_mus = numel(TARGET_MUSCLES);

% Full node ordering: [ROIs, STN, muscles]
NODE_LABELS = [ROI_NAMES, {'STN'}, MUSCLE_LABELS];
n_node = numel(NODE_LABELS);
IDX_M1c = find(strcmp(NODE_LABELS,'M1_contra'));
IDX_STN = find(strcmp(NODE_LABELS,'STN'));
IDX_MUS = (n_roi+2) : n_node;   % muscle node indices

%% ===== EMG PREPROCESSING =====
BP_LOW = 20; BP_HIGH = 180; BP_ORDER = 4;
NOTCH_FREQS = [50 100]; NOTCH_BW = 2; NOTCH_ORDER = 2;
EMG_RECTIFY = true;   % full-wave rectify band-passed EMG before connectivity.
                      % Demodulates motor-unit firing so the common ~15-30 Hz
                      % descending drive appears in the EMG spectrum -> the
                      % standard input for corticomuscular coherence
                      % (Mima & Hallett 1999; Boonstra & Breakspear 2012).

%% ===== ANALYSIS WINDOWS (s) =====
WIN_MOVE = [0, 1.000];  WIN_PULL = [0, 0.875];  WIN_BASE = [-0.5, 0.000];
COND_NAMES = {'Movement','Pull','Baseline'};
n_cond = numel(COND_NAMES);

%% ===== SPECTRAL SETTINGS =====
NW = 3;  K = 2*NW - 1;          % multitaper time-bandwidth, tapers
FREQ_MIN = 2;  FREQ_MAX = 100;  % display range
ALPHA = 0.05;

% Bands chosen to capture (i) low-frequency corticomuscular common drive
% (2-7 Hz), dominant during dynamic movement; (ii) alpha; (iii) beta
% (cortico-subthalamic); (iv) high-beta/low-gamma (STN-EMG, Marsden 2001).
BANDS = struct( ...
    'name', {'low (2-7)','alpha (8-12)','beta (13-30)','high-beta (30-45)'}, ...
    'lim',  {[2 7],       [8 12],        [13 30],       [30 45]});
n_band = numel(BANDS);
BETA_IDX = find(strcmp({BANDS.name},'beta (13-30)'));
LOW_IDX  = find(strcmp({BANDS.name},'low (2-7)'));

%% ===== INIT =====
remove_biosig_stubs();
eeglab nogui;
if ~exist(RESULTS_DIR,'dir'), mkdir(RESULTS_DIR); end

% Load STN contralateral channel indices from the canonical epoching output
assert(exist(EPOCHS_MAT,'file')>0, 'Missing %s', EPOCHS_MAT);
tmp = load(EPOCHS_MAT,'EPOCHS'); EPOCHS = tmp.EPOCHS; clear tmp;

% Per-subject, per-condition storage
% Band-averaged full-network matrices: [n_node x n_node x n_band]
coh_band  = nan(n_node, n_node, n_band, n_cond, n_subj);
pcoh_band = nan(n_node, n_node, n_band, n_cond, n_subj);
gc_band   = nan(n_node, n_node, n_band, n_cond, n_subj);
% Triad spectra [n_freq x 3pairs x n_mus] per condition (freq axis per cond)
triad = struct();
for c = 1:n_cond
    triad(c).f = [];
    triad(c).coh  = cell(n_subj, n_mus);   % each [nf x 3]  pairs: CM1-STN, STN-MUS, M1-MUS
    triad(c).pcoh = cell(n_subj, n_mus);
    triad(c).gc   = cell(n_subj, n_mus);   % each [nf x 6]  directed pairs (see PAIR_DIR)
    triad(c).gc_thr = cell(n_subj, n_mus); % surrogate 95% threshold [nf x 6] or []
end
dof_subj = nan(n_subj, n_cond);            % K * n_trials (coherence DOF)

% Triad pair definitions (nodes ordered [M1c, STN, MUS] = 1,2,3)
PAIR_UND = [1 2; 2 3; 1 3];                % undirected pairs
PAIR_UND_LBL = {'M1-STN','STN-MUS','M1-MUS'};
PAIR_DIR = [2 1; 1 2; 3 2; 2 3; 3 1; 1 3]; % [to from]: STN->M1,M1->STN,MUS->STN,STN->MUS,MUS->M1,M1->MUS
PAIR_DIR_LBL = {'STN->M1','M1->STN','MUS->STN','STN->MUS','MUS->M1','M1->STN_dummy'};
PAIR_DIR_LBL{6} = 'M1->MUS';

%% ===== SUBJECT LOOP =====
for s = 1:n_subj
    subj = SUBJECTS{s};
    fprintf('\n=================  Subject %s  (%d/%d)  =================\n', subj, s, n_subj);

    eeg_path = fullfile(BASE_PATH, subj, EEG_SUB);
    emg_path = fullfile(BASE_PATH, subj, EMG_SUB);
    lfp_path = fullfile(BASE_PATH, subj, LFP_SUB);
    evt_path = fullfile(BASE_PATH, subj, EVENTS_SUB);
    eeg_files = dir(fullfile(eeg_path,'*_manual.set'));
    emg_files = dir(fullfile(emg_path,'*.set'));
    lfp_files = dir(fullfile(lfp_path,'*_wEv.set'));
    if isempty(eeg_files)||isempty(emg_files)||isempty(lfp_files)
        warning('Missing files for %s — skipping.', subj); continue;
    end

    % Block skipping (corrupted)
    if isKey(SKIP_BLOCKS, subj)
        kb = SKIP_BLOCKS(subj);
        keep = setdiff(1:numel(eeg_files), kb);
        eeg_files = eeg_files(keep);
        if numel(emg_files)>=max(keep), emg_files = emg_files(keep); end
        if numel(lfp_files)>=max(keep), lfp_files = lfp_files(keep); end
        fprintf('  Skipping block(s) %s (corrupted).\n', mat2str(kb));
    end

    use_dominant = ismember(subj, DOMINANT_SUBJECTS);
    contra_right = ismember(subj, LEFT_HAND_MOVERS);   % contra hemisphere = right?
    if contra_right
        ROIc = ROI_R; ROIi = ROI_L;     % contra = right, ipsi = left
    else
        ROIc = ROI_L; ROIi = ROI_R;     % contra = left,  ipsi = right
    end

    % STN contralateral channel index
    if isfield(EPOCHS, subj) && isfield(EPOCHS.(subj),'ch_contra')
        stn_ch = EPOCHS.(subj).ch_contra;
    else
        stn_ch = 1;
    end

    % Accumulators across blocks (cells of per-condition 3D arrays appended in trials)
    accE = {[],[]};   % EEG  onset / pull : [nch x nsamp x ntr]
    accL = {[],[]};   % STN  onset / pull : [1   x nsamp x ntr]
    accM = {[],[]};   % EMG  onset / pull : [nmus x nsamp x ntr]
    eeg_labels_ref = [];
    fs = NaN;

    nblk = min([numel(eeg_files), numel(emg_files), numel(lfp_files)]);
    for f = 1:nblk
        % ---- load ----
        EEG = pop_loadset('filename', eeg_files(f).name, 'filepath', eeg_path);
        EMG = pop_loadset('filename', emg_files(f).name, 'filepath', emg_path);
        LFP = pop_loadset('filename', lfp_files(f).name, 'filepath', lfp_path);
        if isnan(fs), fs = EEG.srate; end
        if isempty(eeg_labels_ref), eeg_labels_ref = {EEG.chanlocs.labels}; end

        % ---- STN: keep contralateral bipolar channel only ----
        if stn_ch > LFP.nbchan, stn_ch = 1; end
        LFP = pop_select(LFP, 'channel', stn_ch);

        % ---- EMG: select contralateral target muscles, BP+notch ----
        ch_labels = lower({EMG.chanlocs.labels});
        % NOTE: 'iod_non_dominant_emg' also contains the substring '_dominant_emg',
        % so the non-dominant (ipsilateral) channel must be explicitly excluded —
        % otherwise the first match is the WRONG arm for dominant-set subjects.
        if use_dominant
            sel = contains(ch_labels,'_dominant_emg') & ~contains(ch_labels,'non_dominant');
        else
            sel = contains(ch_labels,'_emg');
        end
        ch_all = find(sel);
        mch = zeros(1,n_mus);
        for m = 1:n_mus
            hit = find(contains(ch_labels(ch_all), TARGET_MUSCLES{m}),1);
            if ~isempty(hit), mch(m) = ch_all(hit); end
        end
        if ~any(mch), warning('  Block %d: no target EMG — skip.', f); continue; end
        [bbp,abp] = butter(BP_ORDER,[BP_LOW BP_HIGH]/(fs/2),'bandpass');
        ncf = cell(numel(NOTCH_FREQS),2);
        for ni=1:numel(NOTCH_FREQS)
            fn=NOTCH_FREQS(ni);
            if fn<fs/2
                bw=NOTCH_BW/(fs/2); fc=fn/(fs/2);
                [ncf{ni,1},ncf{ni,2}]=butter(NOTCH_ORDER,[fc-bw/2 fc+bw/2],'stop');
            end
        end
        emg_bp = zeros(n_mus, EMG.pnts);
        for m=1:n_mus
            if mch(m)==0, continue; end
            sig=double(EMG.data(mch(m),:)); sig=sig-mean(sig);
            sig=filtfilt(bbp,abp,sig);
            for ni=1:numel(NOTCH_FREQS)
                if ~isempty(ncf{ni,1}), sig=filtfilt(ncf{ni,1},ncf{ni,2},sig); end
            end
            if EMG_RECTIFY, sig = abs(sig); end    % full-wave rectification
            emg_bp(m,:)=sig;
        end
        EMGp = EMG; EMGp.data=emg_bp; EMGp.nbchan=n_mus; EMGp.pnts=size(emg_bp,2);
        EMGp.chanlocs = EMGp.chanlocs(1:n_mus);
        for m=1:n_mus, EMGp.chanlocs(m).labels = MUSCLE_LABELS{m}; end

        % ---- EMG events: use the NATIVE kinematic events (EMG-file timeline) ----
        % CRITICAL: the kinematic events were detected on the wrist velocity IN the
        % EMG_KIN file (A02_events_segmentation.m) and stored in samples of THAT file.
        % The Preprocessed LFP/EEG were start-cropped during preprocessing, so their
        % embedded event latencies live in a DIFFERENT (cropped) timeline. Injecting
        % LFP-timeline latencies into the uncropped EMG file misaligns EMG by several
        % seconds and destroys corticomuscular coherence. We therefore epoch EMG with
        % its own native events and EEG/LFP with their own embedded events — each in
        % its own correct timeline, all marking the same physical movement events.
        [~, emg_base] = fileparts(emg_files(f).name);
        kin_mat = fullfile(evt_path, [emg_base '_kinematic_block.mat']);
        if ~exist(kin_mat,'file')
            warning('  Block %d: native kinematic events not found (%s) — skip.', f, kin_mat); continue;
        end
        KM = load(kin_mat, 'events');
        emg_ev = build_kin_events(KM.events);   % A_T*/D_T* in EMG samples
        if isempty(emg_ev), warning('  Block %d: empty native events — skip.', f); continue; end
        EMGp.event = emg_ev; EMGp = eeg_checkset(EMGp);

        % ---- epoch all three on A (onset) and D (pull), each in its own timeline ----
        % Each modality is epoched with its own events; epochs are then matched
        % across modalities by TRIAL NUMBER (the _T<n> index), not by position, so
        % that any epoch dropped in one modality cannot misalign the others.
        for cc = 1:2   % 1=onset, 2=pull
            if cc==1, pat='^A_T\d+$'; win=[-0.5 1.0]; else, pat='^D_T\d+$'; win=[-0.3 0.875]; end
            try
                [Edata, trE] = epoch_with_trials(EEG,  pat, win);
                [Ldata, trL] = epoch_with_trials(LFP,  pat, win);
                [Mdata, trM] = epoch_with_trials(EMGp, pat, win);
            catch ME
                warning('  Block %d %s: epoch failed (%s)', f, COND_NAMES{cc}, ME.message); continue;
            end
            if isempty(trE)||isempty(trL)||isempty(trM), continue; end
            common = intersect(intersect(trE, trL), trM);   % trial numbers in all 3
            if isempty(common), continue; end
            [~, iE] = ismember(common, trE);
            [~, iL] = ismember(common, trL);
            [~, iM] = ismember(common, trM);
            ns = min([size(Edata,2), size(Ldata,2), size(Mdata,2)]);
            accE{cc} = cat(3, accE{cc}, double(Edata(:,1:ns,iE)));
            accL{cc} = cat(3, accL{cc}, double(Ldata(:,1:ns,iL)));
            accM{cc} = cat(3, accM{cc}, double(Mdata(:,1:ns,iM)));
        end
    end % block

    if isempty(accE{1}), warning('Subject %s: no onset epochs — skipping.', subj); continue; end

    % ---- build ROI signals (z-score each channel per trial, then average) ----
    lab = lower(eeg_labels_ref);
    roi_chan = cell(1,n_roi);
    roi_chan{1} = pick_idx(lab, ROIc.M1);
    roi_chan{2} = pick_idx(lab, ROIc.PM);
    roi_chan{3} = pick_idx(lab, ROIc.S1par);
    roi_chan{4} = pick_idx(lab, ROI_MID);
    roi_chan{5} = pick_idx(lab, ROIi.M1);
    for r=1:n_roi
        if isempty(roi_chan{r})
            warning('  ROI %s: no channels matched for %s.', ROI_NAMES{r}, subj);
        end
    end

    % time vectors for the two epoch types
    t_on = (0:size(accE{1},2)-1)/fs - 0.5;
    if ~isempty(accE{2}), t_pu = (0:size(accE{2},2)-1)/fs - 0.3; else, t_pu=[]; end

    % ---- per condition ----
    for c = 1:n_cond
        switch c
            case 1, cc=1; mask = t_on>=WIN_MOVE(1) & t_on<WIN_MOVE(2);
            case 2, cc=2; if isempty(t_pu), continue; end
                          mask = t_pu>=WIN_PULL(1) & t_pu<WIN_PULL(2);
            case 3, cc=1; mask = t_on>=WIN_BASE(1) & t_on<WIN_BASE(2);
        end
        E = accE{cc}; L = accL{cc}; M = accM{cc};
        if isempty(E), continue; end
        ntr = size(E,3);
        if ntr < 3, warning('  %s %s: only %d trials — skipping.', subj, COND_NAMES{c}, ntr); continue; end

        % ROI signals: [n_roi x nsamp x ntr]
        nsw = sum(mask);
        roi_sig = zeros(n_roi, nsw, ntr);
        for r=1:n_roi
            ci = roi_chan{r};
            if isempty(ci), roi_sig(r,:,:) = 0; continue; end
            blk = E(ci, mask, :);                          % [nc x nsw x ntr]
            % z-score each channel per trial over the window, then mean
            mu = mean(blk,2); sd = std(blk,0,2); sd(sd==0)=1;
            blkz = (blk - mu)./sd;
            roi_sig(r,:,:) = mean(blkz,1);
        end
        stn_sig = L(1, mask, :);                            % [1 x nsw x ntr]
        emg_sig = M(:, mask, :);                            % [nmus x nsw x ntr]

        X = cat(1, roi_sig, stn_sig, emg_sig);              % [n_node x nsw x ntr]

        % ---- CSD + connectivity (full network) ----
        [S, ff] = cs_csd_multitaper(X, fs, NW, K);
        [C, fh]  = cs_coherence(S, ff);
        [PC, ~]  = cs_partial_coherence(S, ff);
        [G, ~]   = cs_granger(S, ff);

        % band-average and store
        coh_band(:,:,:,c,s)  = cs_bandavg(C,  fh, BANDS);
        pcoh_band(:,:,:,c,s) = cs_bandavg(PC, fh, BANDS);
        gc_band(:,:,:,c,s)   = cs_bandavg(G,  fh, BANDS);
        dof_subj(s,c) = K * ntr;

        % ---- focused triads [M1c, STN, muscle] ----
        if isempty(triad(c).f), triad(c).f = fh; end
        for m = 1:n_mus
            idx = [IDX_M1c, IDX_STN, IDX_MUS(m)];
            St  = S(idx, idx, :);
            Ct  = cs_coherence(St, ff);
            Pt  = cs_partial_coherence(St, ff);
            Gt  = cs_granger(St, ff);
            % undirected pairs
            cohp = zeros(numel(fh),3); pcp = zeros(numel(fh),3);
            for p=1:3
                cohp(:,p) = squeeze(Ct(PAIR_UND(p,1),PAIR_UND(p,2),:));
                pcp(:,p)  = squeeze(Pt(PAIR_UND(p,1),PAIR_UND(p,2),:));
            end
            gcp = zeros(numel(fh),6);
            for p=1:6
                gcp(:,p) = squeeze(Gt(PAIR_DIR(p,1),PAIR_DIR(p,2),:));
            end
            triad(c).coh{s,m}  = cohp;
            triad(c).pcoh{s,m} = pcp;
            triad(c).gc{s,m}   = gcp;

            % ---- surrogate Granger threshold (trial-shuffle) ----
            % Computed for Movement only (the condition shown with thresholds);
            % observed Movement GC above uses identical estimator settings.
            if N_SURR > 0 && c == 1
                Xt = X(idx,:,:);
                surr = zeros(numel(fh),6,N_SURR);
                for q=1:N_SURR
                    Xs = Xt;
                    for nd=1:3
                        Xs(nd,:,:) = Xt(nd,:,randperm(ntr));   % break cross-node trial pairing
                    end
                    Ss = cs_csd_multitaper(Xs, fs, NW, K);
                    Gs = cs_granger(Ss, ff);
                    for p=1:6, surr(:,p,q)=squeeze(Gs(PAIR_DIR(p,1),PAIR_DIR(p,2),:)); end
                end
                triad(c).gc_thr{s,m} = prctile(surr, 95, 3);   % [nf x 6]
            end
        end
        fprintf('  %-9s: %d trials, %d nodes, nf=%d done\n', COND_NAMES{c}, ntr, n_node, numel(fh));
    end
end % subject

%% ===== GRAND AVERAGE =====
GA.coh  = nanmean(coh_band,5);   GA.coh_n  = sum(~isnan(coh_band(1,2,1,:,:)),5);
GA.pcoh = nanmean(pcoh_band,5);
GA.gc   = nanmean(gc_band,5);
ci_line = nan(1,n_cond);
for c=1:n_cond
    d = dof_subj(:,c); d = d(~isnan(d));
    if ~isempty(d), ci_line(c) = 1 - ALPHA^(1/(median(d)-1)); end
end

%% ===== SAVE =====
RES.subjects=SUBJECTS; RES.node_labels=NODE_LABELS; RES.cond=COND_NAMES;
RES.bands=BANDS; RES.NW=NW; RES.K=K; RES.alpha=ALPHA;
RES.coh_band=coh_band; RES.pcoh_band=pcoh_band; RES.gc_band=gc_band;
RES.GA=GA; RES.ci_line=ci_line; RES.dof_subj=dof_subj;
RES.triad=triad; RES.pair_und=PAIR_UND_LBL; RES.pair_dir=PAIR_DIR_LBL;
RES.muscles=MUSCLE_LABELS; RES.roi_names=ROI_NAMES;
save(fullfile(RESULTS_DIR,'CortexSTNMuscle_results.mat'),'RES','-v7.3');
fprintf('\nSaved results -> %s\n', fullfile(RESULTS_DIR,'CortexSTNMuscle_results.mat'));

%% ===== FIGURES =====
if SAVE_FIGS
    make_figures(RESULTS_DIR);
end

fprintf('\nDONE.\n');

%% ========================= LOCAL FUNCTIONS =============================
function remove_biosig_stubs()
    if ~isempty(which('eeglab'))
        bs = fullfile(fileparts(which('eeglab')),'plugins','Biosig3.8.4','biosig','maybe-missing');
        if exist(bs,'dir'), rmpath(bs); end
    end
end

function t = unique_types(ev, pat)
    types = {ev.type};
    t = unique(types(~cellfun(@isempty, regexp(types, pat))));
end

function [data, trnum] = epoch_with_trials(EEGset, pat, win)
% Epoch EEGset on all event types matching pat, returning the data
% [nch x nsamp x ntrials] and the trial number (the _T<n> index) that each
% epoch is time-locked to. Enables trial-number-based alignment across modalities.
    types = unique_types(EEGset.event, pat);
    if isempty(types), data = []; trnum = []; return; end
    ep = pop_epoch(EEGset, types, win);
    ntr = ep.trials; trnum = nan(1, ntr);
    for i = 1:ntr
        et = ep.epoch(i).eventtype;     el = ep.epoch(i).eventlatency;
        if ~iscell(et), et = {et}; end
        if iscell(el), el = cell2mat(el); end
        [~, j] = min(abs(el));          % time-locking event is at latency 0
        tok = regexp(et{j}, '_T(\d+)$', 'tokens', 'once');
        if ~isempty(tok), trnum(i) = str2double(tok{1}); end
    end
    data = ep.data;
end

function ev = build_kin_events(events)
% Build an EEGLAB event struct (A_T*/D_T* ... F_T*) from the per-trial
% kinematic 'events' struct saved by A02_events_segmentation.m. Latencies are
% in samples of the EMG_KIN file (the native timeline of these events).
    ev = struct('type',{},'latency',{},'urevent',{});
    letters = {'A','B','C','D','E','F'};
    n = numel(events); k = 0;
    for tr = 1:n
        for e = 1:numel(letters)
            if ~isfield(events, letters{e}), continue; end
            lat = events(tr).(letters{e});
            if ~isempty(lat) && ~isnan(lat)
                k = k + 1;
                ev(k).type    = sprintf('%s_T%d', letters{e}, tr);
                ev(k).latency = double(lat);
                ev(k).urevent = k;
            end
        end
    end
end

function idx = pick_idx(lab_lower, wanted)
    idx = [];
    for k=1:numel(wanted)
        h = find(strcmp(lab_lower, lower(wanted{k})),1);
        if ~isempty(h), idx(end+1)=h; end %#ok<AGROW>
    end
end
