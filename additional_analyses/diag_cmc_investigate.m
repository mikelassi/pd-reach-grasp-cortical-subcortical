%% diag_cmc_investigate.m
% Why is corticomuscular / subthalamo-muscular coherence weak?
% Decisive controls across all 8 subjects:
%   (1) EMG_KIN vs LFP file length (cross-file timeline match)
%   (2) movement-locked EMG envelope (events correctly aligned + activation)
%   (3) INTERMUSCULAR coherence (positive control: common drive, should be high)
%   (4) full-spectrum STN-EMG and M1-EMG coherence (look beyond 13-30 Hz)
%   (5) compare windows: onset-transient [0,0.3] vs reach [0,1.0] vs pull/hold
% Uses the validated cs_* estimators.
clear; clc;
base='F:\Projects\Parkinson_ReachGrasp\Reprocessing';
subs={'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
DOM={'wue02','wue03','wue10','wue11'};                 % dominant EMG set; contra=RIGHT
muscles={'iod','triceps','deltoid'}; ML={'IOD','Tri','Del'};
fs=400; [bbp,abp]=butter(4,[20 180]/(fs/2),'bandpass'); [bn,an]=butter(2,[49 51]/(fs/2),'stop');
M1R={'C4','C2','C6','FCC4h','FCC2h','FCC6h','CCP4h','CCP2h','CCP6h'};
M1L={'C3','C1','C5','FCC3h','FCC1h','FCC5h','CCP3h','CCP1h','CCP5h'};
rp=fullfile(fileparts(which('eeglab')),'plugins','Biosig3.8.4','biosig','maybe-missing'); if exist(rp,'dir'),rmpath(rp);end
eeglab nogui;

NW=3; K=5;
env_acc=[]; tenv=[];
IM=cell(3,1); SE=cell(3,1); ME=cell(3,1);      % intermuscular, STN-EMG, M1-EMG spectra (per muscle/pair)
fkeep=[];
fprintf('%-7s %8s %8s\n','subj','EMGpnts','LFPpnts');
for s=1:numel(subs)
  subj=subs{s};
  ep=fullfile(base,subj,'Preprocessed','EEG'); mp=fullfile(base,subj,'01_Extracted','EMG_KIN'); lp=fullfile(base,subj,'Preprocessed','LFP');
  ef=dir(fullfile(ep,'*_manual.set')); mf=dir(fullfile(mp,'*.set')); lf=dir(fullfile(lp,'*_wEv.set'));
  if isempty(ef)||isempty(mf)||isempty(lf), continue; end
  if strcmp(subj,'wue06'), keep=setdiff(1:numel(ef),3); ef=ef(keep); if numel(mf)>=3,mf=mf(keep);end; if numel(lf)>=3,lf=lf(keep);end; end
  isdom=ismember(subj,DOM); if isdom, M1set=M1R; else, M1set=M1L; end
  Etmp=load(fullfile(base,'RESULTS_final','Epochs','Epochs_allSubjects.mat'),'EPOCHS');
  stn_ch=Etmp.EPOCHS.(subj).ch_contra;
  accM1=[]; accSTN=[]; accE=cell(3,1);
  for f=1:min([numel(ef) numel(mf) numel(lf)])
    EEG=pop_loadset('filename',ef(f).name,'filepath',ep);
    EMG=pop_loadset('filename',mf(f).name,'filepath',mp);
    LFP=pop_loadset('filename',lf(f).name,'filepath',lp);
    if s<=numel(subs) && f==1, fprintf('%-7s %8d %8d\n',subj,EMG.pnts,LFP.pnts); end
    LFPc=pop_select(LFP,'channel',stn_ch);
    lab=lower({EMG.chanlocs.labels});
    if isdom, sel=contains(lab,'_dominant_emg')&~contains(lab,'non_dominant'); else, sel=contains(lab,'_emg'); end
    chs=find(sel); mch=zeros(1,3);
    for m=1:3, h=find(contains(lab(chs),muscles{m}),1); if ~isempty(h), mch(m)=chs(h); end, end
    emg=zeros(3,EMG.pnts);
    for m=1:3, if mch(m)==0,continue;end; g=double(EMG.data(mch(m),:)); g=g-mean(g); g=filtfilt(bbp,abp,g); g=filtfilt(bn,an,g); emg(m,:)=abs(g); end
    EMGp=EMG; EMGp.data=emg; EMGp.nbchan=3; EMGp.chanlocs=EMGp.chanlocs(1:3); EMGp=eeg_checkset(EMGp);
    kin=LFP.event(~cellfun(@isempty,regexp({LFP.event.type},'^A_T\d+$'))); EMGp.event=kin; EMGp=eeg_checkset(EMGp);
    tA=unique({kin.type}); tAe=unique({EEG.event(~cellfun(@isempty,regexp({EEG.event.type},'^A_T\d+$'))).type});
    Eep=pop_epoch(EEG,tAe,[-0.5 1.0]); Lep=pop_epoch(LFPc,tA,[-0.5 1.0]); Mep=pop_epoch(EMGp,tA,[-0.5 1.0]);
    nt=min([Eep.trials Lep.trials Mep.trials]); ns=min([size(Eep.data,2) size(Lep.data,2) size(Mep.data,2)]);
    labE=lower({EEG.chanlocs.labels}); m1i=[]; for c=M1set, h=find(strcmp(labE,lower(c{1})),1); if ~isempty(h),m1i(end+1)=h;end; end
    zr=@(B)(B-mean(B,2))./max(std(B,0,2),eps);
    m1=zeros(1,ns,nt); for t=1:nt, m1(1,:,t)=mean(zr(squeeze(Eep.data(m1i,1:ns,t))),1); end
    accM1=cat(3,accM1,m1); accSTN=cat(3,accSTN,Lep.data(1,1:ns,1:nt));
    for m=1:3, accE{m}=cat(3,accE{m}, Mep.data(m,1:ns,1:nt)); end
  end
  ns=size(accM1,2); t=(0:ns-1)/fs-0.5; if isempty(tenv),tenv=t;end
  % EMG envelope time course (mean rectified IOD across trials), normalised
  e=squeeze(mean(accE{1},3)); env_acc=[env_acc; e/max(e)];
  % movement window [0,1.0]
  mask=t>=0 & t<1.0;
  pairs=[1 2;1 3;2 3];
  for p=1:3
    X=[accE{pairs(p,1)}(1,mask,:); accE{pairs(p,2)}(1,mask,:)]; [S,ff]=cs_csd_multitaper(X,fs,NW,K); [C,fh]=cs_coherence(S,ff);
    IM{p}=[IM{p}, squeeze(C(1,2,:))]; if isempty(fkeep),fkeep=fh;end
  end
  for m=1:3
    Xs=[accSTN(1,mask,:); accE{m}(1,mask,:)]; [S,ff]=cs_csd_multitaper(Xs,fs,NW,K); C=cs_coherence(S,ff); SE{m}=[SE{m}, squeeze(C(1,2,:))];
    Xm=[accM1(1,mask,:);  accE{m}(1,mask,:)]; [S,ff]=cs_csd_multitaper(Xm,fs,NW,K); C=cs_coherence(S,ff); ME{m}=[ME{m}, squeeze(C(1,2,:))];
  end
end
ci=1-0.05^(1/(K*25-1));   % approx DOF
fprintf('\napprox 95%% CI (K*~25 trials): %.4f\n',ci);

b1330=fkeep>=13&fkeep<=30; b3045=fkeep>=30&fkeep<=45; blow=fkeep>=2&fkeep<=12;
fprintf('\n=== INTERMUSCULAR coherence (positive control), grand-avg peak ===\n');
pl={'IOD-Tri','IOD-Del','Tri-Del'};
for p=1:3, m=mean(IM{p},2); fprintf('  %-8s: low(2-12)=%.3f  beta(13-30)=%.3f  hibeta(30-45)=%.3f\n',pl{p},max(m(blow)),max(m(b1330)),max(m(b3045))); end
fprintf('\n=== STN-EMG coherence grand-avg peak ===\n');
for m=1:3, v=mean(SE{m},2); fprintf('  %-8s: low=%.3f beta=%.3f hibeta=%.3f\n',ML{m},max(v(blow)),max(v(b1330)),max(v(b3045))); end
fprintf('\n=== M1-EMG coherence grand-avg peak ===\n');
for m=1:3, v=mean(ME{m},2); fprintf('  %-8s: low=%.3f beta=%.3f hibeta=%.3f\n',ML{m},max(v(blow)),max(v(b1330)),max(v(b3045))); end

% EMG envelope time course
fprintf('\n=== EMG(IOD) envelope locked to onset (mean over subj) ===\n');
me=mean(env_acc,1); tt=tenv;
pre=me(tt>=-0.5&tt<0); post=me(tt>=0&tt<0.5);
fprintf('  pre-onset mean=%.3f  post-onset mean=%.3f  ratio=%.2f\n',mean(pre),mean(post),mean(post)/mean(pre));
[~,ip]=max(me); fprintf('  envelope peak at t=%.2f s\n',tt(ip));
save(fullfile(base,'RESULTS_final','CortexSTNMuscle','diag_cmc.mat'),'IM','SE','ME','fkeep','env_acc','tenv','pl','ML');
fprintf('\nsaved diag_cmc.mat\n');
