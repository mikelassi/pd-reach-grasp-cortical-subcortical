%% diag_hold.m -- CMC across movement phases: reach(A), grasp-hold(C), pull(D)
% Expectation from literature: corticomuscular/STN-muscular coherence is
% strongest during sustained grip (around contact C), weakest during phasic reach.
clear; clc;
base='F:\Projects\Parkinson_ReachGrasp\Reprocessing';
subs={'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
DOM={'wue02','wue03','wue10','wue11'};
muscles={'iod','triceps','deltoid'}; ML={'IOD','Tri','Del'};
fs=400; [bbp,abp]=butter(4,[20 180]/(fs/2),'bandpass'); [bn,an]=butter(2,[49 51]/(fs/2),'stop');
M1R={'C4','C2','C6','FCC4h','FCC2h','FCC6h','CCP4h','CCP2h','CCP6h'};
M1L={'C3','C1','C5','FCC3h','FCC1h','FCC5h','CCP3h','CCP1h','CCP5h'};
rp=fullfile(fileparts(which('eeglab')),'plugins','Biosig3.8.4','biosig','maybe-missing'); if exist(rp,'dir'),rmpath(rp);end
eeglab nogui;
EP=load(fullfile(base,'RESULTS_final','Epochs','Epochs_allSubjects.mat'),'EPOCHS').EPOCHS;
NW=3;K=5;
phases={'A','reach',[0 0.5]; 'C','hold',[-0.1 0.5]; 'D','pull',[0 0.5]};
% storage: coh{phase}.SE{m}, .ME{m}  [nf x nsubj]
Z=struct(); for ph=1:3, Z(ph).SE=cell(3,1); Z(ph).ME=cell(3,1); end
fh=[];
for s=1:numel(subs)
  subj=subs{s}; isdom=ismember(subj,DOM); if isdom,M1set=M1R;else,M1set=M1L;end
  ep=fullfile(base,subj,'Preprocessed','EEG'); mp=fullfile(base,subj,'01_Extracted','EMG_KIN');
  lp=fullfile(base,subj,'Preprocessed','LFP'); evp=fullfile(base,subj,'02_Kinematics','Events');
  ef=dir(fullfile(ep,'*_manual.set')); mf=dir(fullfile(mp,'*.set')); lf=dir(fullfile(lp,'*_wEv.set'));
  if strcmp(subj,'wue06'), keep=setdiff(1:numel(ef),3); ef=ef(keep); mf=mf(keep); lf=lf(keep); end
  stn_ch=EP.(subj).ch_contra;
  for ph=1:3
    plet=phases{ph,1}; win=phases{ph,3}; pat=['^' plet '_T\d+$'];
    accM1=[]; accSTN=[]; accE=cell(3,1);
    for f=1:min([numel(ef) numel(mf) numel(lf)])
      EEG=pop_loadset('filename',ef(f).name,'filepath',ep);
      EMG=pop_loadset('filename',mf(f).name,'filepath',mp);
      LFP=pop_loadset('filename',lf(f).name,'filepath',lp); LFPc=pop_select(LFP,'channel',stn_ch);
      [~,b2]=fileparts(mf(f).name); km=load(fullfile(evp,[b2 '_kinematic_block.mat']),'events');
      ev=struct('type',{},'latency',{},'urevent',{}); k=0; L={'A','B','C','D','E','F'};
      for tr=1:numel(km.events), for e=1:6, if isfield(km.events,L{e}), la=km.events(tr).(L{e}); if ~isempty(la)&&~isnan(la), k=k+1; ev(k).type=sprintf('%s_T%d',L{e},tr); ev(k).latency=double(la); ev(k).urevent=k; end, end, end, end
      lab=lower({EMG.chanlocs.labels});
      if isdom,sel=contains(lab,'_dominant_emg')&~contains(lab,'non_dominant');else,sel=contains(lab,'_emg');end
      chs=find(sel); mch=zeros(1,3); for m=1:3,h=find(contains(lab(chs),muscles{m}),1);if ~isempty(h),mch(m)=chs(h);end;end
      emg=zeros(3,EMG.pnts); for m=1:3,if mch(m)==0,continue;end;g=double(EMG.data(mch(m),:));g=g-mean(g);g=filtfilt(bbp,abp,g);g=filtfilt(bn,an,g);emg(m,:)=abs(g);end
      EMGp=EMG;EMGp.data=emg;EMGp.nbchan=3;EMGp.chanlocs=EMGp.chanlocs(1:3);EMGp.event=ev;EMGp=eeg_checkset(EMGp);
      tE=ut(EEG.event,pat); tL=ut(LFPc.event,pat); tM=ut(EMGp.event,pat);
      if isempty(tE)||isempty(tL)||isempty(tM), continue; end
      try, Eep=pop_epoch(EEG,tE,win); Lep=pop_epoch(LFPc,tL,win); Mep=pop_epoch(EMGp,tM,win); catch, continue; end
      trE=tn(Eep); trL=tn(Lep); trM=tn(Mep);
      comm=intersect(intersect(trE,trL),trM); if isempty(comm),continue;end
      [~,iE]=ismember(comm,trE);[~,iL]=ismember(comm,trL);[~,iM]=ismember(comm,trM);
      ns=min([size(Eep.data,2) size(Lep.data,2) size(Mep.data,2)]);
      labE=lower({EEG.chanlocs.labels}); m1i=[]; for c=M1set,h=find(strcmp(labE,lower(c{1})),1);if ~isempty(h),m1i(end+1)=h;end;end
      zr=@(B)(B-mean(B,2))./max(std(B,0,2),eps);
      m1=zeros(1,ns,numel(comm)); for q=1:numel(comm), m1(1,:,q)=mean(zr(squeeze(Eep.data(m1i,1:ns,iE(q)))),1); end
      accM1=cat(3,accM1,m1); accSTN=cat(3,accSTN,Lep.data(1,1:ns,iL));
      for m=1:3,accE{m}=cat(3,accE{m},Mep.data(m,1:ns,iM));end
    end
    if isempty(accM1)||size(accM1,3)<3, continue; end
    for m=1:3
      Xs=[accSTN;accE{m}]; [S,ff]=cs_csd_multitaper(Xs,fs,NW,K); [C,fhh]=cs_coherence(S,ff); Z(ph).SE{m}=[Z(ph).SE{m},squeeze(C(1,2,:))]; if isempty(fh),fh=fhh;end
      Xm=[accM1;accE{m}];  [S,ff]=cs_csd_multitaper(Xm,fs,NW,K); C=cs_coherence(S,ff); Z(ph).ME{m}=[Z(ph).ME{m},squeeze(C(1,2,:))];
    end
  end
end
ci=1-0.05^(1/(K*27-1)); fprintf('approx CI=%.4f\n\n',ci);
bl=fh>=2&fh<=7; bb=fh>=13&fh<=30; bh=fh>=30&fh<=45;
for ph=1:3
  fprintf('=== phase %s (%s) grand-avg PEAK low/beta/hibeta ===\n',phases{ph,2},mat2str(phases{ph,3}));
  for m=1:3, v=mean(Z(ph).ME{m},2); fprintf('  M1 -%-4s: %.3f / %.3f / %.3f\n',ML{m},max(v(bl)),max(v(bb)),max(v(bh))); end
  for m=1:3, v=mean(Z(ph).SE{m},2); fprintf('  STN-%-4s: %.3f / %.3f / %.3f\n',ML{m},max(v(bl)),max(v(bb)),max(v(bh))); end
end
function t=ut(ev,pat), ty={ev.type}; t=unique(ty(~cellfun(@isempty,regexp(ty,pat)))); end
function trn=tn(ep), n=ep.trials; trn=nan(1,n); for i=1:n, et=ep.epoch(i).eventtype; el=ep.epoch(i).eventlatency; if ~iscell(et),et={et};end; if iscell(el),el=cell2mat(el);end; [~,j]=min(abs(el)); tok=regexp(et{j},'_T(\d+)$','tokens','once'); if ~isempty(tok),trn(i)=str2double(tok{1});end; end, end
