%% diag_native_group.m
% Group-level test of the NATIVE-event EMG alignment + full-spectrum CMC.
%   - crop offset per subject = (Preprocessed-LFP A latency) - (native EMG A latency)
%   - EMG envelope movement-locking (post/pre ratio) per subject, native events
%   - grand-avg STN-EMG / M1-EMG / intermuscular coherence, full spectrum, native
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
SE=cell(3,1); ME=cell(3,1); IM=cell(3,1); fh=[]; NW=3;K=5;
fprintf('%-7s %10s %10s %8s\n','subj','cropOffset','envRatio','envPk(s)');
for s=1:numel(subs)
  subj=subs{s}; isdom=ismember(subj,DOM); if isdom,M1set=M1R;else,M1set=M1L;end
  ep=fullfile(base,subj,'Preprocessed','EEG'); mp=fullfile(base,subj,'01_Extracted','EMG_KIN');
  lp=fullfile(base,subj,'Preprocessed','LFP'); evp=fullfile(base,subj,'02_Kinematics','Events');
  ef=dir(fullfile(ep,'*_manual.set')); mf=dir(fullfile(mp,'*.set')); lf=dir(fullfile(lp,'*_wEv.set'));
  if strcmp(subj,'wue06'), keep=setdiff(1:numel(ef),3); ef=ef(keep); mf=mf(keep); lf=lf(keep); end
  stn_ch=EP.(subj).ch_contra;
  accM1=[]; accSTN=[]; accE=cell(3,1); offs=[]; envr=[]; envpk=[];
  for f=1:min([numel(ef) numel(mf) numel(lf)])
    EEG=pop_loadset('filename',ef(f).name,'filepath',ep);
    EMG=pop_loadset('filename',mf(f).name,'filepath',mp);
    LFP=pop_loadset('filename',lf(f).name,'filepath',lp); LFPc=pop_select(LFP,'channel',stn_ch);
    [~,bn2]=fileparts(mf(f).name); km=load(fullfile(evp,[bn2 '_kinematic_block.mat']));
    % native EMG events
    evA=[];evD=[];
    for tr=1:numel(km.events)
      a=km.events(tr).A; d=km.events(tr).D;
      if ~isempty(a)&&~isnan(a), clear evt; evt.type=sprintf('A_T%d',tr);evt.latency=a;evt.urevent=tr;evA=[evA,evt];end
      if ~isempty(d)&&~isnan(d), clear evt2; evt2.type=sprintf('D_T%d',tr);evt2.latency=d;evt2.urevent=100+tr;evD=[evD,evt2];end
    end
    % crop offset: LFP A_T1 latency - native A_T1 latency
    la=LFP.event(strcmp({LFP.event.type},'A_T1')); na=km.events(1).A;
    if ~isempty(la)&&~isempty(na), offs(end+1)=la(1).latency-na; end
    lab=lower({EMG.chanlocs.labels});
    if isdom,sel=contains(lab,'_dominant_emg')&~contains(lab,'non_dominant');else,sel=contains(lab,'_emg');end
    chs=find(sel); mch=zeros(1,3); for m=1:3,h=find(contains(lab(chs),muscles{m}),1);if ~isempty(h),mch(m)=chs(h);end;end
    emg=zeros(3,EMG.pnts);
    for m=1:3,if mch(m)==0,continue;end;g=double(EMG.data(mch(m),:));g=g-mean(g);g=filtfilt(bbp,abp,g);g=filtfilt(bn,an,g);emg(m,:)=abs(g);end
    EMGp=EMG;EMGp.data=emg;EMGp.nbchan=3;EMGp.chanlocs=EMGp.chanlocs(1:3);EMGp.event=evA;EMGp=eeg_checkset(EMGp);
    tA=unique({evA.type}); tAe=unique({EEG.event(~cellfun(@isempty,regexp({EEG.event.type},'^A_T\d+$'))).type});
    tAl=unique({LFPc.event(~cellfun(@isempty,regexp({LFPc.event.type},'^A_T\d+$'))).type});
    Mep=pop_epoch(EMGp,tA,[-0.5 1.0]); Eep=pop_epoch(EEG,tAe,[-0.5 1.0]); Lep=pop_epoch(LFPc,tAl,[-0.5 1.0]);
    ns=min([size(Mep.data,2) size(Eep.data,2) size(Lep.data,2)]); nt=min([Mep.trials Eep.trials Lep.trials]);
    t=(0:ns-1)/fs-0.5;
    e=squeeze(mean(Mep.data(1,1:ns,:),3)); envr(end+1)=mean(e(t>=0&t<0.5))/mean(e(t>=-0.5&t<0)); [~,ix]=max(e); envpk(end+1)=t(ix);
    labE=lower({EEG.chanlocs.labels}); m1i=[]; for c=M1set,h=find(strcmp(labE,lower(c{1})),1);if ~isempty(h),m1i(end+1)=h;end;end
    zr=@(B)(B-mean(B,2))./max(std(B,0,2),eps);
    m1=zeros(1,ns,nt);for tr=1:nt,m1(1,:,tr)=mean(zr(squeeze(Eep.data(m1i,1:ns,tr))),1);end
    accM1=cat(3,accM1,m1); accSTN=cat(3,accSTN,Lep.data(1,1:ns,1:nt));
    for m=1:3,accE{m}=cat(3,accE{m},Mep.data(m,1:ns,1:nt));end
  end
  fprintf('%-7s %10.1f %10.2f %8.2f\n',subj,mean(offs),mean(envr),mean(envpk));
  ns=size(accM1,2); t=(0:ns-1)/fs-0.5; mask=t>=0&t<1.0;
  for m=1:3
    Xs=[accSTN(1,mask,:);accE{m}(1,mask,:)];[S,ff]=cs_csd_multitaper(Xs,fs,NW,K);[C,fhh]=cs_coherence(S,ff);SE{m}=[SE{m},squeeze(C(1,2,:))];if isempty(fh),fh=fhh;end
    Xm=[accM1(1,mask,:);accE{m}(1,mask,:)];[S,ff]=cs_csd_multitaper(Xm,fs,NW,K);C=cs_coherence(S,ff);ME{m}=[ME{m},squeeze(C(1,2,:))];
  end
  pr=[1 2;1 3;2 3];
  for p=1:3,X=[accE{pr(p,1)}(1,mask,:);accE{pr(p,2)}(1,mask,:)];[S,ff]=cs_csd_multitaper(X,fs,NW,K);C=cs_coherence(S,ff);IM{p}=[IM{p},squeeze(C(1,2,:))];end
end
ci=1-0.05^(1/(K*27-1));
fprintf('\napprox CI=%.4f\n',ci);
blow=fh>=2&fh<=12; bb=fh>=13&fh<=30; bhi=fh>=30&fh<=45;
fprintf('\n(native events) grand-avg coherence PEAK  low(2-12)/beta(13-30)/hibeta(30-45):\n');
for m=1:3,v=mean(SE{m},2);fprintf('  STN-%-4s: %.3f / %.3f / %.3f\n',ML{m},max(v(blow)),max(v(bb)),max(v(bhi)));end
for m=1:3,v=mean(ME{m},2);fprintf('  M1 -%-4s: %.3f / %.3f / %.3f\n',ML{m},max(v(blow)),max(v(bb)),max(v(bhi)));end
pl={'IOD-Tri','IOD-Del','Tri-Del'};
for p=1:3,v=mean(IM{p},2);fprintf('  %-8s: %.3f / %.3f / %.3f\n',pl{p},max(v(blow)),max(v(bb)),max(v(bhi)));end
save(fullfile(base,'RESULTS_final','CortexSTNMuscle','diag_native.mat'),'SE','ME','IM','fh','ML','pl');
