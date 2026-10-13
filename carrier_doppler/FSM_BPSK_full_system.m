%% REAL DSSS-BPSK: FREQUENCY-DOMAIN FSM carrier and Doppler-rate estimation
% Self-contained, MATLAB R2018b+, Signal Processing Toolbox.
% CORE: frame FFT -> shifted spectral cross-products Xp*conj(Xm) ->
% COMPLEX frame average -> frequency smoothing -> cyclic spectrum metric.
% No signal squaring and no decision-directed data supplied to receiver.
% Nonzero rate uses polynomial-phase compensation before the frame FFTs.
% Default known-PN preprocessing; exact sync baseline and separate mismatch
% scenarios. These are not a full unknown-code acquisition receiver.
clear; clc; close all; rng(20261014,'twister');

%% Configuration: full is a long statistical experiment; quick is a smoke run
run_mode='full';                 % 'full': 1000/condition; 'quick': 10/condition
cfg.Rb=1e5; cfg.Rc=31e5; cfg.fs=24.8e6; cfg.fc=1e5;
cfg.T=.1; cfg.scfFs=8e5; cfg.L=128; cfg.hop=64;
cfg.confirmCandidates=5;          % exact FSM confirmation before local refinement
cfg.smoothingBins=17;            % spectral smoothing width, 17*6250 Hz
cfg.fRange=[50e3 150e3];         % fixed received-frequency prior, NO truth input
cfg.rateRange=[-1000 1000]; cfg.rateStep=100;
cfg.fTolerance=20; cfg.rateTolerance=200;
cfg.fdTest=45000; cfg.rateTest=900;
cfg.useKnownPN=true;             % false: separate RAW FSM baseline, no PN wipeoff
cfg.rxPhaseChip=0; cfg.rxRatePPM=0; cfg.jerk=0;
cfg.rawDisplaySeconds=60e-6;
EbN0_dB=[-10 -8 -7 -6 -4];
duration_list=[.001 .01 .05 .1];
fd_edges=[-45000 -30000 -15000 0 15000 30000 45000];
rate_edges=[-900 -600 -300 0 300 600 900];
if strcmp(run_mode,'full'); numMonte=1000; else; numMonte=10; end
numAux=numMonte;
assert(any(strcmp(run_mode,{'full','quick'})));
assert(exist('resample','file')~=0,'Signal Processing Toolbox required.');
assert(cfg.fs/cfg.scfFs==round(cfg.fs/cfg.scfFs));
assert(mod(cfg.smoothingBins,2)==1 && cfg.smoothingBins<cfg.L);
assert(cfg.fRange(2)+(cfg.smoothingBins-1)/2*cfg.scfFs/cfg.L<cfg.scfFs/2);
code=make_gold();
outputDir=fullfile(pwd,['FSM_frequency_results_' datestr(now,'yyyymmdd_HHMMSS')]);
baseOutput=outputDir; suffix=1;
while exist(outputDir,'dir')
    outputDir=[baseOutput '_' num2str(suffix)]; suffix=suffix+1;
end
mkdir(outputDir);
fprintf('Mode %s: %d independent trials/condition; long runs may take hours.\n',run_mode,numMonte);
fprintf('RAW fs %.2f MHz; FSM fs %.2f MHz; FFT df %.1f Hz; frame %.3f ms.\n', ...
    cfg.fs/1e6,cfg.scfFs/1e6,cfg.scfFs/cfg.L,1e3*cfg.L/cfg.scfFs);
fprintf('-7 dB Eb/N0 -> real sample SNR %.3f dB.\n',-7+10*log10(2*cfg.Rb/cfg.fs));
fprintf('Core is rate-compensated ordinary FSM cyclic spectrum; no BPSK/QPSK classifier.\n');

%% A. Noiseless range AND boundary scan, off-grid frequencies and rates
scan_f=[50000 50037.7 100173.7 149963.3 150000];
scan_a=[-1000 -873.1 0 731.4 1000];
scan=[];
for f=scan_f
    for a=scan_a
        [row,~]=one_trial(cfg,code,Inf,f-cfg.fc,a,false);
        scan=[scan;row]; %#ok<AGROW>
    end
end
scan_table=rows_table(scan); writetable(scan_table,fullfile(outputDir,'noiseless_scan.csv'));
fprintf('Noiseless scan maximum final |f error| %.4f Hz, |rate error| %.4f Hz/s.\n', ...
    max(abs(scan(:,9))),max(abs(scan(:,10))));

%% B. Independent Eb/N0 sweep: preserve every coarse and fine result
all_rows=cell(numel(EbN0_dB),1); summary_rows=[]; demo=struct();
for k=1:numel(EbN0_dB)
    rows=zeros(numMonte,13);
    for mc=1:numMonte
        fd=(2*rand-1)*cfg.fdTest; a=(2*rand-1)*cfg.rateTest;
        keep=EbN0_dB(k)==-7 && mc==1;
        [rows(mc,:),detail]=one_trial(cfg,code,EbN0_dB(k),fd,a,keep);
        if keep; demo=detail; end
        if mod(mc,25)==0 || mc==numMonte
            fprintf('Eb/N0 %g: %d/%d\n',EbN0_dB(k),mc,numMonte);
        end
    end
    all_rows{k}=rows;
    writetable(rows_table(rows),fullfile(outputDir,sprintf('trials_EbN0_%g.csv',EbN0_dB(k))));
    summary_rows=[summary_rows;EbN0_dB(k) summarize(rows,cfg)]; %#ok<AGROW>
    % Checkpoint every condition; interrupted runs retain completed results.
    save(fullfile(outputDir,'checkpoint.mat'),'cfg','EbN0_dB','all_rows','summary_rows','scan','demo');
end
summary=summary_table(summary_rows,'EbN0_dB'); disp(summary);
writetable(summary,fullfile(outputDir,'EbN0_summary.csv'));
j=find(EbN0_dB==-7,1); target_rows=all_rows{j};

%% C. Group by TRUE frequency offset/rate; include counts AND Wilson intervals
frequency_groups=group_stats(target_rows,target_rows(:,1)-cfg.fc,fd_edges,cfg);
rate_groups=group_stats(target_rows,target_rows(:,2),rate_edges,cfg);
writetable(frequency_groups,fullfile(outputDir,'frequency_groups_minus7.csv'));
writetable(rate_groups,fullfile(outputDir,'rate_groups_minus7.csv'));

%% D. Observation duration comparison; reuse the identical 100 ms baseline
length_rows=[];
for duration=duration_list
    c=cfg; c.T=duration;
    if duration==cfg.T
        rows=target_rows;
    else
        rows=zeros(numAux,13);
        for mc=1:numAux
            [rows(mc,:),~]=one_trial(c,code,-7,(2*rand-1)*c.fdTest,(2*rand-1)*c.rateTest,false);
        end
        writetable(rows_table(rows),fullfile(outputDir,sprintf('duration_%g_ms_trials.csv',1000*duration)));
    end
    length_rows=[length_rows;duration*1000 summarize(rows,c)]; %#ok<AGROW>
    fprintf('Duration %.1f ms complete.\n',duration*1000);
end
length_summary=summary_table(length_rows,'Duration_ms');
writetable(length_summary,fullfile(outputDir,'duration_summary.csv'));

%% E. One-factor-at-a-time sensitivity at -7 dB; never pool with ideal baseline
% Relative receiver-code phase, code-rate mismatch, nonlinear phase model.
scenario_names={'IdealSync','Phase_0p125_chip','Phase_0p5_chip','Phase_1_chip', ...
    'Rate_1_ppm','Rate_5_ppm','Rate_25_ppm','Jerk_1000','Jerk_5000','RawFSM_BandlimitedNoPN'};
phase_chip=[0 .125 .5 1 0 0 0 0 0 0];
rate_ppm=[0 0 0 0 1 5 25 0 0 0];
jerk_values=[0 0 0 0 0 0 0 1000 5000 0];
scenario_values=[]; scenario_data=cell(numel(scenario_names),1);
for s=1:numel(scenario_names)
    c=cfg; c.rxPhaseChip=phase_chip(s); c.rxRatePPM=rate_ppm(s); c.jerk=jerk_values(s);
    if s==numel(scenario_names); c.useKnownPN=false; end
    if s==1; rows=target_rows;
    else
        rows=zeros(numAux,13);
        for mc=1:numAux
            [rows(mc,:),~]=one_trial(c,code,-7,(2*rand-1)*c.fdTest,(2*rand-1)*c.rateTest,false);
        end
        writetable(rows_table(rows),fullfile(outputDir,[scenario_names{s} '_trials.csv']));
    end
    scenario_data{s}=rows;
    scenario_values=[scenario_values;s summarize(rows,c)]; %#ok<AGROW>
    fprintf('Sensitivity %s complete.\n',scenario_names{s});
end
scenario_summary=summary_table(scenario_values,'ScenarioIndex');
scenario_summary.Scenario=scenario_names(:);
writetable(scenario_summary,fullfile(outputDir,'sensitivity_summary.csv'));

%% F. Figures: data/channel, ACTUAL cyclic spectrum, coarse/fine statistics
figs=gobjects(9,1);
figs(1)=figure('Name','Signal chain: random data to real AWGN received','Color','w');
subplot(3,2,1); stairs(demo.t_us,demo.data); title('Random data'); ylim([-1.3 1.3]); grid on;
subplot(3,2,2); stairs(0:30,code); title('Known Gold chips'); grid on;
subplot(3,2,3); plot(demo.t_us,demo.spread); title('Spread and sampled +/-1'); grid on;
subplot(3,2,4); plot(demo.t_us,demo.tx); title('Transmitted real BPSK'); grid on;
subplot(3,2,5); plot(demo.t_us,demo.clean); title('Doppler channel before noise'); grid on;
subplot(3,2,6); plot(demo.t_us,demo.rx); title('-7 dB Eb/N0 (not sample SNR)'); grid on;
for ax=findall(figs(1),'Type','axes').'; xlabel(ax,'Time (us), except chip panel'); end

figs(2)=figure('Name','Ordinary FSM cyclic spectrum and frequency smoothing','Color','w');
D=demo.est.detail; local_alpha=D.alpha_axis;
subplot(2,2,[1 2]); surf(D.f_axis/1e3,local_alpha/1e3,abs(D.Smap),'EdgeColor','none');
xlabel('Spectral frequency f (kHz)'); ylabel('Cyclic frequency alpha=2fc (kHz)');
zlabel('|S_x^{alpha}(f)| (power/Hz)'); title('FFT cross-products, complex frame averaging; rate compensated');
view(45,35); colorbar; grid on;
subplot(2,2,3); plot(D.f_axis/1e3,abs(D.Sfine),'LineWidth',1.2); hold on;
plot(D.f_axis/1e3,abs(D.SfineSmoothed),'LineWidth',1.2);
xlabel('Spectral frequency f (kHz)'); ylabel('Cyclic spectrum magnitude');
legend('Before frequency smoothing','After frequency smoothing'); grid on;
subplot(2,2,4); plot(local_alpha/2,abs(D.Smap(:,cfg.L/2+1)).^2,'o-');
hold on; xline(demo.truthF,'k--'); xline(demo.est.fFine,'r--');
xlabel('Received carrier candidate alpha/2 (Hz)'); ylabel('FSM score |S(0)|^2'); grid on;

figs(3)=figure('Name','Coarse joint FSM search and fitted Doppler','Color','w');
subplot(2,2,[1 2]); surf(D.coarseF/1e3,D.coarseRates,D.coarseMetric,'EdgeColor','none');
xlabel('Received carrier (kHz)'); ylabel('Rate (Hz/s)'); zlabel('Normalized coarse FSM score');
title('Global fixed-prior search, display zoom near selected peak'); view(45,35); colorbar;
subplot(2,2,3); tl=linspace(0,cfg.T,500); tr=tl-demo.rawReference;
ftrue=cfg.fc+demo.fd+demo.rate*tr+.5*cfg.jerk*tr.^2;
fhat=demo.est.fFine+demo.est.aFine*(tl-demo.est.reference);
plot(tl*1000,ftrue-cfg.fc,'k',tl*1000,fhat-cfg.fc,'r--');
xlabel('Time (ms)'); ylabel('Doppler (Hz)'); legend('True','Estimated linear model'); grid on;
subplot(2,2,4); plot(tl*1000,fhat-ftrue); xlabel('Time (ms)'); ylabel('Model frequency error (Hz)'); grid on;

figs(4)=figure('Name','Coarse vs fine success and confidence','Color','w');
subplot(1,2,1); errorbar(summary.EbN0_dB,100*summary.CoarseJoint, ...
    100*(summary.CoarseJoint-summary.CoarseLo),100*(summary.CoarseHi-summary.CoarseJoint),'o-'); hold on;
errorbar(summary.EbN0_dB,100*summary.FineJoint,100*(summary.FineJoint-summary.FineLo), ...
    100*(summary.FineHi-summary.FineJoint),'s-'); yline(90,'k--'); ylim([0 105]);
xlabel('Information Eb/N0 (dB)'); ylabel('Joint success (%)'); legend('Coarse','Fine','90% target'); grid on;
subplot(1,2,2); plot(summary.EbN0_dB,100*summary.StrictFrequencyWin,'o-', ...
    summary.EbN0_dB,100*summary.StrictRateWin,'s-');
xlabel('Eb/N0 (dB)'); ylabel('Strict fine improvement (%)');
legend('|Fine f error| < |Coarse f error|','|Fine rate error| < |Coarse rate error|'); grid on;
title('This win fraction is NOT estimation success');

figs(5)=figure('Name','Coarse/fine RMSE and biases','Color','w');
subplot(2,2,1); semilogy(EbN0_dB,max(summary.CoarseF_RMSE,eps),'o-',EbN0_dB,max(summary.FineF_RMSE,eps),'s-');
ylabel('Frequency RMSE (Hz)'); legend('Coarse','Fine'); grid on;
subplot(2,2,2); semilogy(EbN0_dB,max(summary.CoarseRate_RMSE,eps),'o-',EbN0_dB,max(summary.FineRate_RMSE,eps),'s-');
ylabel('Rate RMSE (Hz/s)'); legend('Coarse','Fine'); grid on;
subplot(2,2,3); errorbar(EbN0_dB,summary.FineF_Bias,summary.FineF_Std,'o-');
ylabel('Fine f bias +/- std (Hz)'); grid on;
subplot(2,2,4); errorbar(EbN0_dB,summary.FineRate_Bias,summary.FineRate_Std,'o-');
ylabel('Fine rate bias +/- std (Hz/s)'); grid on;
for ax=findall(figs(5),'Type','axes').'; xlabel(ax,'Eb/N0 (dB)'); end

figs(6)=figure('Name','Frequency/rate stratified performance at -7 dB','Color','w');
subplot(2,2,1); plot_groups(frequency_groups); xlabel('True Doppler bin center (Hz)'); ylabel('Fine joint success (%)');
subplot(2,2,2); plot_groups(rate_groups); xlabel('True rate bin center (Hz/s)'); ylabel('Fine joint success (%)');
subplot(2,2,3); bar(frequency_groups.Center,frequency_groups.Count); xlabel('Doppler bin center (Hz)'); ylabel('Independent trials'); grid on;
subplot(2,2,4); bar(rate_groups.Center,rate_groups.Count); xlabel('Rate bin center (Hz/s)'); ylabel('Independent trials'); grid on;

figs(7)=figure('Name','Duration tradeoff and sensitivity','Color','w');
subplot(2,2,1); errorbar(length_summary.Duration_ms,100*length_summary.FineJoint, ...
    100*(length_summary.FineJoint-length_summary.FineLo),100*(length_summary.FineHi-length_summary.FineJoint),'o-');
xlabel('Observation (ms)'); ylabel('Joint success (%)'); ylim([0 105]); grid on;
subplot(2,2,2); semilogy(length_summary.Duration_ms,max(length_summary.FineF_RMSE,eps),'o-');
xlabel('Observation (ms)'); ylabel('Frequency RMSE (Hz)'); grid on;
subplot(2,2,[3 4]); bar(100*scenario_summary.FineJoint);
set(gca,'XTick',1:numel(scenario_names),'XTickLabel',scenario_names); xtickangle(35);
ylabel('Fine joint success (%)'); ylim([0 105]); grid on;

figs(8)=figure('Name','Noiseless parameter scan: finite-record bias','Color','w');
subplot(1,2,1); scatter3(scan(:,1)/1e3,scan(:,2),scan(:,9),40,scan(:,9),'filled');
xlabel('True carrier (kHz)'); ylabel('True rate (Hz/s)'); zlabel('Fine f bias (Hz)'); grid on; colorbar;
subplot(1,2,2); scatter3(scan(:,1)/1e3,scan(:,2),scan(:,10),40,scan(:,10),'filled');
xlabel('True carrier (kHz)'); ylabel('True rate (Hz/s)'); zlabel('Fine rate bias (Hz/s)'); grid on; colorbar;

figs(9)=figure('Name','Independent errors at -7 dB: keep failed trials','Color','w');
subplot(2,2,1); histogram(target_rows(:,7),25); hold on; histogram(target_rows(:,9),25);
xlabel('Frequency error (Hz)'); legend('Coarse','Fine'); grid on;
subplot(2,2,2); histogram(target_rows(:,8),25); hold on; histogram(target_rows(:,10),25);
xlabel('Rate error (Hz/s)'); legend('Coarse','Fine'); grid on;
subplot(2,2,3); scatter(target_rows(:,1)/1e3,target_rows(:,5)/1e3,20,'filled');
xlabel('True carrier (kHz)'); ylabel('Fine estimate (kHz)'); grid on;
subplot(2,2,4); scatter(target_rows(:,2),target_rows(:,6),20,'filled');
xlabel('True rate (Hz/s)'); ylabel('Fine estimate (Hz/s)'); grid on;
for i=1:numel(figs)
    savefig(figs(i),fullfile(outputDir,sprintf('figure_%02d.fig',i)));
    print(figs(i),fullfile(outputDir,sprintf('figure_%02d.png',i)),'-dpng','-r150');
end
save(fullfile(outputDir,'results.mat'),'cfg','run_mode','summary','scan_table','all_rows', ...
    'frequency_groups','rate_groups','length_summary','scenario_summary','scenario_data','demo','code');
fprintf('Saved completed results: %s\n',outputDir);
if summary.FineLo(j)>0.9 && strcmp(run_mode,'full')
    fprintf('This run supports >90%% under IDEAL known-code timing and stated thresholds only.\n');
else
    fprintf('This run does not establish the requested full statistical claim.\n');
end

%% Independent transmitter, channel and receiver; truth never enters estimator
function [row,demo]=one_trial(c,code,eb,fd,a,keep)
    nb=round(c.Rb*c.T); Ns=c.fs/c.Rc; N=nb*31*Ns;
    assert(Ns==round(Ns) && nb>0);
    rawTime=(0:N-1)/c.fs; reference=(N-1)/(2*c.fs); t=rawTime-reference;
    bits=2*randi([0 1],1,nb)-1;
    data=repelem(bits,31*Ns); pn=repelem(repmat(code,1,nb),Ns);
    spread=data.*pn;
    channelPhase=2*pi*fd*t+pi*a*t.^2+pi/3*c.jerk*t.^3+2*pi*rand;
    clean=spread.*cos(2*pi*c.fc*t+channelPhase);
    if isinf(eb); rx=clean;
    else
        variance=mean(clean.^2)*c.fs/(2*c.Rb*10^(eb/10));
        rx=clean+sqrt(variance)*randn(1,N);
    end
    if c.useKnownPN
        ci=floor((0:N-1)*(1+c.rxRatePPM*1e-6)/Ns+c.rxPhaseChip);
        receiverPN=code(mod(ci,31)+1);
        processed=rx.*receiverPN;
    else; processed=rx; end
    % REAL resampling, no downconversion, no extra 90 kHz cutoff.
    y=resample(processed(:),1,round(c.fs/c.scfFs));
    est=fsm_estimate(y,c,keep);
    delta=est.reference-reference;
    trueF=c.fc+fd+a*delta+.5*c.jerk*delta^2;
    trueA=a+c.jerk*delta;
    row=[trueF,trueA,est.fCoarse,est.aCoarse,est.fFine,est.aFine, ...
        est.fCoarse-trueF,est.aCoarse-trueA,est.fFine-trueF,est.aFine-trueA, ...
        est.converged,est.atBoundary,est.frameCount];
    demo=struct();
    if keep
        ns=min(N,round(c.rawDisplaySeconds*c.fs));
        demo=struct('t_us',rawTime(1:ns)*1e6,'data',data(1:ns),'spread',spread(1:ns), ...
            'tx',spread(1:ns).*cos(2*pi*c.fc*t(1:ns)),'clean',clean(1:ns),'rx',rx(1:ns), ...
            'fd',fd,'rate',a,'rawReference',reference,'truthF',trueF,'est',est);
    end
end

%% Actual FSM in the frequency domain
function est=fsm_estimate(y,c,keep)
    n=numel(y); starts=0:c.hop:n-c.L; assert(numel(starts)>=3);
    indices=bsxfun(@plus,(0:c.L-1)',starts)+1;
    w=0.5-0.5*cos(2*pi*(0:c.L-1)'/(c.L-1));
    frames=bsxfun(@times,y(indices),w);
    ref=(n-1)/(2*c.scfFs);
    tc=(indices-1)/c.scfFs-ref;
    ts=starts/c.scfFs-ref;
    band=-(c.smoothingBins-1)/2:(c.smoothingBins-1)/2;
    rates=unique([c.rateRange(1):c.rateStep:c.rateRange(2),c.rateRange(2)]);
    df=c.scfFs/c.L; X=fft(frames,[],1);
    M=numel(starts); nt=2^nextpow2(2*M);
    deltaAlpha=(-nt/2:nt/2-1)*(c.scfFs/c.hop)/nt;
    best=-Inf; f0=NaN; a0=NaN;
    coarseF=[]; coarseScore=[]; shortlist=[];
    % Bank of S^{alpha_b}(f) spectral cross-products. A temporal FFT fills
    % the cyclic-frequency interval around each alpha_b=2*k*df.
    for k=floor(c.fRange(1)/df):ceil(c.fRange(2)/df)
        fb=k*df;
        p=mean(X(mod(band+k,c.L)+1,:).*conj(X(mod(band-k,c.L)+1,:)),1);
        p=p.*exp(-1j*4*pi*fb*ts);
        bank=bsxfun(@times,p,exp(-1j*2*pi*rates(:)*ts.^2));
        power=abs(fftshift(fft(bank,nt,2),2)).^2;
        candidates=fb+deltaAlpha/2;
        valid=deltaAlpha>=-df & deltaAlpha<df & candidates>=c.fRange(1) & candidates<=c.fRange(2);
        ids=find(valid); if isempty(ids); continue; end
        P=power(:,valid); [value,linear]=max(P(:)); [ir,jf]=ind2sub(size(P),linear);
        shortlist=[shortlist;value,candidates(ids(jf)),rates(ir)]; %#ok<AGROW>
        if value>best; best=value; f0=candidates(ids(jf)); a0=rates(ir); end
        if keep
            coarseF=[coarseF,candidates(valid)]; coarseScore=[coarseScore,P]; %#ok<AGROW>
        end
    end
    % Short-time cyclic sampling can produce adjacent FFT-branch aliases.
    % Rescore multiple candidates using exact FSM; never use ground truth.
    [~,ranking]=sort(shortlist(:,1),'descend');
    ranking=ranking(1:min(c.confirmCandidates,numel(ranking)));
    verified=-Inf;
    for j=ranking(:).'
        value=fsm_score(frames,tc,shortlist(j,2),shortlist(j,3),c,w,band);
        if value>verified
            verified=value; f0=shortlist(j,2); a0=shortlist(j,3);
        end
    end
    T=n/c.scfFs;
    initial=max(fsm_score(frames,tc,f0,a0,c,w,band),realmin);
    fun=@(p) exact_objective(p,frames,tc,f0,a0,T,c,w,band,initial);
    [p,~,flag]=fminsearch(fun,[0 0],optimset('Display','off','TolX',1e-5, ...
        'TolFun',1e-8,'MaxIter',200));
    est.fCoarse=f0; est.aCoarse=a0; est.fFine=f0+p(1)/T; est.aFine=a0+p(2)/T^2;
    est.reference=ref; est.converged=flag>0; est.frameCount=M;
    est.atBoundary=min(abs(est.fFine-c.fRange))<1/T || ...
        min(abs(est.aFine-c.rateRange))<1/T^2;
    if keep
        faxis=(-c.L/2:c.L/2-1)'*df;
        S=fsm_spectrum(frames,tc,est.fFine,est.aFine,c,w);
        kernel=ones(c.smoothingBins,1)/c.smoothingBins;
        fineF=linspace(max(c.fRange(1),est.fFine-50),min(c.fRange(2),est.fFine+50),41);
        Smap=zeros(numel(fineF),c.L);
        for j=1:numel(fineF)
            spec=fsm_spectrum(frames,tc,fineF(j),est.aFine,c,w);
            Smap(j,:)=conv(spec,kernel,'same').';
        end
        [coarseF,order]=sort(coarseF); coarseScore=coarseScore(:,order);
        zoom=abs(coarseF-f0)<100;
        est.detail=struct('f_axis',faxis,'alpha_axis',2*fineF,'Smap',Smap, ...
            'Sfine',S,'SfineSmoothed',conv(S,kernel,'same'), ...
            'coarseF',coarseF(zoom),'coarseRates',rates, ...
            'coarseMetric',coarseScore(:,zoom)/max(coarseScore(:)));
    end
end

function spec=fsm_spectrum(frames,tc,f,a,c,w)
    % alpha=2*f; polynomial phase accounts for cyclic-frequency drift 2*a.
    % At a=0, exactly the ordinary shifted-FFT FSM cyclic spectrum.
    phase=2*pi*f*tc+pi*a*tc.^2;
    Xp=fft(frames.*exp(-1j*phase),[],1);
    Xm=fft(frames.*exp( 1j*phase),[],1);
    spec=fftshift(mean(Xp.*conj(Xm),2))/(c.scfFs*sum(w.^2));
end

function J=fsm_score(frames,tc,f,a,c,w,band)
    spec=fsm_spectrum(frames,tc,f,a,c,w);
    % Frequency-smoothed complex SCF at spectral f=0; absolute value last.
    J=abs(mean(spec(c.L/2+1+band)))^2;
end

function value=exact_objective(p,frames,tc,f0,a0,T,c,w,band,normalizer)
    f=f0+p(1)/T; a=a0+p(2)/T^2;
    if f<c.fRange(1)||f>c.fRange(2)||a<c.rateRange(1)||a>c.rateRange(2)
        value=1e6; return;
    end
    value=-fsm_score(frames,tc,f,a,c,w,band)/normalizer;
end

function code=make_gold()
    a=mseq([0 2]); b=mseq([0 2 3 4]);
    code=2*double(xor(logical(a),circshift(logical(b),[0 1])))-1;
end
function bits=mseq(taps)
    s=ones(1,36); for j=1:31; s(j+5)=mod(sum(s(j+taps)),2); end
    assert(isequal(s(32:36),s(1:5)) && sum(s(1:31))==16); bits=s(1:31);
end
function t=rows_table(r)
    t=array2table(r,'VariableNames',{'TrueF_Hz','TrueRate_HzPerS','CoarseF_Hz', ...
        'CoarseRate_HzPerS','FineF_Hz','FineRate_HzPerS','CoarseF_Error', ...
        'CoarseRate_Error','FineF_Error','FineRate_Error','Converged','AtBoundary','FrameCount'});
end
function s=summarize(r,c)
    cf=abs(r(:,7))<=c.fTolerance; ca=abs(r(:,8))<=c.rateTolerance;
    ff=abs(r(:,9))<=c.fTolerance; fa=abs(r(:,10))<=c.rateTolerance;
    cj=cf&ca; fj=ff&fa&r(:,11)>0;
    [cl,ch]=wilson(sum(cj),size(r,1)); [fl,fh]=wilson(sum(fj),size(r,1));
    s=[size(r,1),mean(cj),cl,ch,mean(fj),fl,fh,mean(cf),mean(ca),mean(ff),mean(fa), ...
        sqrt(mean(r(:,7).^2)),sqrt(mean(r(:,8).^2)),sqrt(mean(r(:,9).^2)),sqrt(mean(r(:,10).^2)), ...
        mean(r(:,9)),mean(r(:,10)),std(r(:,9)),std(r(:,10)), ...
        mean(abs(r(:,9))<abs(r(:,7))),mean(abs(r(:,10))<abs(r(:,8))),mean(r(:,11)>0),mean(r(:,12)>0), ...
        mean(r(:,7)),mean(r(:,8)),std(r(:,7)),std(r(:,8))];
end
function t=summary_table(r,firstName)
    names={firstName,'Count','CoarseJoint','CoarseLo','CoarseHi','FineJoint','FineLo','FineHi', ...
        'CoarseF_Success','CoarseRate_Success','FineF_Success','FineRate_Success', ...
        'CoarseF_RMSE','CoarseRate_RMSE','FineF_RMSE','FineRate_RMSE', ...
        'FineF_Bias','FineRate_Bias','FineF_Std','FineRate_Std','StrictFrequencyWin', ...
        'StrictRateWin','ConvergenceFraction','BoundaryFraction','CoarseF_Bias','CoarseRate_Bias', ...
        'CoarseF_Std','CoarseRate_Std'};
    t=array2table(r,'VariableNames',names);
end
function t=group_stats(r,x,edges,c)
    values=nan(numel(edges)-1,12);
    for j=1:numel(edges)-1
        use=x>=edges(j)&x<edges(j+1);
        if j==numel(edges)-1; use=x>=edges(j)&x<=edges(j+1); end
        n=sum(use); center=(edges(j)+edges(j+1))/2;
        if n==0; values(j,:)=[edges(j),edges(j+1),center,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]; continue; end
        rr=r(use,:); ok=abs(rr(:,9))<=c.fTolerance & abs(rr(:,10))<=c.rateTolerance & rr(:,11)>0;
        [lo,hi]=wilson(sum(ok),n);
        coarse_ok=abs(rr(:,7))<=c.fTolerance & abs(rr(:,8))<=c.rateTolerance;
        values(j,:)=[edges(j),edges(j+1),center,n,mean(ok),lo,hi,sqrt(mean(rr(:,9).^2)), ...
            sqrt(mean(rr(:,10).^2)),mean(rr(:,9)),mean(rr(:,10)),mean(coarse_ok)];
    end
    t=array2table(values,'VariableNames',{'Lower','Upper','Center','Count','FineJoint','Lo','Hi','FrequencyRMSE','RateRMSE','FrequencyBias','RateBias','CoarseJoint'});
end
function plot_groups(t)
    errorbar(t.Center,100*t.FineJoint,100*(t.FineJoint-t.Lo),100*(t.Hi-t.FineJoint),'o-');
    hold on; yline(90,'k--'); ylim([0 105]); grid on;
end
function [lo,hi]=wilson(success,n)
    z=1.96; p=success/n; den=1+z^2/n;
    center=(p+z^2/(2*n))/den; half=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/den;
    lo=max(0,center-half); hi=min(1,center+half);
end
