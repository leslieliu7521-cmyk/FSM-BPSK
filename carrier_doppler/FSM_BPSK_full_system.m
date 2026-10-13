%% Full real-passband DSSS-BPSK carrier/Doppler simulation
% MATLAB R2018b+ and Signal Processing Toolbox (resample, fir1).
% Self-contained: no external mseq.m or estimator file is required.
% Known PN, exact code/bit timing; unknown random data; AWGN; linear Doppler.
% Includes NO code-rate Doppler, multipath, jerk or oscillator drift.
% 100 ms results apply to these assumptions, not complete acquisition.
% The included Gold pair uses an explicit convention; it does not claim
% to reproduce the missing original mseq.m convention.
clear; clc; close all; rng(20261014,'twister');

%% 1. User parameters
Rb=1e5; Rc=31*Rb; fs=24.8e6; Ns=fs/Rc;
Tobs=0.1;                         % 100 ms; do not merely extend a plotting axis
fc_nominal=1e5;                   % independently known receiver IF, Hz
EbN0_dB=[-10 -7 -4];             % INFORMATION-bit Eb/N0, not sample SNR
numMonte=100;                    % independent bits, phase, Doppler and noise
plot_EbN0=-7;
fd_test_max=45000;               % random midpoint Doppler within +/-45 kHz
rate_test_max=900;               % random linear frequency rate within +/-900 Hz/s
frequency_tolerance_hz=20;
rate_tolerance_hz_per_s=200;
saveResults=true;               % MAT, CSV, PNG and FIG in timestamped directory
cfg.outFs=4e5; cfg.bandwidthHz=9e4;
cfg.fRange=[-5e4 5e4];           % fixed receiver prior, never centered at truth
cfg.rateRange=[-1000 1000];
cfg.rateStep=100;                % suitable scale <=1/Tobs^2; refine continuously
assert(exist('resample','file')~=0 && exist('fir1','file')~=0, ...
    'Install/enable Signal Processing Toolbox.');
assert(Ns==round(Ns) && fs/Rb==round(fs/Rb));
assert(abs(Rb*Tobs-round(Rb*Tobs))<1e-8 && Tobs>=0.001);
assert(any(EbN0_dB==plot_EbN0));
assert(fd_test_max<max(cfg.fRange) && rate_test_max<max(cfg.rateRange));
assert(fc_nominal-fd_test_max-rate_test_max*Tobs/2>0);
assert(cfg.bandwidthHz<cfg.outFs/2);

%% 2. Gold spreading code, chip sampling and time convention
% Primitive preferred pair x^5+x^2+1 and x^5+x^4+x^3+x^2+1.
% Recurrence/seed explicitly defined below; shift convention is explicit.
m1=local_mseq([0 2]); m2=local_mseq([0 2 3 4]);
gold_bits=xor(logical(m1),circshift(logical(m2),[0 1]));
gold=2*double(gold_bits)-1; L=numel(gold);
assert(L==31 && Rc==L*Rb);
numBits=round(Rb*Tobs); N=numBits*L*Ns;
t0=(0:N-1)/fs; record_mid=(N-1)/(2*fs); tc=t0-record_mid;
% Exactly synchronized receiver PN contains NO information data.
pn_chips=repmat(gold,1,numBits);
pn_sample=repelem(pn_chips,Ns);
sample_snr_dB=EbN0_dB+10*log10(2*Rb/fs);
fprintf('N=%d, T=%.3f ms, bits=%d, samples/chip=%d\n',N,1e3*Tobs,numBits,Ns);
fprintf('Eb/N0=-7 dB corresponds to sample SNR %.3f dB (real AWGN).\n', ...
    -7+10*log10(2*Rb/fs));
fprintf('Independent MC: %d runs per Eb/N0; this may take several minutes.\n',numMonte);

%% 3. Noiseless diagnostic, including off-grid truth
[~,~,~,clean0]=transmit_record(numBits,gold,Ns,tc,fc_nominal,17321.7,413.9,0.83);
e0=receiver_estimate(clean0,pn_sample,fs,fc_nominal,cfg,false);
truth0=17321.7+413.9*(e0.reference_time_s-record_mid);
fprintf('Noiseless check: frequency bias %.6f Hz; rate bias %.6f Hz/s\n', ...
    e0.doppler_mid_hz-truth0,e0.rate_hz_per_s-413.9);
clear clean0;

%% 4. Full independent Monte Carlo channel and receiver
K=numel(EbN0_dB); error_frequency=zeros(K,numMonte);
error_rate=zeros(K,numMonte); true_frequency=zeros(K,numMonte);
true_rate=zeros(K,numMonte); hat_frequency=zeros(K,numMonte);
hat_rate=zeros(K,numMonte); solver_ok=false(K,numMonte);
demo=struct();
for k=1:K
    for mc=1:numMonte
        fd_mid=(2*rand-1)*fd_test_max;
        rate=(2*rand-1)*rate_test_max;
        phase=2*pi*rand;
        % Random bits -> spreading -> rectangular chip sampling -> REAL BPSK.
        % Phase includes Doppler from the propagation channel. Equivalent to
        % received passband cos(2*pi*(fc+fd_mid)*tc+pi*rate*tc^2+phase).
        [bits,spread_chips,spread_sample,clean]=transmit_record( ...
            numBits,gold,Ns,tc,fc_nominal,fd_mid,rate,phase);
        signal_power=mean(clean.^2);
        % Real AWGN with two-sided PSD N0/2; Eb=signal_power/Rb.
        noise_variance=signal_power*fs/(2*Rb*10^(EbN0_dB(k)/10));
        noise=sqrt(noise_variance)*randn(1,N);
        rx=clean+noise;
        keep=(EbN0_dB(k)==plot_EbN0 && mc==1);
        est=receiver_estimate(rx,pn_sample,fs,fc_nominal,cfg,keep);
        % Always compare frequency at the same observation reference time.
        fd_at_reference=fd_mid+rate*(est.reference_time_s-record_mid);
        true_frequency(k,mc)=fc_nominal+fd_at_reference;
        true_rate(k,mc)=rate;
        hat_frequency(k,mc)=est.carrier_mid_hz;
        hat_rate(k,mc)=est.rate_hz_per_s;
        error_frequency(k,mc)=est.doppler_mid_hz-fd_at_reference;
        error_rate(k,mc)=est.rate_hz_per_s-rate;
        solver_ok(k,mc)=est.optimizer_converged;
        if keep
            demo.bits=bits(1:min(12,numBits)); demo.gold=gold;
            nshow=min(N,round(60e-6*fs));
            demo.t_us=t0(1:nshow)*1e6;
            demo.data_samples=bits(floor((0:nshow-1)/(fs/Rb))+1);
            demo.spread_chips=spread_chips(1:min(186,numel(spread_chips)));
            demo.pn=pn_sample(1:nshow); demo.spread=spread_sample(1:nshow);
            demo.tx=spread_sample(1:nshow).*cos(2*pi*fc_nominal*tc(1:nshow));
            demo.clean=clean(1:nshow); demo.noise=noise(1:nshow); demo.rx=rx(1:nshow);
            demo.fd_mid=fd_mid; demo.rate=rate; demo.phase=phase;
            demo.noise_variance=noise_variance;
            demo.measured_snr_dB=10*log10(signal_power/mean(noise.^2));
            demo.est=est; demo.record_mid=record_mid;
        end
        if mod(mc,10)==0 || mc==numMonte
            fprintf('Eb/N0=%g dB, trial %d/%d\n',EbN0_dB(k),mc,numMonte);
        end
    end
end
clear rx clean noise spread_sample spread_chips bits;

%% 5. Success, RMSE, bias, standard deviation and confidence intervals
pass_f=abs(error_frequency)<=frequency_tolerance_hz;
pass_r=abs(error_rate)<=rate_tolerance_hz_per_s;
% Optimizer failure is counted as failure, not silently removed.
pass_joint=pass_f & pass_r & solver_ok;
rate_f=mean(pass_f,2); rate_r=mean(pass_r,2); rate_joint=mean(pass_joint,2);
[ci_low,ci_high]=wilson_interval(sum(pass_joint,2),numMonte);
rmse_f=sqrt(mean(error_frequency.^2,2)); rmse_r=sqrt(mean(error_rate.^2,2));
bias_f=mean(error_frequency,2); bias_r=mean(error_rate,2);
std_f=std(error_frequency,0,2); std_r=std(error_rate,0,2);
summary_table=table(EbN0_dB(:),sample_snr_dB(:),rate_f,rate_r,rate_joint, ...
    ci_low,ci_high,rmse_f,rmse_r,bias_f,bias_r,std_f,std_r,sum(solver_ok,2), ...
    'VariableNames',{'EbN0_dB','SampleSNR_dB','CarrierSuccess','RateSuccess', ...
    'JointSuccess','JointCI95Lower','JointCI95Upper','CarrierRMSE_Hz', ...
    'RateRMSE_HzPerS','CarrierBias_Hz','RateBias_HzPerS','CarrierStd_Hz', ...
    'RateStd_HzPerS','ConvergedTrials'});
disp(summary_table);
j=find(EbN0_dB==plot_EbN0,1);
fprintf('At %g dB Eb/N0: joint success %.2f%%; 95%% CI [%.2f%% %.2f%%].\n', ...
    plot_EbN0,100*rate_joint(j),100*ci_low(j),100*ci_high(j));
fprintf('Acceptance: carrier <=%g Hz AND rate <=%g Hz/s; nonconvergence fails.\n', ...
    frequency_tolerance_hz,rate_tolerance_hz_per_s);
if ci_low(j)>0.9
    fprintf('This run supports >90%% under the configured model and thresholds.\n');
else
    fprintf('This run does NOT establish a confidence lower bound above 90%%.\n');
end

%% 6. Intuitive signal-chain figures
figs=gobjects(7,1);
figs(1)=figure('Name','1 Random data - spreading - sampling - BPSK','Color','w');
subplot(3,2,1); stairs(demo.t_us,demo.data_samples,'LineWidth',1.2); ylim([-1.3 1.3]);
xlabel('Time (us)'); ylabel('Data'); title('Random information bits'); grid on;
subplot(3,2,2); stairs(0:30,demo.gold,'LineWidth',1.2); ylim([-1.3 1.3]);
xlabel('Chip index'); title('Known 31-chip Gold sequence'); grid on;
subplot(3,2,3); stairs((0:numel(demo.spread_chips)-1)/Rc*1e6,demo.spread_chips);
xlabel('Time (us)'); title('Spread chip sequence: data x PN'); grid on;
subplot(3,2,4); plot(demo.t_us,demo.spread); ylim([-1.3 1.3]);
xlabel('Time (us)'); title(sprintf('Chip sampling: %g samples/chip',Ns)); grid on;
subplot(3,2,5); plot(demo.t_us,demo.tx); xlabel('Time (us)');
ylabel('Amplitude'); title('Transmitted real BPSK at nominal carrier'); grid on;
subplot(3,2,6); plot(demo.t_us,demo.rx); xlabel('Time (us)'); ylabel('Amplitude');
title(sprintf('Received: Eb/N0=%g dB, sample SNR=%.2f dB',plot_EbN0,demo.measured_snr_dB)); grid on;

figs(2)=figure('Name','2 Channel and synchronized despreading','Color','w');
subplot(3,1,1); plot(demo.t_us,demo.tx); hold on;
plot(demo.t_us,demo.clean); plot(demo.t_us,demo.rx);
xlabel('Time (us)'); legend('Transmitted BPSK','Doppler channel (clean)','AWGN received'); grid on;
title('Real received waveform: sample noise power is much larger than signal');
subplot(3,1,2); d=demo.est.detail; show=min(numel(d.y),800);
plot(d.time_s(1:show)*1e3,real(d.y(1:show))); hold on;
plot(d.time_s(1:show)*1e3,imag(d.y(1:show)));
xlabel('Time (ms)'); legend('I','Q'); title('After PN wipeoff, downconversion and low-pass filtering'); grid on;
subplot(3,1,3); plot(d.time_s(1:show)*1e3,real(d.y(1:show).^2)); hold on;
plot(d.time_s(1:show)*1e3,imag(d.y(1:show).^2));
xlabel('Time (ms)'); legend('Squared I','Squared Q'); title('Complex squaring removes unknown +/-1 data; noise remains'); grid on;

figs(3)=figure('Name','3 True and estimated linear Doppler','Color','w');
tl=linspace(0,(N-1)/fs,500);
ftrue=demo.fd_mid+demo.rate*(tl-record_mid);
fhat=demo.est.doppler_mid_hz+demo.est.rate_hz_per_s*(tl-demo.est.reference_time_s);
subplot(2,2,1); plot(tl*1e3,(fc_nominal+ftrue)/1e3,'k','LineWidth',1.5); hold on;
plot(tl*1e3,(fc_nominal+fhat)/1e3,'r--','LineWidth',1.3);
xlabel('Time (ms)'); ylabel('Carrier (kHz)'); legend('True','Estimated linear model'); grid on;
subplot(2,2,2); plot(tl*1e3,ftrue,'k','LineWidth',1.5); hold on;
plot(tl*1e3,fhat,'r--'); xlabel('Time (ms)'); ylabel('Doppler (Hz)'); grid on;
legend('True','Estimated linear model');
subplot(2,2,3); plot(tl*1e3,fhat-ftrue,'LineWidth',1.3);
xlabel('Time (ms)'); ylabel('Frequency error (Hz)'); grid on;
subplot(2,2,4); bar([demo.rate,demo.est.rate_hz_per_s]);
set(gca,'XTickLabel',{'True','Estimated'}); ylabel('Frequency rate (Hz/s)'); grid on;
title(sprintf('True rate %.2f; estimated %.2f Hz/s',demo.rate,demo.est.rate_hz_per_s));

figs(4)=figure('Name','4 Joint frequency-rate search (not a cyclic spectrum)','Color','w');
D=demo.est.detail; normalized=D.metric/max(D.metric(:));
subplot(2,2,[1 2]); surf(D.f_axis/1e3,D.rates,normalized,'EdgeColor','none');
xlabel('Doppler (kHz)'); ylabel('Rate (Hz/s)'); zlabel('Normalized search score');
title('Coarse joint search near the selected peak'); view(45,35); colorbar; grid on;
subplot(2,2,3); plot(D.f_axis,max(normalized,[],1),'LineWidth',1.2); hold on;
fdref=demo.fd_mid+demo.rate*(demo.est.reference_time_s-record_mid);
xline(fdref,'k--','True'); xline(demo.est.doppler_mid_hz,'r--','Estimated');
xlabel('Doppler (Hz)'); ylabel('Max score over rates'); grid on;
subplot(2,2,4); plot(D.rates,max(normalized,[],2),'o-'); hold on;
xline(demo.rate,'k--','True'); xline(demo.est.rate_hz_per_s,'r--','Estimated');
xlabel('Rate (Hz/s)'); ylabel('Max score over frequencies'); grid on;

figs(5)=figure('Name','5 Estimation success and confidence','Color','w');
subplot(1,2,1); plot(EbN0_dB,100*rate_f,'o-',EbN0_dB,100*rate_r,'s-','LineWidth',1.3);
hold on; yline(90,'k--'); ylim([0 105]); xlabel('Eb/N0 (dB)'); ylabel('Success (%)');
legend('Carrier','Rate','90% target','Location','best'); grid on;
title(sprintf('Thresholds: %g Hz; %g Hz/s',frequency_tolerance_hz,rate_tolerance_hz_per_s));
subplot(1,2,2); errorbar(EbN0_dB,100*rate_joint,100*(rate_joint-ci_low), ...
    100*(ci_high-rate_joint),'o-','LineWidth',1.3); hold on; yline(90,'k--');
ylim([0 105]); xlabel('Eb/N0 (dB)'); ylabel('Joint success (%)'); grid on;
title(sprintf('%d independent trials; Wilson 95%% intervals',numMonte));

figs(6)=figure('Name','6 RMSE bias and spread','Color','w');
subplot(2,2,1); semilogy(EbN0_dB,max(rmse_f,eps),'o-');
xlabel('Eb/N0 (dB)'); ylabel('Carrier RMSE (Hz)'); grid on;
subplot(2,2,2); semilogy(EbN0_dB,max(rmse_r,eps),'o-');
xlabel('Eb/N0 (dB)'); ylabel('Rate RMSE (Hz/s)'); grid on;
subplot(2,2,3); errorbar(EbN0_dB,bias_f,std_f,'o-');
xlabel('Eb/N0 (dB)'); ylabel('Carrier bias +/- sample std (Hz)'); grid on;
subplot(2,2,4); errorbar(EbN0_dB,bias_r,std_r,'o-');
xlabel('Eb/N0 (dB)'); ylabel('Rate bias +/- sample std (Hz/s)'); grid on;

figs(7)=figure('Name','7 Independent trials at -7 dB Eb/N0','Color','w');
subplot(2,2,1); histogram(error_frequency(j,:),20); hold on;
xline(-frequency_tolerance_hz,'r--'); xline(frequency_tolerance_hz,'r--');
xlabel('Carrier error (Hz)'); ylabel('Count'); grid on;
subplot(2,2,2); histogram(error_rate(j,:),20); hold on;
xline(-rate_tolerance_hz_per_s,'r--'); xline(rate_tolerance_hz_per_s,'r--');
xlabel('Rate error (Hz/s)'); ylabel('Count'); grid on;
subplot(2,2,3); scatter(true_frequency(j,:)/1e3,hat_frequency(j,:)/1e3,24,'filled'); hold on;
limits=[min(true_frequency(j,:)),max(true_frequency(j,:))]/1e3;
plot(limits,limits,'k--'); xlabel('True carrier (kHz)'); ylabel('Estimated (kHz)'); grid on;
subplot(2,2,4); scatter(true_rate(j,:),hat_rate(j,:),24,'filled'); hold on;
plot([-rate_test_max rate_test_max],[-rate_test_max rate_test_max],'k--');
xlabel('True rate (Hz/s)'); ylabel('Estimated rate (Hz/s)'); grid on;

%% 7. Save the current run; do not substitute previous Python outcomes
if saveResults
    outputDir=fullfile(pwd,['FSM_BPSK_results_' datestr(now,'yyyymmdd_HHMMSS')]);
    if ~exist(outputDir,'dir'); mkdir(outputDir); end
    writetable(summary_table,fullfile(outputDir,'summary.csv'));
    parameters=struct('Rb',Rb,'Rc',Rc,'fs',fs,'Tobs',Tobs,'fc_nominal',fc_nominal, ...
        'EbN0_dB',EbN0_dB,'numMonte',numMonte,'frequency_tolerance_hz', ...
        frequency_tolerance_hz,'rate_tolerance_hz_per_s',rate_tolerance_hz_per_s,'cfg',cfg);
    save(fullfile(outputDir,'results.mat'),'parameters','summary_table', ...
        'error_frequency','error_rate','true_frequency','true_rate', ...
        'hat_frequency','hat_rate','solver_ok','demo','gold','e0');
    for i=1:numel(figs)
        savefig(figs(i),fullfile(outputDir,sprintf('figure_%02d.fig',i)));
        print(figs(i),fullfile(outputDir,sprintf('figure_%02d.png',i)),'-dpng','-r150');
    end
    fprintf('Saved this MATLAB run to:\n%s\n',outputDir);
end

%% Local functions
function bits=local_mseq(taps)
    % s[n+5]=XOR of s[n+tap]; degree-five recurrence, all-one seed.
    s=ones(1,36);
    for n=1:31; s(n+5)=mod(sum(s(n+taps)),2); end
    assert(isequal(s(32:36),s(1:5)),'LFSR period check failed.');
    bits=s(1:31);
    assert(sum(bits)==16,'Expected maximal-length sequence balance.');
end

function [bits,chips,samples,clean]=transmit_record(nb,gold,Ns,t,fc,fd,rate,phase)
    bits=2*randi([0 1],1,nb)-1;
    chips=repelem(bits,numel(gold)).*repmat(gold,1,nb);
    samples=repelem(chips,Ns);
    % Transmitter: samples .* cos(carrier_phase). The propagation channel
    % adds Doppler phase to the underlying carrier, not cos*cos mixing.
    carrier_phase=2*pi*fc*t;
    channel_phase=2*pi*fd*t+pi*rate*t.^2+phase;
    clean=samples.*cos(carrier_phase+channel_phase);
end

function est=receiver_estimate(x,pn,fs,fc,cfg,keep_detail)
    x=x(:); pn=pn(:); n=numel(x); tr=(0:n-1)'/fs;
    dec=round(fs/cfg.outFs); assert(abs(fs/dec-cfg.outFs)<1e-6);
    % Known PN only: data bits are never supplied to the receiver.
    base=resample(2*x.*pn.*exp(-1j*2*pi*fc*tr),1,dec);
    h=fir1(128,cfg.bandwidthHz/(cfg.outFs/2)); base=conv(base,h,'same');
    valid=65:numel(base)-64; y=base(valid); tm=(valid(:)-1)/cfg.outFs;
    reference=(tm(1)+tm(end))/2; t=tm-reference; T=numel(y)/cfg.outFs;
    assert(2*max(abs(cfg.fRange))+max(abs(cfg.rateRange))*T<cfg.outFs/2, ...
        'Squared chirp aliasing: increase baseband sample rate.');
    assert(cfg.rateStep<=1/T^2*1.01,'Rate grid too coarse for observation time.');
    window=0.5-0.5*cos(2*pi*(0:numel(y)-1)'/(numel(y)-1));
    q=y.^2.*window; nfft=2^nextpow2(2*numel(q));
    faxis=(-nfft/2:nfft/2-1)'*cfg.outFs/nfft/2;
    select=find(faxis>=cfg.fRange(1)&faxis<=cfg.fRange(2));
    rates=unique([cfg.rateRange(1):cfg.rateStep:cfg.rateRange(2),cfg.rateRange(2)]);
    best=-Inf; f0=NaN; a0=NaN;
    if keep_detail; scores=zeros(numel(rates),numel(select)); end
    for r=1:numel(rates)
        p=abs(fftshift(fft(q.*exp(-1j*2*pi*rates(r)*t.^2),nfft))).^2;
        if keep_detail; scores(r,:)=p(select).'; end
        [v,k]=max(p(select));
        if v>best; best=v; f0=faxis(select(k)); a0=rates(r); end
    end
    fun=@(p) receiver_score(p,q,t,f0,a0,T,cfg);
    [p,~,flag]=fminsearch(fun,[0 0],optimset('Display','off', ...
        'TolX',1e-6,'TolFun',1e-7,'MaxIter',300));
    est.doppler_mid_hz=f0+p(1)/T; est.rate_hz_per_s=a0+p(2)/T^2;
    est.carrier_mid_hz=fc+est.doppler_mid_hz; est.reference_time_s=reference;
    est.optimizer_converged=flag>0; est.coarse_doppler_hz=f0; est.coarse_rate=a0;
    if keep_detail
        % Display only a local slice; full fixed-prior search was performed.
        zoom=abs(faxis(select)-f0)<=100;
        est.detail=struct('y',y,'time_s',tm,'f_axis',faxis(select(zoom)).', ...
            'rates',rates,'metric',scores(:,zoom));
    end
end

function v=receiver_score(p,q,t,f0,a0,T,cfg)
    f=f0+p(1)/T; a=a0+p(2)/T^2;
    if f<cfg.fRange(1)||f>cfg.fRange(2)||a<cfg.rateRange(1)||a>cfg.rateRange(2)
        v=Inf; return;
    end
    v=-abs(sum(q.*exp(-1j*(4*pi*f*t+2*pi*a*t.^2))))^2 / max(sum(abs(q))^2,eps);
end

function [lo,hi]=wilson_interval(success,n)
    z=1.96; p=success/n; den=1+z^2/n;
    c=(p+z^2/(2*n))/den;
    half=z*sqrt(p.*(1-p)/n+z^2/(4*n^2))/den;
    lo=max(0,c-half); hi=min(1,c+half);
end
