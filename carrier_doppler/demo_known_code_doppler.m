% Known-code linear Doppler experiment, -7 dB Eb/N0, real passband.
% Run in the directory containing estimate_known_code_doppler.m.
% Requires Signal Processing Toolbox. Exact PN synchronization is assumed.
% Random 31-chip code is a test fixture, NOT a replacement for original Gold code.
clear; clc; rng(20261014,'twister');
fs=24.8e6; Rb=1e5; Rc=3.1e6; T=.1; fIF=1e5; EbN0_dB=-7;
N=round(fs*T); t=((0:N-1)-(N-1)/2)/fs;
cfg.outFs=4e5; cfg.bandwidthHz=9e4;
cfg.fRange=[-5e4 5e4]; cfg.rateRange=[-1000 1000]; cfg.rateStep=100;
numMonte=100; errors=zeros(numMonte,2);
for k=1:numMonte
    fd=(2*rand-1)*45000; rate=(2*rand-1)*900; phase=2*pi*rand;
    code=2*randi([0 1],1,31)-1;
    pn=code(mod(floor((0:N-1)/(fs/Rc)),31)+1);
    bits=2*randi([0 1],1,round(Rb*T))-1;
    data=bits(floor((0:N-1)/(fs/Rb))+1);
    clean=pn.*data.*cos(2*pi*(fIF+fd)*t+pi*rate*t.^2+phase);
    variance=mean(clean.^2)*fs/(2*Rb*10^(EbN0_dB/10));
    x=clean+sqrt(variance)*randn(1,N);
    est=estimate_known_code_doppler(x,pn,fs,fIF,cfg);
    % Function references its actual retained observation midpoint.
    expected_fd=fd+rate*(est.reference_time_s-(N-1)/(2*fs));
    errors(k,:)=[est.doppler_mid_hz-expected_fd,est.rate_hz_per_s-rate];
    if mod(k,10)==0; fprintf('trial %d/%d\n',k,numMonte); end
end
fTol=20; rateTol=200;
ok=abs(errors(:,1))<=fTol & abs(errors(:,2))<=rateTol;
p=mean(ok); z=1.96; n=numMonte;
c=(p+z^2/(2*n))/(1+z^2/n); half=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n);
fprintf('Joint success %.2f%%; Wilson 95%% [%.2f%% %.2f%%]\n', ...
    100*p,100*(c-half),100*(c+half));
fprintf('Frequency RMSE %.3f Hz, rate RMSE %.3f Hz/s\n', ...
    sqrt(mean(errors(:,1).^2)),sqrt(mean(errors(:,2).^2)));
% Original signal integration:
% pn must be the sampled Gold spreading sequence WITHOUT data bits.
% Provide the original mseq.m to reproduce the original sequence convention.
% Extend original t1 and regenerate all dependent arrays consistently.
% Replace constant-carrier phase with a linear-Doppler phase when testing.
% The estimator never receives true fd, rate, phase or data bits.
