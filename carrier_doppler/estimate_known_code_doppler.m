function est=estimate_known_code_doppler(x,pn,fs,fIF,cfg)
% Real DSSS-BPSK: synchronized known PN wipeoff, baseband conversion,
% squared unknown-data wipeoff, joint residual frequency/chirp-rate search.
% pn contains ONLY the sampled spreading code, never unknown data bits.
% Requires Signal Processing Toolbox. Linear Doppler and exact PN timing.
% cfg fields: outFs, fRange, rateRange, rateStep (Hz, Hz/s), bandwidthHz.
% fIF must be independently known; no true received-frequency input.
    x=double(x(:)); pn=double(pn(:));
    assert(isreal(x) && numel(x)==numel(pn));
    assert(all(abs(pn)==1),'PN samples must be +/-1.');
    D=round(fs/cfg.outFs); assert(abs(fs/D-cfg.outFs)<1e-6);
    N=numel(x); tr=(0:N-1)'/fs;
    base=resample(2*x.*pn.*exp(-1j*2*pi*fIF*tr),1,D);
    % Include anti-alias filtering before decimation via resample.
    h=fir1(128,cfg.bandwidthHz/(cfg.outFs/2));
    base=conv(base,h,'same');
    valid=65:numel(base)-64;
    assert(numel(valid)>100,'Insufficient observation after edge trimming.');
    y=base(valid); tm=(valid(:)-1)/cfg.outFs;
    centerTime=(tm(1)+tm(end))/2; t=tm-centerTime;
    T=numel(y)/cfg.outFs;
    assert(2*max(abs(cfg.fRange))+max(abs(cfg.rateRange))*T<cfg.outFs/2, ...
        'Squared chirp aliases: increase outFs or reduce search range.');
    window=0.5-0.5*cos(2*pi*(0:numel(y)-1)'/(numel(y)-1));
    q=y.^2.*window; nfft=2^nextpow2(2*numel(q));
    axis=(-nfft/2:nfft/2-1)'*cfg.outFs/nfft/2;
    idx=find(axis>=cfg.fRange(1) & axis<=cfg.fRange(2));
    rates=unique([cfg.rateRange(1):cfg.rateStep:cfg.rateRange(2),cfg.rateRange(2)]);
    best=-Inf; f0=NaN; a0=NaN;
    for a=rates
        power=abs(fftshift(fft(q.*exp(-1j*2*pi*a*t.^2),nfft))).^2;
        [v,j]=max(power(idx));
        if v>best; best=v; f0=axis(idx(j)); a0=a; end
    end
    objective=@(p) score(p,q,t,f0,a0,T,cfg);
    [p,~,flag]=fminsearch(objective,[0 0],optimset('Display','off', ...
        'TolX',1e-6,'TolFun',1e-7,'MaxIter',300));
    est.doppler_mid_hz=f0+p(1)/T;
    est.rate_hz_per_s=a0+p(2)/T^2;
    est.carrier_mid_hz=fIF+est.doppler_mid_hz;
    est.reference_time_s=centerTime;
    est.carrier_start_hz=est.carrier_mid_hz-est.rate_hz_per_s*centerTime;
    est.optimizer_converged=flag>0;
end
function v=score(p,q,t,f0,a0,T,cfg)
    f=f0+p(1)/T; a=a0+p(2)/T^2;
    if f<cfg.fRange(1)||f>cfg.fRange(2)||a<cfg.rateRange(1)||a>cfg.rateRange(2)
        v=Inf; return;
    end
    v=-abs(sum(q.*exp(-1j*(4*pi*f*t+2*pi*a*t.^2))))^2;
end
