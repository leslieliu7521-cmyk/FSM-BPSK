%% fixed700第一层f0重复性测试：默认50条独立记录
% 基于fixed700单记录版本另存，原文件保留；MATLAB R2018b+及Signal Processing Toolbox。
% 固定信息比特Eb/N0=-7dB、700Hz/s、观测中点100kHz、100ms；原PN同步假设不变。
% 每次独立生成、独立搜索后统计，不平均多次搜索指标。不作全先验性能外推。
% run('FSM_BPSK_fixed700_MC.m')
% 先检查流程：mc_options=struct('numTrials',3); run('FSM_BPSK_fixed700_MC.m');
clearvars -except mc_options; clc;
if ~exist('mc_options','var'); mc_options=struct(); end
opts=struct('numTrials',50,'baseSeed',20261015,'makeSummaryPlot',true, ...
    'singleRecordPlots',false,'outputRoot',pwd,'overview',false);
names=fieldnames(mc_options);
for ii=1:numel(names)
    assert(any(strcmp(names{ii},{'numTrials','baseSeed','makeSummaryPlot','singleRecordPlots','outputRoot'})), ...
        '未知选项：%s',names{ii});
    opts.(names{ii})=mc_options.(names{ii});
end
assert(isscalar(opts.numTrials) && isfinite(opts.numTrials) && opts.numTrials>=1 && opts.numTrials==floor(opts.numTrials));
assert(isscalar(opts.baseSeed) && isfinite(opts.baseSeed) && opts.baseSeed>=0 && opts.baseSeed==floor(opts.baseSeed));
assert(opts.baseSeed+opts.numTrials-1<=2^32-1,'种子超出mt19937ar支持范围。');
assert(exist('resample','file')~=0 && exist('freqz','file')~=0,'需要Signal Processing Toolbox。');
sourceText=fileread([mfilename('fullpath') '.m']);
%% 一、固定科学配置，保持原第一层先验及评分
cfg.Rb=1e5; cfg.Rc=31e5; cfg.fs=24.8e6; cfg.fc=1e5;
cfg.T=.1; cfg.longT=.1; cfg.scfFs=8e5; cfg.L=128; cfg.hop=64;
cfg.smoothingBins=17; cfg.fRange=[50e3 150e3];
cfg.rateRange=[-1000 1000]; cfg.rateStep=100;
cfg.fTolerance=20; cfg.rateTolerance=200; cfg.fdTest=45000; cfg.rateTest=900;
cfg.useKnownPN=true; cfg.rxPhaseChip=0; cfg.rxRatePPM=0; cfg.jerk=0;
cfg.filterN=10; cfg.filterBeta=5; cfg.branchPeaks=3; cfg.confirmCandidates=10;
cfg.multiStarts=0; cfg.dedupFrequencyFactor=.5; cfg.dedupRateFactor=2;
cfg.maxRefineIter=200;

%% 二、每条记录独立种子；预分配全部试验，异常也保留行
outputDir=unique_directory(fullfile(opts.outputRoot,['FSM_fixed700_MC_' datestr(now,'yyyymmdd_HHMMSS')]));
mkdir(outputDir); code=make_gold(); N=opts.numTrials;
seeds=opts.baseSeed+(0:N-1)';
trials=table((1:N)',seeds,false(N,1),false(N,1),false(N,1), ...
    NaN(N,1),NaN(N,1),NaN(N,1),NaN(N,1),NaN(N,1),NaN(N,1),NaN(N,1), ...
    NaN(N,1),NaN(N,1),zeros(N,1),repmat({''},N,1), ...
    'VariableNames',{'TrialID','Seed','Completed','Valid','Success20Hz', ...
    'ReferenceTime_s','TrueF0_Hz','EstimatedF0_Hz','F0Error_Hz','TrueRate_HzPerS', ...
    'EstimatedRate_HzPerS','RateError_HzPerS','ActualRawSampleSNR_dB','ComputeTime_s', ...
    'ContinuousRefinementCalls','ErrorMessage'});
trials.RuntimeException=false(N,1);
% 固定真值列即使该次运行异常也保留；仅供事后统计，不传入接收机。
trials.TrueF0_Hz(:)=cfg.fc; trials.TrueRate_HzPerS(:)=700;
trials.ReferenceTime_s(:)=(round(cfg.T*cfg.scfFs)-1)/(2*cfg.scfFs);
state=struct('version','fixed700_MC_1','cfg',cfg,'options',opts,'sourceText',sourceText, ...
    'seeds',seeds,'trials',trials,'complete',false,'representative',[], ...
    'note','固定真值重复性测试；700在变化率网格上，100kHz在第一层分支中心频率网格上。');
wholeTimer=tic;
for trial=1:N
    trialTimer=tic;
    try
        % 可选诊断只在预先指定的第1条记录计算，不按真实误差挑选。
        trialOpts=opts; trialOpts.singleRecordPlots=opts.singleRecordPlots && trial==1;
        q=coordinate_case(cfg,code,seeds(trial),trialOpts);
        trials.ReferenceTime_s(trial)=q.referenceTime_s;
        trials.TrueF0_Hz(trial)=q.F_ref; trials.EstimatedF0_Hz(trial)=q.F_coarse;
        trials.TrueRate_HzPerS(trial)=q.trueRate;
        trials.EstimatedRate_HzPerS(trial)=q.estimatedRate;
        trials.ActualRawSampleSNR_dB(trial)=q.actualRawSampleSNR_dB;
        trials.ContinuousRefinementCalls(trial)=q.search.continuousRefinementCalls;
        valid=q.search.coarseValid && isfinite(q.F_coarse) && isfinite(q.estimatedRate);
        trials.Valid(trial)=valid;
        if valid
            trials.F0Error_Hz(trial)=q.F_coarse-q.F_ref;
            trials.RateError_HzPerS(trial)=q.estimatedRate-q.trueRate;
        else
            trials.ErrorMessage{trial}='无效/非有限估计：误差按NaN保留，计为失败。';
        end
        trials.Success20Hz(trial)=valid && abs(trials.F0Error_Hz(trial))<=20;
        if trial==1 && opts.singleRecordPlots && valid; state.representative=q; end
    catch ME
        % 运行异常独立列出，并占用分母；不补抽样或把异常删除为未试验。
        trials.RuntimeException(trial)=true;
        trials.ErrorMessage{trial}=[ME.identifier ': ' ME.message];
        warning('FSM:TrialException','Trial %d异常，计为无效/失败：%s',trial,ME.message);
    end
    trials.ComputeTime_s(trial)=toc(trialTimer); trials.Completed(trial)=true;
    state.trials=trials; state.completedTrials=trial;
    % 每条完成后保存必要标量；没有完整波形。中断文件不冒充完整统计。
    state.elapsedLoop_s=toc(wholeTimer);
    save(fullfile(outputDir,'checkpoint.mat'),'state','-v7.3');
    fprintf('fixed700 %d/%d：f0误差=%g Hz，有效=%d，用时=%.3f s\n', ...
        trial,N,trials.F0Error_Hz(trial),trials.Valid(trial),trials.ComputeTime_s(trial));
    clear q;
end
%% 三、全分母统计：无效不删除，不只统计成功样本
state.totalElapsed_s=toc(wholeTimer); state.complete=all(trials.Completed);
summary=mc_summary(trials,state.totalElapsed_s,cfg);
state.summary=summary; state.trials=trials;
% 最终MAT先于CSV/图像导出；CSV失败也不会丢失数值。
save(fullfile(outputDir,'results.mat'),'state','-v7.3');
writetable(trials,fullfile(outputDir,'trials.csv'));
writetable(summary,fullfile(outputDir,'summary.csv'));
disp(summary);
fprintf(['目标信息比特Eb/N0=-7 dB；目标原始采样SNR=%.9g dB。\n' ...
    '无效误差按NaN存储：只要任一误差非有限，完整样本RMSE=Inf，偏差/样本标准差=NaN。\n' ...
    '失败数包含无效及有效但超20Hz门限的记录。没有仅成功样本统计。\n' ...
    '700 Hz/s在原变化率网格上；100kHz在原第一层分支中心频率网格上。\n' ...
    '固定真值重复性测试，不代表覆盖所有频率/变化率，也不外推总体准确率。\n'], ...
    -7+10*log10(2*cfg.Rb/cfg.fs));
%% 四、只画一幅默认汇总图；诊断图必须明确启用
plotReport=struct('summary',[],'singleRecord',{{}},'errors',{{}});
if opts.makeSummaryPlot
    try
        h=figure('Color','w','Name','fixed700第一层：逐次f0误差');
        ax=axes('Parent',h,'Position',[.09 .16 .58 .72]);
        plot(ax,trials.TrialID,trials.F0Error_Hz,'b.-'); hold(ax,'on');
        plot(ax,[.5 N+.5],[20 20],'r--'); plot(ax,[.5 N+.5],[-20 -20],'r--');
        xlim(ax,[.5 N+.5]); grid(ax,'on');
        xlabel(ax,'独立试验编号'); ylabel(ax,'中点f0估计误差 (Hz)');
        title(ax,'fixed700第一层；Eb/N0=-7dB，100ms');
        invalid=find(~trials.Valid);
        if ~isempty(invalid)
            plot(ax,invalid,zeros(size(invalid)),'rx','MarkerSize',9);
            text(ax,.02,.02,'红叉仅标识无效试验编号，纵坐标0不表示误差为0。', ...
                'Units','normalized','Interpreter','none');
        end
        msg=sprintf(['N=%d\n成功=%d (|e|<=20Hz)\n成功率=%.2f%%\nWilson双侧95%%\n[%.2f%%, %.2f%%]\nRMSE=%.6g Hz\n偏差=%.6g Hz\n样本标准差=%.6g Hz\n无效=%d，失败=%d\n变化率RMSE=%.6g Hz/s\n总用时=%.3f s\n仅固定真值重复性测试'], ...
            N,summary.SuccessCount,100*summary.SuccessRate,100*summary.WilsonLo95,100*summary.WilsonHi95, ...
            summary.F0RMSE_Hz,summary.F0Bias_Hz,summary.F0SampleStd_Hz, ...
            summary.InvalidCount,summary.FailureCount,summary.RateRMSE_HzPerS,summary.TotalElapsed_s);
        annotation(h,'textbox',[.70 .12 .29 .8],'String',msg,'Interpreter','none','EdgeColor','none');
        plotReport.summary=export_figures({h},outputDir,false);
    catch ME; plotReport.errors{end+1}=ME.message; warning('FSM:SummaryPlot','数值已保存：%s',ME.message); end
end
if opts.singleRecordPlots && ~isempty(state.representative)
    for kind=1:4
        try
            h=coordinate_figure(state.representative,kind);
            folder=fullfile(outputDir,sprintf('optional_trial1_plot%d',kind)); mkdir(folder);
            plotReport.singleRecord{end+1}=export_figures({h},folder,false);
        catch ME; plotReport.errors{end+1}=ME.message; end
    end
end
save(fullfile(outputDir,'plot_report.mat'),'plotReport');
fprintf('结果目录：%s\n',outputDir);

%% 五、Wilson双侧95%%区间及完整样本误差规则
function summary=mc_summary(T,totalSeconds,c)
    n=height(T); k=sum(T.Success20Hz); p=k/n; z=1.959963984540054;
    den=1+z^2/n; mid=(p+z^2/(2*n))/den;
    half=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/den;
    [rmse,bias,sd]=complete_error_stats(T.F0Error_Hz);
    [rateRMSE,rateBias,rateSD]=complete_error_stats(T.RateError_HzPerS);
    summary=table(n,k,p,max(0,mid-half),min(1,mid+half),rmse,bias,sd, ...
        sum(~T.Valid),n-k,sum(T.RuntimeException), ...
        rateRMSE,rateBias,rateSD,totalSeconds,totalSeconds/n,sum(T.ComputeTime_s), ...
        mean(T.ComputeTime_s),-7,-7+10*log10(2*c.Rb/c.fs), ...
        sum(T.ContinuousRefinementCalls), ...
        'VariableNames',{'NumTrials','SuccessCount','SuccessRate','WilsonLo95','WilsonHi95', ...
        'F0RMSE_Hz','F0Bias_Hz','F0SampleStd_Hz','InvalidCount','FailureCount','RuntimeExceptionCount', ...
        'RateRMSE_HzPerS','RateBias_HzPerS','RateSampleStd_HzPerS','TotalElapsed_s', ...
        'MeanElapsedPerTrial_s','TotalCompute_s','MeanComputePerTrial_s','TargetEbN0_dB', ...
        'TargetRawSampleSNR_dB','ContinuousRefinementCalls'});
end
function [rmse,bias,sd]=complete_error_stats(e)
    if any(~isfinite(e)); rmse=Inf; bias=NaN; sd=NaN; return; end
    rmse=sqrt(mean(e.^2)); bias=mean(e);
    if numel(e)>1; sd=std(e,0); else; sd=NaN; end
end

%% 六、保留fixed700信号与第一层函数；仅在单记录诊断关闭时跳过绘图数组计算
function out=coordinate_case(c,code,seed,opts)
    b=make_base(c,code,seed);
    eb=-7; fd=0; a=700; label='fixed700：仅第一层，单条含噪记录';
    % 先按原重采样长度/中心截取规则确定观测段中点。原始与重采样
    % 采样网格中点相差约0.6us；将相位参考设在接收观测段中点，
    % 避免仅在原始记录中点fd=0而接收评价中点有微小偏差。
    [up,down]=rat(c.scfFs/c.fs,1e-12);
    nExpected=ceil(numel(b.rawTime)*up/down); n=round(c.T*c.scfFs);
    firstExpected=floor((nExpected-n)/2)+1;
    observationReference=(firstExpected-1)/c.scfFs+(n-1)/(2*c.scfFs);
    b.reference=observationReference; b.t=b.rawTime-b.reference;
    phase=2*pi*fd*b.t+pi*a*b.t.^2+pi/3*c.jerk*b.t.^3+b.phase;
    clean=b.data.*b.pn.*cos(2*pi*c.fc*b.t+phase);
    % 实AWGN双边谱密度N0/2；这里-7是信息比特Eb/N0而非采样SNR。
    signalPower=mean(clean.^2);
    noiseVariance=signalPower*c.fs/(2*c.Rb*10^(-7/10));
    noise=sqrt(noiseVariance)*b.noise;
    rx=clean+noise;
    targetSampleSNR_dB=eb+10*log10(2*c.Rb/c.fs);
    actualNoisePower=mean(noise.^2);
    actualSampleSNR_dB=10*log10(signalPower/actualNoisePower);
    if c.useKnownPN
        Ns=c.fs/c.Rc;
        ci=floor((0:numel(rx)-1)*(1+c.rxRatePPM*1e-6)/Ns+c.rxPhaseChip);
        pnRx=code(mod(ci,31)+1); clean=clean.*pnRx; rx=rx.*pnRx;
        stage='解扩后、重采样及中心截取后';
    else; stage='解扩前（未解扩）、重采样及中心截取后'; end
    [up,down]=rat(c.scfFs/c.fs,1e-12);
    yCleanLong=resample(clean(:),up,down,c.filterN,c.filterBeta);
    yLong=resample(rx(:),up,down,c.filterN,c.filterBeta);
    n=round(c.T*c.scfFs); first=floor((numel(yLong)-n)/2)+1;
    assert(first>=1 && first+n-1<=numel(yLong));
    assert(numel(yLong)==nExpected && first==firstExpected,'实际重采样长度与中点预计算不一致。');
    y=yLong(first:first+n-1); yClean=yCleanLong(first:first+n-1);
    r=receiver_config(c);
    e=fsm_receiver(y,r);
    % 无效估计作为试验结果保留，不删除试验。
    if ~isfield(e,'continuousRefinementCalls'); e.continuousRefinementCalls=0; end
    reference=(first-1)/c.scfFs+e.reference; delta=reference-b.reference;
    F_ref=c.fc+fd+a*delta+.5*c.jerk*delta^2; trueRate=a+c.jerk*delta;
    assert(abs(delta)<1e-12 && abs(F_ref-c.fc)<1e-8,'f0参考时刻不一致。');
    assert(isempty(e.refinements) && e.continuousRefinementCalls==0,'不允许执行第二层。');
    % 默认Monte Carlo只返回必要标量，不做PSD/局部扫描、不保存波形。
    out=struct('seed',seed,'F_ref',F_ref,'F_coarse',e.fCoarse, ...
        'trueRate',trueRate,'estimatedRate',e.aCoarse,'referenceTime_s',reference, ...
        'actualRawSampleSNR_dB',actualSampleSNR_dB,'search',e);
    if ~opts.singleRecordPlots || ~e.coarseValid; return; end
    [frames,tc,~,w]=prepare_frames(y,r);
    band=-(c.smoothingBins-1)/2:(c.smoothingBins-1)/2;
    % 局部网格只由接收机结果和固定先验决定，绝不由真值决定。
    f_candidate=unique([linspace(max(c.fRange(1),e.fFine-50),min(c.fRange(2),e.fFine+50),41),e.fFine]);
    local=coordinate_grid(frames,tc,w,r,f_candidate,e.aFine);
    raw=fsm_spectrum(frames,tc,e.fFine,e.aFine,r,w);
    smooth=conv(raw,ones(c.smoothingBins,1)/c.smoothingBins,'same');
    exact=fsm_score(frames,tc,e.fFine,e.aFine,r,w,band);
    assert(abs(abs(smooth(c.L/2+1))^2-exact)<=1e-10*max(exact,realmin), ...
        '切片中心评分与原fsm_score不一致。');
    overview=[];
    if opts.overview
        % 仅51个固定先验载波点，稀疏概览不用于估计，也不保证命中窄峰。
        overview=coordinate_grid(frames,tc,w,r,linspace(c.fRange(1),c.fRange(2),51),e.aFine);
    end
    Pclean=pwelch(yClean,w,c.L-c.hop,c.L,c.scfFs,'twosided');
    Pnoisy=pwelch(y,w,c.L-c.hop,c.L,c.scfFs,'twosided');
    out=struct('label',label,'stage',stage,'seed',seed,'ebN0_dB',eb,'cfg',c, ...
        'fc_nominal',c.fc,'F_ref',F_ref,'F_coarse',e.fCoarse,'F_hat',e.fFine, ...
        'trueRate',trueRate,'coarseRate',e.aCoarse,'estimatedRate',e.aFine, ...
        'referenceTime_s',reference,'transmitReferenceTime_s',b.reference, ...
        'search',e,'local',local,'overview',overview,'rawSlice',raw, ...
        'smoothSlice',smooth,'exactScore',exact,'frameCount',size(frames,2), ...
        'sampleCount',n, ...
        'PSDclean',fftshift(Pclean),'PSDnoisy',fftshift(Pnoisy));
    % pwelch双边输出为[0,fs)，fftshift后使用明确的[-fs/2,fs/2)坐标。
    out.PSD_f_Hz=(-c.L/2:c.L/2-1)'*c.scfFs/c.L;
    out.targetEbN0_dB=eb;
    out.targetRawSampleSNR_dB=targetSampleSNR_dB;
    out.actualRawSampleSNR_dB=actualSampleSNR_dB;
    out.rawSignalPower=signalPower; out.rawNoiseVariance=noiseVariance;
    out.actualRawNoisePower=actualNoisePower;
    out.f0Reference='观测段中点，从原始记录起点计时';
    out.f0Error_Hz=e.fCoarse-F_ref;
    out.rateError_HzPerS=e.aCoarse-trueRate;
    out.continuousRefinementCalls=e.continuousRefinementCalls;
    out.zeroCase=false;
    % 特例只评价测得峰坐标，不对错误结果强行移动坐标或保证通过。
    out.zeroCoordinatePass=false && abs(local.surfacePeak_f_spec_Hz)<=c.scfFs/c.L ...
        && abs(local.scorePeak_alpha_Hz-2*c.fc)<=2/max(c.T,eps);
end
function g=coordinate_grid(frames,tc,w,c,f_candidate,a_hat)
    f_spec=(-c.L/2:c.L/2-1)*c.scfFs/c.L;
    S=complex(zeros(numel(f_candidate),c.L)); raw=S;
    for ii=1:numel(f_candidate)
        z=fsm_spectrum(frames,tc,f_candidate(ii),a_hat,c,w);
        raw(ii,:)=z.';
        S(ii,:)=conv(z,ones(c.smoothingBins,1)/c.smoothingBins,'same').';
    end
    center=c.L/2+1; J=abs(S(:,center)).^2;
    [peak,ind]=max(abs(S(:))); [ia,isp]=ind2sub(size(S),ind);
    [~,ij]=max(J);
    g=struct('f_spec_Hz',f_spec,'alpha_Hz',2*f_candidate,'f_candidate_Hz',f_candidate, ...
        'rawComplex',raw,'smoothComplex',S,'scoreAtZero',J, ...
        'surfacePeakAmplitude',peak,'surfacePeak_f_spec_Hz',f_spec(isp), ...
        'surfacePeak_alpha_Hz',2*f_candidate(ia),'scorePeak_alpha_Hz',2*f_candidate(ij), ...
        'alphaGridSteps_Hz',diff(2*f_candidate),'estimatedRate',a_hat);
end

%% 六、只读取已经完成的基础含噪记录；不把方法对照重复计成独立试验
function text=coordinate_report(d)
    text=sprintf('循环谱坐标诊断；每场景1条记录，非Monte Carlo统计。\n%s\n新生成记录=%d，新接收机运行=%d\n', ...
        d.reuseNote,d.actualNewRecords,d.actualNewReceiverRuns);
    for k=1:numel(d.cases)
        q=d.cases{k}; g=q.local; c=q.cfg;
        text=[text sprintf(['目标信息比特Eb/N0=%.6g dB；目标原始采样SNR=%.9g dB；实际原始采样SNR=%.9g dB\n' ...
            '原始解扩前信号功率=%.12g；设置噪声方差=%.12g；实测噪声功率=%.12g\n' ...
            '观测段中点f0误差=%.12g Hz；变化率误差=%.12g Hz/s；连续精搜索调用=%d\n' ...
            '最终输出等于第一层；旧fFine/aFine字段仅用于绘图兼容。\n'], ...
            q.targetEbN0_dB,q.targetRawSampleSNR_dB,q.actualRawSampleSNR_dB, ...
            q.rawSignalPower,q.rawNoiseVariance,q.actualRawNoisePower, ...
            q.f0Error_Hz,q.rateError_HzPerS,q.continuousRefinementCalls)]; %#ok<AGROW>
        text=[text sprintf(['\n场景%d：%s，seed=%d，Eb/N0=%g dB\n' ...
            'fc_nominal=%.12g Hz；F_ref=%.12g Hz；第一层=%.12g Hz；最终=%.12g Hz\n' ...
            '2fc_nominal=%.12g Hz；2F_ref=%.12g Hz；2F_coarse=%.12g Hz；2F_hat=%.12g Hz\n' ...
            '真实/第一层/输出变化率=%.12g / %.12g / %.12g Hz/s\n' ...
            '接收参考时刻=%.12g s（从原始长记录起点计）；发射参考=%.12g s\n' ...
            '局部平滑幅度曲面最大点：(f_spec,alpha)=(%.12g,%.12g) Hz\n' ...
            'f_spec=0评分切片网格最大点alpha=%.12g Hz，F_candidate=%.12g Hz\n' ...
            '谱频率间隔=%.12g Hz；alpha步长范围=[%.12g,%.12g] Hz（含估计点可非均匀）\n' ...
            '观测=%.12g s，N=%d，FFT=%d，hop=%d，帧数=%d\n' ...
            '频率平滑=%d点，名义宽度=%.12g Hz，首末中心跨度=%.12g Hz\n' ...
            '最终原指标J=%.12g，固定a_hat局部网格最大J=%.12g；实际试验数=1\n'], ...
            k,q.label,q.seed,q.ebN0_dB,q.fc_nominal,q.F_ref,q.F_coarse,q.F_hat, ...
            2*q.fc_nominal,2*q.F_ref,2*q.F_coarse,2*q.F_hat, ...
            q.trueRate,q.coarseRate,q.estimatedRate,q.referenceTime_s,q.transmitReferenceTime_s, ...
            g.surfacePeak_f_spec_Hz,g.surfacePeak_alpha_Hz,g.scorePeak_alpha_Hz,g.scorePeak_alpha_Hz/2, ...
            c.scfFs/c.L,min(g.alphaGridSteps_Hz),max(g.alphaGridSteps_Hz), ...
            c.T,q.sampleCount,c.L,c.hop,q.frameCount,c.smoothingBins, ...
            c.smoothingBins*c.scfFs/c.L,(c.smoothingBins-1)*c.scfFs/c.L,q.exactScore,max(g.scoreAtZero))]; %#ok<AGROW>
        if q.zeroCase
            text=[text sprintf('零频偏案例坐标容差检查=%d（仅此单例，非统计；补偿值仍来自接收机）。\n',q.zeroCoordinatePass)]; %#ok<AGROW>
        end
        text=[text sprintf(['曲面峰对所有f_spec取最大；估计指标只在f_spec=0取复数平滑后模平方。\n' ...
            '两者不同不能用移动坐标掩盖；局部扫描固定a_hat，原估计还搜索变化率。\n' ...
            '平滑卷积边缘使用same零填充，边缘值仅供显示；中心评分不受边缘填充影响。\n'])]; %#ok<AGROW>
    end
    text=[text sprintf(['\nPSD不保证离散载波尖线，PSD峰不作为载波估计。局部曲面不是全局捕获证据。\n' ...
        '本文件仅为fixed700单记录诊断，不代表完成性能统计实验。\n'])];
end

%% 八、分图绘制：每张图都有独立异常边界，真值只用于辅助标线
function h=coordinate_figure(q,kind)
    g=q.local; c=q.cfg;
    h=figure('Name',[q.label ' / ' num2str(kind)],'Color','w');
    caption=sprintf('Eb/N0=%.3g dB；F_ref/F_hat=%.9g/%.9g Hz；a_true/a_hat=%.6g/%.6g Hz/s；t_ref=%.9g s',q.ebN0_dB,q.F_ref,q.F_hat,q.trueRate,q.estimatedRate,q.referenceTime_s);
    switch kind
        case 1
            plot(q.PSD_f_Hz/1e3,10*log10(max(q.PSDclean,realmin)),'b'); hold on;
            plot(q.PSD_f_Hz/1e3,10*log10(max(q.PSDnoisy,realmin)),'Color',[.7 .3 .1]);
            legend('同记录干净信号','同记录含噪信号','Location','best');
            for sign=[-1 1]
                coordinate_line(sign*q.fc_nominal/1e3,'标称','k--');
                coordinate_line(sign*q.F_ref/1e3,'真值参考（仅评价）','g:');
            end
            xlabel('双边频率 (kHz)'); ylabel('PSD (信号单位^2/Hz，dB)'); grid on;
            title({[q.label '：' q.stage],caption,'BPSK不保证离散载波尖线；PSD峰不是载波估计'},'Interpreter','none');
        case {2,5}
            if kind==5; g=q.overview; end
            surf(g.f_spec_Hz/1e3,g.alpha_Hz/1e3,abs(g.smoothComplex),'EdgeColor','none'); hold on;
            center=c.L/2+1;
            z=abs(q.smoothSlice(center));
            plot3(0,2*q.F_hat/1e3,z,'ro','MarkerFaceColor','r');
            % 真值标记落在底平面，避免把插值值伪称实测真值谱值。
            plot3(0,2*q.F_ref/1e3,0,'gx','MarkerSize',10,'LineWidth',2);
            plot3(g.surfacePeak_f_spec_Hz/1e3,g.surfacePeak_alpha_Hz/1e3,g.surfacePeakAmplitude,'kd');
            xlabel('谱频率 f_spec (kHz)'); ylabel('循环频率 alpha (kHz)');
            zlabel('|复数频率平滑后的循环谱| (信号单位^2/Hz)'); view(45,30); colorbar;
            if kind==2; name='以接收机估计为中心的局部展示；不是全局捕获证据';
            else; name='固定先验稀疏概览：仅显示，不参与估计，可能漏过窄峰'; end
            title({name,caption,'红圈：估计；绿叉：真值坐标参考（底平面）；黑菱形：实测曲面峰'},'Interpreter','none');
            % 局部轴固定为实际计算范围；范围外真值用文字，不扩展数据。
            ylim([min(g.alpha_Hz) max(g.alpha_Hz)]/1e3);
            coordinate_outside_note(2*q.F_ref,g.alpha_Hz,'真值2F_ref');
        case 3
            plot(g.f_spec_Hz/1e3,abs(q.rawSlice),'b.-'); hold on;
            plot(g.f_spec_Hz/1e3,abs(q.smoothSlice),'r.-');
            legend('平滑前幅度','复数平滑后幅度','Location','best');
            coordinate_line(0,'f_spec=0','k:');
            xlabel('谱频率 f_spec (kHz)'); ylabel('循环谱幅度 (信号单位^2/Hz)'); grid on;
            title({sprintf('固定 alpha=2F_hat=%.9g kHz',2*q.F_hat/1e3),caption, ...
                '正负载频相关对的中心约为0；这里不是alpha=0'},'Interpreter','none');
        case 4
            subplot(2,1,1); plot(g.alpha_Hz/1e3,g.scoreAtZero,'b.-'); hold on;
            coordinate_marks(q,g.alpha_Hz,2); grid on;
            xlabel('循环频率 alpha (kHz)'); ylabel('J=|复数平滑S(0)|^2 (单位^4/Hz^2)');
            title({'固定 f_spec=0：与fsm_score同顺序、同评价位置',caption},'Interpreter','none');
            subplot(2,1,2); plot(g.f_candidate_Hz/1e3,g.scoreAtZero,'b.-'); hold on;
            coordinate_marks(q,g.f_candidate_Hz,1); grid on;
            xlabel('接收载波候选 F_candidate=alpha/2 (kHz)'); ylabel('同一评分 J (单位^4/Hz^2)');
            title('只重标横轴，不移动数据；局部扫描固定接收机估计变化率');
    end
end
function coordinate_marks(q,axisHz,factor)
    xlim([min(axisHz) max(axisHz)]/1e3);
    values=factor*[q.fc_nominal q.F_ref q.F_hat];
    names={'标称','真值参考（仅评价）','第一层输出'}; styles={'k--','g:','r-.'};
    for k=1:3
        if values(k)>=min(axisHz) && values(k)<=max(axisHz)
            coordinate_line(values(k)/1e3,names{k},styles{k});
        else
            text(.02,.95-.12*(k-1),sprintf('%s %.9g kHz 在局部范围外',names{k},values(k)/1e3), ...
                'Units','normalized','Interpreter','none');
        end
    end
end
function coordinate_outside_note(value,axisHz,name)
    if value<min(axisHz) || value>max(axisHz)
        text(.02,.95,sprintf('%s=%.9g kHz 在局部范围外',name,value/1e3), ...
            'Units','normalized','Interpreter','none');
    end
end
function coordinate_line(value,label,style)
    % 不依赖xline；普通plot支持R2018b，固定y范围防止标线改变纵轴。
    yl=ylim; plot([value value],yl,style,'HandleVisibility','off'); ylim(yl);
    text(value,yl(2),label,'VerticalAlignment','top','Rotation',90,'Interpreter','none');
end

%% 九、保留原接收计算函数及必要辅助函数（不包含全实验总控）
function b=make_base(c,code,seed)
    s=RandStream('mt19937ar','Seed',seed);
    nb=round(c.Rb*c.longT); Ns=c.fs/c.Rc; N=nb*31*Ns;
    b.bits=2*randi(s,[0 1],1,nb)-1;
    b.phase=2*pi*rand(s); b.fd=(2*rand(s)-1)*c.fdTest; b.rate=(2*rand(s)-1)*c.rateTest;
    b.noise=randn(s,1,N); % 所有配对条件共用这一条标准噪声，绝不另抽短噪声
    b.pn=repelem(repmat(code,1,nb),Ns); b.data=repelem(b.bits,31*Ns);
    b.rawTime=(0:N-1)/c.fs; b.reference=(N-1)/(2*c.fs); b.t=b.rawTime-b.reference;
    b.seed=seed; b.rngState=struct('Type',s.Type,'State',s.State,'Seed',seed);
end

function r=receiver_config(c)
    % 明确的接口白名单：任何真实发送参数都不在此结构中。
    names={'scfFs','L','hop','smoothingBins','fRange','rateRange','rateStep', ...
        'branchPeaks','confirmCandidates','multiStarts','dedupFrequencyFactor', ...
        'dedupRateFactor','maxRefineIter'};
    r=struct(); for i=1:numel(names); r.(names{i})=c.(names{i}); end
end

function e=fsm_receiver(y,c)
    [frames,tc,ts,w,n]=prepare_frames(y,c);
    e=empty_search(n,c); e.frameCount=size(frames,2);
    if isempty(frames) || any(~isfinite(frames(:))); return; end
    T=n/c.scfFs; dF=c.dedupFrequencyFactor/T; dA=c.dedupRateFactor/T^2;
    e.dedupScale=[dF dA];
    band=-(c.smoothingBins-1)/2:(c.smoothingBins-1)/2;
    X=fft(frames,[],1); df=c.scfFs/c.L;
    nt=2^nextpow2(2*numel(ts));
    deltaAlpha=(-nt/2:nt/2-1)*(c.scfFs/c.hop)/nt;
    rates=unique([c.rateRange(1):c.rateStep:c.rateRange(2),c.rateRange(2)]);
    generated=zeros(0,4);
    for k=floor(c.fRange(1)/df):ceil(c.fRange(2)/df)
        fb=k*df;
        cross=X(mod(band+k,c.L)+1,:).*conj(X(mod(band-k,c.L)+1,:));
        p=mean(cross,1).*exp(-1j*4*pi*fb*ts); % 先复数频率平均，保留相位
        bank=bsxfun(@times,p,exp(-1j*2*pi*rates(:)*ts.^2));
        power=abs(fftshift(fft(bank,nt,2),2)).^2;
        f=fb+deltaAlpha/2;
        valid=deltaAlpha>=-df & deltaAlpha<df & f>=c.fRange(1) & f<=c.fRange(2);
        branch=branch_candidates(power(:,valid),f(valid),rates,c.branchPeaks,dF,dA,k);
        generated=[generated;branch]; %#ok<AGROW>
    end
    e.generated=generated; % 列：f(Hz),a(Hz/s),近似分数,FFT分支号
    if isempty(generated); return; end
    [~,order]=sort(generated(:,3),'descend'); sorted=generated(order,:);
    keep=distinct_indices(sorted(:,1:2),dF,dA);
    distinct=sorted(keep,:);
    K=size(distinct,1); confirm=min(K,c.confirmCandidates);
    % 全部候选对照指全部“生成后去重”的不同候选；原始重复候选也保存以便审计。
    C=[distinct,NaN(K,1),(1:K)',NaN(K,1),zeros(K,1)];
    e.candidateColumns={'Frequency_Hz','Rate_HzPerS','ApproxScore','Branch', ...
        'ExactScore','ApproxRank','ExactRank','Confirmed'};
    for j=1:confirm
        C(j,5)=fsm_score(frames,tc,C(j,1),C(j,2),c,w,band); C(j,8)=1;
    end
    available=find(C(:,8)>0 & isfinite(C(:,5)) & C(:,5)>=0);
    if isempty(available); e.candidates=C; return; end
    [~,r]=sort(C(available,5),'descend'); ranking=available(r);
    C(ranking,7)=(1:numel(ranking))'; e.candidates=C; e.exactRanking=ranking;
    best=ranking(1); e.fCoarse=C(best,1); e.aCoarse=C(best,2);
    e.coarseScore=C(best,5); e.coarseValid=true;
    % 粗估计定义：准确FSM确认后的最佳离散候选，不含连续fminsearch结果。

    %% 十一、fixed700只用第一层：不创建目标函数、不调用fminsearch
    % fFine/aFine仅作旧绘图兼容字段，值严格等于第一层估计，绝非精搜索。
    e.fFine=e.fCoarse; e.aFine=e.aCoarse; e.finalScore=e.coarseScore;
    e.finalValid=e.coarseValid; e.refinements=zeros(0,8);
    e.continuousRefinementCalls=0; e.selectedStart=0; e.selectedExitFlag=NaN;
    e.fallback=false; e.fallbackReason='SecondLayerDisabled_FirstLayerOutput';
    ftol=min(dF,.01*diff(c.fRange)); atol=min(dA,.01*diff(c.rateRange));
    e.boundary=min(abs(e.fCoarse-c.fRange))<=ftol || min(abs(e.aCoarse-c.rateRange))<=atol;
end

function rows=branch_candidates(P,f,rates,count,dF,dA,branch)
    rows=zeros(0,4); if isempty(P); return; end
    P(~isfinite(P))=-Inf;
    for j=1:count
        [score,linear]=max(P(:)); if ~isfinite(score); break; end
        [a,b]=ind2sub(size(P),linear); fp=f(b); ap=rates(a);
        rows(end+1,:)=[fp ap score branch]; %#ok<AGROW>
        % 椭圆距离而非固定bin距离；T改变后依然采用同一循环相位容差尺度。
        distance=bsxfun(@plus,((rates(:)-ap)/dA).^2,((f-fp)/dF).^2);
        P(distance<1)=-Inf;
    end
end

function indices=distinct_indices(points,dF,dA)
    indices=zeros(0,1);
    for j=1:size(points,1)
        if isempty(indices); indices=j; continue; end
        distance=((points(indices,1)-points(j,1))/dF).^2+((points(indices,2)-points(j,2))/dA).^2;
        if all(distance>=1); indices(end+1,1)=j; end %#ok<AGROW>
    end
end

function yes=in_prior(point,c)
    yes=numel(point)==2 && all(isfinite(point)) && point(1)>=c.fRange(1) ...
        && point(1)<=c.fRange(2) && point(2)>=c.rateRange(1) && point(2)<=c.rateRange(2);
end

function e=empty_search(n,c)
    e=struct('fCoarse',NaN,'aCoarse',NaN,'coarseScore',NaN,'coarseValid',false, ...
        'fFine',NaN,'aFine',NaN,'finalScore',NaN,'finalValid',false,'reference',(n-1)/(2*c.scfFs), ...
        'fallback',true,'fallbackReason','NoValidCandidates','selectedStart',0,'selectedExitFlag',NaN, ...
        'boundary',false,'frameCount',0,'generated',zeros(0,4),'candidates',zeros(0,8), ...
        'refinements',zeros(0,8),'dedupScale',[NaN NaN],'exactRanking',[], ...
        'candidateColumns',{{}},'refinementColumns',{{}},'refineMessages',{{}});
    e.candidateColumns={'Frequency_Hz','Rate_HzPerS','ApproxScore','Branch', ...
        'ExactScore','ApproxRank','ExactRank','Confirmed'};
    e.refinementColumns={'StartF_Hz','StartRate_HzPerS','OutputF_Hz','OutputRate_HzPerS', ...
        'ExactScore','ExitFlag','OutputValid','Converged'};
end

function [frames,tc,ts,w,n]=prepare_frames(y,c)
    n=numel(y); starts=0:c.hop:n-c.L;
    w=.5-.5*cos(2*pi*(0:c.L-1)'/(c.L-1));
    if numel(starts)<3; frames=[]; tc=[]; ts=[]; return; end
    indices=bsxfun(@plus,(0:c.L-1)',starts)+1;
    frames=bsxfun(@times,y(indices),w);
    reference=(n-1)/(2*c.scfFs);
    tc=(indices-1)/c.scfFs-reference; ts=starts/c.scfFs-reference;
end

function spec=fsm_spectrum(frames,tc,f_candidate,a,c,w)
    % alpha=2*f_candidate；FFT输出轴为f_spec。相位psi=2*pi*f*t+pi*a*t^2；交叉乘积中的循环相位是2*psi。
    % a=0退化为普通移频FSM；a非零是多项式相位补偿的FSM，不声称原chirp平稳。
    phase=2*pi*f_candidate*tc+pi*a*tc.^2;
    Xp=fft(frames.*exp(-1j*phase),[],1);
    Xm=fft(frames.*exp( 1j*phase),[],1);
    spec=fftshift(mean(Xp.*conj(Xm),2))/(c.scfFs*sum(w.^2));
end

function J=fsm_score(frames,tc,f,a,c,w,band)
    spec=fsm_spectrum(frames,tc,f,a,c,w);
    % 先复数频率平滑，最后取模；17点平滑在f=0处的值作为统一精确指标。
    J=abs(mean(spec(c.L/2+1+band)))^2;
end

function report=export_figures(figs,folder,forceFailure)
    report=struct('pngSuccess',0,'pngFailure',0,'figSuccess',0,'figFailure',0,'details',{{}});
    for j=1:numel(figs)
        h=figs{j}; png=fullfile(folder,sprintf('figure_%02d.png',j));
        fig=fullfile(folder,sprintf('figure_%02d.fig',j)); message='';
        if isempty(h) || ~isgraphics(h,'figure')
            report.pngFailure=report.pngFailure+1; report.figFailure=report.figFailure+1;
            report.details{end+1}=sprintf('图%d句柄不存在或已经关闭。',j); continue;
        end
        try
            number=get(h,'Number'); % 窗口可能在isgraphics检查后被用户关闭
        catch ME
            report.pngFailure=report.pngFailure+1; report.figFailure=report.figFailure+1;
            report.details{end+1}=ME.message; continue;
        end
        % print回退明确使用图窗编号，不传失效的对象句柄
        try
            if forceFailure; error('FSM:Injected','有意模拟FIG导出失败。'); end
            savefig(h,fig); report.figSuccess=report.figSuccess+1;
        catch ME; report.figFailure=report.figFailure+1; message=[message ' FIG:' ME.message]; end
        pngOK=false; method='none';
        if exist('exportgraphics','file')~=0 || exist('exportgraphics','builtin')~=0
            try
                if forceFailure; error('FSM:Injected','有意模拟exportgraphics失败。'); end
                exportgraphics(h,png,'Resolution',150); pngOK=true; method='exportgraphics';
            catch ME; message=[message ' exportgraphics:' ME.message]; end
        end
        if ~pngOK
            try
                assert(isgraphics(h,'figure'),'图窗在导出期间已关闭。');
                if forceFailure; error('FSM:Injected','有意模拟print失败。'); end
                print(number,png,'-dpng','-r150'); pngOK=true; method='print_by_number';
            catch ME; message=[message ' print:' ME.message]; end
        end
        if pngOK; report.pngSuccess=report.pngSuccess+1;
        else; report.pngFailure=report.pngFailure+1; end
        report.details{end+1}=sprintf('图%d PNG=%d method=%s: %s',j,pngOK,method,message);
    end
end

function folder=unique_directory(base)
    folder=base; index=1;
    while exist(folder,'dir') || exist(folder,'file')
        folder=[base '_' num2str(index)]; index=index+1;
    end
end

function code=make_gold()
    a=mseq([0 2]); b=mseq([0 2 3 4]);
    code=2*double(xor(logical(a),circshift(logical(b),[0 1])))-1;
end

function bits=mseq(taps)
    s=ones(1,36); for i=1:31; s(i+5)=mod(sum(s(i+taps)),2); end
    assert(isequal(s(32:36),s(1:5)) && sum(s(1:31))==16);
    bits=s(1:31);
end

