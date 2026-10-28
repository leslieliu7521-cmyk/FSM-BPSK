%% FSM第一层固定窗长基线：真实变化率700Hz/s（默认精简150次含噪估计）
% 原v4/700程序映射：第一层=分支离散候选→去重→前10候选精确FSM确认→最佳离散点。
% 第二层=多起点连续优化；本文件完全移除第二层，仅输出第一层，不是隐藏曲线。
% 核心始终为频移FFT→谱交叉乘积→复数跨帧平均→复数频率平滑→最后取模。
% 接收机只获得接收数据、固定先验和算法参数；真值仅用于信道与评价。
% 直接运行：run('FSM_BPSK_stage1_700.m')；需要MATLAB R2018b+及Signal Processing Toolbox。
% 自定义：stage1_options=struct('trials',100,'makePlot',true); run(...)
% 恢复：stage1_options=struct('resumeFile','实际目录/checkpoint.mat'); run(...)
% 不保存高采样率波形；默认不画图/不保存FIG，不做真值辅助窗长选择。
clearvars -except stage1_options; clc;
if ~exist('stage1_options','var'); stage1_options=struct(); end
sourceText=fileread([mfilename('fullpath') '.m']);
version='FSM_stage1_700_compact_20261028_1';

%% 一、保持原信号参数、Eb/N0约定、固定先验、离散搜索定义
cfg=struct('Rb',1e5,'Rc',31e5,'fs',24.8e6,'fc',1e5,'longT',.1, ...
    'scfFs',8e5,'L',128,'hop',64,'smoothingBins',17,'fRange',[50e3 150e3], ...
    'rateRange',[-1000 1000],'rateStep',100,'branchPeaks',3,'confirmCandidates',10, ...
    'dedupFrequencyFactor',.5,'dedupRateFactor',2,'fTolerance',20,'rateTolerance',200, ...
    'fdTest',45000,'filterN',10,'filterBeta',5);
opts=struct('mode','compact','trials',50,'ebN0',[-10 -7 -4],'durations',.1,'rates',700, ...
    'seed',20261014,'windowComparison',false,'rateComparison',false, ...
    'syncComparison',false,'bandwidthComparison',false,'candidateComparison',false, ...
    'makePlot',false,'runNoiseless',true,'verifyOptimization',true, ...
    'checkpointEveryTrials',10,'maxNewBaseTrials',Inf,'failureDiagnosticsPerCondition',2, ...
    'saveFrameSequence',false);
keys=fieldnames(stage1_options);
for ii=1:numel(keys)
    assert(isfield(opts,keys{ii}) || strcmp(keys{ii},'resumeFile'),'未知选项：%s',keys{ii});
    opts.(keys{ii})=stage1_options.(keys{ii});
end
assert(any(strcmp(opts.mode,{'compact','debug','full','selftest'})));
if ~isfield(stage1_options,'trials')
    if strcmp(opts.mode,'debug'); opts.trials=100;
    elseif strcmp(opts.mode,'full'); opts.trials=1000; end
end
assert(isfinite(opts.trials) && opts.trials>=1 && opts.trials==floor(opts.trials));
if strcmp(opts.mode,'full'); assert(opts.trials==1000,'full每条件必须1000次。'); end
assert(all(isfinite(opts.ebN0)) && ~isempty(opts.ebN0));
assert(all(ismember(opts.durations,[.001 .01 .05 .1])) && ~isempty(opts.durations));
assert(all(isfinite(opts.rates)) && ~isempty(opts.rates));
assert(opts.checkpointEveryTrials>=1 && opts.checkpointEveryTrials==floor(opts.checkpointEveryTrials));
assert(opts.maxNewBaseTrials>=1 && (isinf(opts.maxNewBaseTrials) || opts.maxNewBaseTrials==floor(opts.maxNewBaseTrials)));
assert(opts.failureDiagnosticsPerCondition>=0 && opts.failureDiagnosticsPerCondition==floor(opts.failureDiagnosticsPerCondition));
assert(exist('resample','file')~=0,'需要Signal Processing Toolbox的resample。');

%% 二、配置/随机种子/完成标记检查点；科学配置不能混合恢复
if isfield(opts,'resumeFile') && ~isempty(opts.resumeFile)
    loaded=load(opts.resumeFile,'state'); state=loaded.state;
    assert(strcmp(state.version,version) && strcmp(state.sourceText,sourceText),'检查点版本或源文件不匹配。');
    outputDir=fileparts(opts.resumeFile); if isempty(outputDir); outputDir=pwd; end
    saved=state.opts;
    supplied=fieldnames(stage1_options);
    for ii=1:numel(supplied)
        key=supplied{ii};
        if ~ismember(key,{'resumeFile','makePlot','checkpointEveryTrials','maxNewBaseTrials'})
            assert(isfield(saved,key) && isequaln(stage1_options.(key),saved.(key)), ...
                '恢复时不能更改科学选项%s，请另建实验。',key);
        end
    end
    for ii={'makePlot','checkpointEveryTrials','maxNewBaseTrials'}
        key=ii{1}; if isfield(stage1_options,key); saved.(key)=opts.(key); end
    end
    opts=saved; cfg=state.cfg; state.opts=opts; rng(state.rngState);
else
    rng(opts.seed,'twister');
    conditions=make_conditions(cfg,opts);
    outputDir=unique_directory(fullfile(pwd,['FSM_stage1_700_' datestr(now,'yyyymmdd_HHMMSS')])); mkdir(outputDir);
    state=struct('version',version,'sourceText',sourceText,'cfg',cfg,'opts',opts, ...
        'conditions',conditions,'seeds',opts.seed+(1:opts.trials),'done',false(numel(conditions),opts.trials), ...
        'records',{cell(numel(conditions),opts.trials)},'noiseless',table(), ...
        'noiselessDone',false,'benchmarks',{{}},'summary',table(),'rngState',rng, ...
        'failureDiagnostics',{cell(numel(conditions),1)},'lastException','', ...
        'wallSeconds',0,'ioSeconds',0,'continuousRefinementCalls',0);
    atomic_save(state,fullfile(outputDir,'checkpoint.mat'));
end
fprintf('输出：%s\n模式%s；%d条件×%d次=%d次含噪估计，跨条件配对，不合并为独立样本。\n', ...
    outputDir,opts.mode,numel(state.conditions),opts.trials,numel(state.conditions)*opts.trials);
fprintf('第二层已移除；f0门限20Hz；原采样SNR(-7dB Eb/N0)=%.3fdB。\n',-7+10*log10(2*cfg.Rb/cfg.fs));
checks=internal_checks(cfg); state.checks=checks;
if strcmp(opts.mode,'selftest')
    state.summary=table(); atomic_save(state,fullfile(outputDir,'results.mat')); return;
end
runClock=tic; oldWall=state.wallSeconds;
context=make_context(cfg); % 时间轴/已知Gold PN只构建一次，不按条件重复
if opts.runNoiseless && ~state.noiselessDone
    [state.noiseless,benchmark]=noiseless_checks(cfg,opts,context);
    state.benchmarks{end+1}=benchmark;
    state.noiselessDone=true; state.rngState=rng;
    atomic_save(state,fullfile(outputDir,'checkpoint.mat'));
end

%% 三、TrialID外循环：同一基础比特/初相/频偏/标准噪声共用于SNR/窗长/率对照
% 配对组内先滤完整100ms长记录再截中心段；先前版本每条件重复生成同seed记录。
% 此处每个TrialID只生成一次，利用线性FIR分别滤干净信号与标准噪声，再按Eb/N0缩放。
% resample(x+sigma*n)=resample(x)+sigma*resample(n)，允许合理浮点误差；首批原生核对。
newTrials=0; preKeys={state.conditions.PrepKey}; groupKeys=unique(preKeys,'stable');
for trialID=1:opts.trials
    if all(state.done(:,trialID)); continue; end
    if newTrials>=opts.maxNewBaseTrials; break; end
    try
        baseClock=tic; base=generate_base(cfg,context,state.seeds(trialID)); baseSeconds=toc(baseClock);
        todoCount=sum(~state.done(:,trialID));
        noiseCache=struct('key',{},'value',{});
        for gg=1:numel(groupKeys)
            ids=find(strcmp(preKeys,groupKeys{gg}) & ~state.done(:,trialID)');
            if isempty(ids); continue; end
            prototype=state.conditions(ids(1));
            prepClock=tic;
            % 真值RateTruth只进入发射信道；接收配置receiver_config不含它。
            clean=base.spread.*cos(2*pi*(cfg.fc+base.fd)*context.t ...
                +pi*prototype.RateTruth*context.t.^2+base.phase);
            signalPower=mean(clean.^2);
            pnRx=local_pn(context,prototype.PNPhaseChip,prototype.LocalPNRatePPM,cfg);
            [up,down]=rat(prototype.ScfFs/cfg.fs,1e-12);
            cleanLong=resample((clean.*pnRx)',up,down,cfg.filterN,cfg.filterBeta);
            noiseKey=sprintf('%g_%g_%g',prototype.ScfFs,prototype.PNPhaseChip,prototype.LocalPNRatePPM);
            ni=find(strcmp({noiseCache.key},noiseKey),1);
            if isempty(ni)
                noiseLong=resample((base.noise.*pnRx)',up,down,cfg.filterN,cfg.filterBeta);
                noiseCache(end+1)=struct('key',noiseKey,'value',noiseLong); %#ok<AGROW>
            else; noiseLong=noiseCache(ni).value; end
            prepSeconds=toc(prepClock);
            for kk=ids(:)'
                rowPrepClock=tic; current=state.conditions(kk); receiverCfg=receiver_config(cfg,current);
                variance=signalPower*cfg.fs/(2*cfg.Rb*10^(current.EbN0/10));
                n=round(current.Duration*current.ScfFs); first=floor((numel(cleanLong)-n)/2)+1;
                y=cleanLong(first:first+n-1)+sqrt(variance)*noiseLong(first:first+n-1);
                localReference=(first-1)/current.ScfFs+(n-1)/(2*current.ScfFs);
                rowPrepSeconds=toc(rowPrepClock);
                searchClock=tic; estimate=stage1_receiver(y,receiverCfg,true); searchSeconds=toc(searchClock);
                % 优化核对仅最初两个基础试验的首个条件；额外开销另记，不算常规搜索计时。
                if opts.verifyOptimization && trialID<=2 && kk==1
                    verifyClock=tic;
                    referencePrep=tic;
                    referenceLong=resample(((clean+sqrt(variance)*base.noise).*pnRx)',up,down,cfg.filterN,cfg.filterBeta);
                    referencePrepSeconds=toc(referencePrep);
                    refClock=tic; referenceEstimate=stage1_reference(referenceLong(first:first+n-1),receiverCfg); refSeconds=toc(refClock);
                    compare_estimates(estimate,referenceEstimate);
                    state.benchmarks{end+1}=struct('Kind','paired_native_reference_check', ...
                        'TrialID',trialID,'Condition',kk,'ReferencePreprocess_s',referencePrepSeconds, ...
                        'ReferenceSearch_s',refSeconds,'CachedSearch_s',searchSeconds,'Verified',true, ...
                        'ValidationWall_s',toc(verifyClock));
                    clear referenceLong;
                end
                evaluateClock=tic;
                metrics=evaluate_record(estimate,base,context,current,cfg,trialID,localReference,n);
                metrics.Search_s=searchSeconds;
                metrics.Generation_s=estimate.generationSeconds; metrics.Confirmation_s=estimate.confirmationSeconds;
                metrics.PreprocessAllocated_s=baseSeconds/todoCount+prepSeconds/numel(ids)+rowPrepSeconds;
                metrics.Evaluation_s=toc(evaluateClock);
                metrics.Compute_s=metrics.PreprocessAllocated_s+metrics.Search_s+metrics.Evaluation_s;
                metrics.NoiseVariance=variance;
                record=struct('metrics',metrics,'search',estimate);
                % 默认仅存一致性摘要。必要时显式开启每帧复谱序列，不存高采样率波形。
                if ~opts.saveFrameSequence; record.search.frameSequence=[]; end
                state.records{kk,trialID}=record; state.done(kk,trialID)=true;
                if ~metrics.F0Success && numel(state.failureDiagnostics{kk})<opts.failureDiagnosticsPerCondition
                    diagnostic=struct('TrialID',trialID,'BaseSeed',base.seed,'Estimate',estimate, ...
                        'F0Error_Hz',metrics.F0Error_Hz,'RateError_HzPerS',metrics.RateError_HzPerS);
                    if isempty(state.failureDiagnostics{kk}); state.failureDiagnostics{kk}={diagnostic};
                    else; state.failureDiagnostics{kk}{end+1}=diagnostic; end
                end
            end
            clear clean cleanLong noiseLong pnRx y;
        end
        newTrials=newTrials+1;
        fprintf('基础试验%d/%d完成；已估计%d/%d。\n',trialID,opts.trials,sum(state.done(:)),numel(state.done));
        if trialID==min(2,opts.trials)
            show_timing_estimate(state); % 实测后才估算本机剩余时间，不预设分钟数
        end
        if mod(newTrials,opts.checkpointEveryTrials)==0 || all(state.done(:))
            state.rngState=rng; state.wallSeconds=oldWall+toc(runClock);
            saveClock=tic; atomic_save(state,fullfile(outputDir,'checkpoint.mat'));
            state.ioSeconds=state.ioSeconds+toc(saveClock);
        end
        clear base noiseCache;
    catch ME
        state.lastException=getReport(ME); state.rngState=rng; state.wallSeconds=oldWall+toc(runClock);
        atomic_save(state,fullfile(outputDir,'checkpoint.mat')); rethrow(ME);
    end
end

%% 四、数值先保存，最后可选一幅汇总PNG；关闭图窗/导出失败不影响统计
state.trials=trial_table(state); state.summary=summarize_conditions(state);
state.complete=all(state.done(:)); state.rngState=rng; state.wallSeconds=oldWall+toc(runClock);
atomic_save(state,fullfile(outputDir,'results.mat'));
atomic_save(state,fullfile(outputDir,'checkpoint.mat'));
state.csvExport=export_tables(state,outputDir);
state.plotExport=struct('success',0,'failure',0,'message','Disabled');
if opts.makePlot; state.plotExport=export_summary_plot(state,outputDir); end
try
    atomic_save(state,fullfile(outputDir,'results.mat'));
    atomic_save(state,fullfile(outputDir,'checkpoint.mat'));
catch ME; warning('FSM:FinalSave','最终报告未更新，导图前MAT仍保留：%s',ME.message); end
show=state.summary;
disp(show(:,{'Condition','Count','F0SuccessRate','F0Lo','F0Hi','F0RMSE_Hz', ...
    'F0Bias_Hz','F0Std_Hz','AbsF0Median_Hz','AbsF0P90_Hz','AbsF0P95_Hz','InvalidCount','F0FailureCount'}));
fprintf('数值CSV成功%d/失败%d；PNG成功%d/失败%d。\n',state.csvExport.success,state.csvExport.failure, ...
    state.plotExport.success,state.plotExport.failure);
fprintf('实际%d次/条件；默认50次仅作初步统计，尾部分位数需谨慎，不承诺总体成功率。\n',opts.trials);

%% 五、条件表：默认仅700Hz/s、100ms、三个Eb/N0；重复700自动去重
function C=make_conditions(cfg,o)
    durations=unique(o.durations(:)','stable'); rates=unique(o.rates(:)','stable'); eb=unique(o.ebN0(:)','stable');
    if o.windowComparison; durations=unique([durations .001 .01 .05 .1],'stable'); end
    if o.rateComparison; rates=unique([rates 650 700 750],'stable'); end
    phases=0; ppms=0; fsout=cfg.scfFs; confirms=cfg.confirmCandidates;
    if o.syncComparison; phases=[0 .125 0]; ppms=[0 0 5]; end
    if o.bandwidthComparison; fsout=[cfg.scfFs 16e5]; end
    if o.candidateComparison; confirms=[5 10 Inf]; end
    C=struct('Name',{},'EbN0',{},'Duration',{},'RateTruth',{},'PNPhaseChip',{}, ...
        'LocalPNRatePPM',{},'ScfFs',{},'ConfirmCandidates',{},'PrepKey',{});
    for a=rates
        for si=1:numel(phases)
            for fs=fsout
                for T=durations
                    for count=confirms
                        for snr=eb
                            name=sprintf('E%g_T%g_A%g_P%g_R%g_Fs%g_K%g',snr,1000*T,a,phases(si),ppms(si),fs,count);
                            key=sprintf('A%g_P%g_R%g_Fs%g',a,phases(si),ppms(si),fs);
                            C(end+1)=struct('Name',name,'EbN0',snr,'Duration',T,'RateTruth',a, ...
                                'PNPhaseChip',phases(si),'LocalPNRatePPM',ppms(si),'ScfFs',fs, ...
                                'ConfirmCandidates',count,'PrepKey',key); %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
    assert(numel(unique({C.Name}))==numel(C),'重复条件会重复计数。');
end
function r=receiver_config(c,d)
    % 白名单没有RateTruth/频偏/比特/初相/PN失配真值，不能利用其筛选候选。
    names={'smoothingBins','fRange','rateRange','rateStep','branchPeaks', ...
        'dedupFrequencyFactor','dedupRateFactor'};
    r=struct(); for j=1:numel(names); r.(names{j})=c.(names{j}); end
    r.scfFs=d.ScfFs; r.L=round(c.L*d.ScfFs/c.scfFs); r.hop=r.L/2;
    r.confirmCandidates=d.ConfirmCandidates;
end

%% 六、基础长记录：原Gold码、矩形采样、实通带与原随机消耗顺序
function x=make_context(c)
    nb=round(c.Rb*c.longT); Ns=c.fs/c.Rc; assert(Ns==round(Ns));
    code=make_gold(); N=nb*31*Ns;
    x.code=code; x.pn=repelem(repmat(code,1,nb),Ns);
    x.reference=(N-1)/(2*c.fs); x.t=(0:N-1)/c.fs-x.reference;
    x.N=N; x.nb=nb; x.Ns=Ns;
end
function b=generate_base(c,x,seed)
    stream=RandStream('mt19937ar','Seed',seed);
    bits=2*randi(stream,[0 1],1,x.nb)-1;
    b.phase=2*pi*rand(stream); b.fd=(2*rand(stream)-1)*c.fdTest;
    unusedRateDraw=rand(stream); %#ok<NASGU> 保留原随机变化率的抽样位置，真实率由条件施加
    b.noise=randn(stream,1,x.N); b.spread=repelem(bits,31*x.Ns).*x.pn; b.seed=seed;
end
function pn=local_pn(x,phaseChip,pnRatePPM,c)
    % 仅本地PN索引失配，不是ADC时钟误差；相位分辨率1/8码片，不声称连续时延。
    if phaseChip==0 && pnRatePPM==0; pn=x.pn; return; end
    ci=floor((0:x.N-1)*(1+pnRatePPM*1e-6)/(c.fs/c.Rc)+phaseChip);
    pn=x.code(mod(ci,31)+1);
end

%% 七、第一层接收机：预计算只复用数据无关项，候选数/网格/评分不变
function e=stage1_receiver(y,c,useCache)
    e=empty_estimate(); generationClock=tic;
    if any(~isfinite(y(:))); return; end
    p=frame_plan(numel(y),c,useCache);
    if isempty(p.indices); return; end
    frames=bsxfun(@times,y(p.indices),p.w); X=fft(frames,[],1);
    generated=zeros(numel(p.branches)*c.branchPeaks,4); cursor=0;
    for j=1:numel(p.branches)
        b=p.branches(j);
        cross=X(b.plus,:).*conj(X(b.minus,:));
        v=mean(cross,1).*b.phase;
        bank=bsxfun(@times,v,p.rateBank);
        power=abs(fftshift(fft(bank,p.nt,2),2)).^2;
        rows=branch_candidates(power(:,b.valid),b.f,p.rates,c.branchPeaks,p.dF,p.dA,b.k);
        generated(cursor+(1:size(rows,1)),:)=rows; cursor=cursor+size(rows,1);
    end
    generated=generated(1:cursor,:); genSeconds=toc(generationClock);
    confirmationClock=tic;
    e=confirm_select(generated,frames,p.tc,p.w,p.band,c,p.T);
    e.generationSeconds=genSeconds; e.confirmationSeconds=toc(confirmationClock);
end
function p=frame_plan(n,c,useCache)
    % 缓存中只有固定窗/采样索引/时间/先验栅格/补偿相位，无接收数据或发送真值。
    persistent cacheKeys cacheValues
    if isempty(cacheKeys); cacheKeys={}; cacheValues={}; end
    key=sprintf('%.17g_', [n c.scfFs c.L c.hop c.smoothingBins c.fRange c.rateRange ...
        c.rateStep c.dedupFrequencyFactor c.dedupRateFactor]);
    if useCache
        found=find(strcmp(cacheKeys,key),1);
        if ~isempty(found); p=cacheValues{found}; return; end
    end
    p=struct(); starts=0:c.hop:n-c.L;
    p.w=.5-.5*cos(2*pi*(0:c.L-1)'/(c.L-1));
    p.indices=[];
    if numel(starts)<3; return; end
    p.indices=bsxfun(@plus,(0:c.L-1)',starts)+1;
    p.tc=(p.indices-1)/c.scfFs-(n-1)/(2*c.scfFs);
    p.ts=starts/c.scfFs-(n-1)/(2*c.scfFs);
    p.T=n/c.scfFs; p.dF=c.dedupFrequencyFactor/p.T; p.dA=c.dedupRateFactor/p.T^2;
    p.band=-(c.smoothingBins-1)/2:(c.smoothingBins-1)/2;
    p.rates=unique([c.rateRange(1):c.rateStep:c.rateRange(2),c.rateRange(2)]);
    p.nt=2^nextpow2(2*numel(starts));
    da=(-p.nt/2:p.nt/2-1)*(c.scfFs/c.hop)/p.nt; df=c.scfFs/c.L;
    p.rateBank=exp(-1j*2*pi*p.rates(:)*p.ts.^2);
    p.branches=struct('k',{},'plus',{},'minus',{},'phase',{},'valid',{},'f',{});
    for k=floor(c.fRange(1)/df):ceil(c.fRange(2)/df)
        fb=k*df; f=fb+da/2;
        valid=da>=-df & da<df & f>=c.fRange(1) & f<=c.fRange(2);
        p.branches(end+1)=struct('k',k,'plus',mod(p.band+k,c.L)+1,'minus',mod(p.band-k,c.L)+1, ...
            'phase',exp(-1j*4*pi*fb*p.ts),'valid',valid,'f',f(valid)); %#ok<AGROW>
    end
    if useCache; cacheKeys{end+1}=key; cacheValues{end+1}=p; end
end
function e=stage1_reference(y,c)
    % 对照实现：逐次重建窗/索引/各分支相位和率补偿，直接对应原v4第一层公式。
    % 只用于首批正确性核对与测速，绝不执行连续参数优化。
    e=empty_estimate(); if any(~isfinite(y(:))); return; end
    n=numel(y); starts=0:c.hop:n-c.L; if numel(starts)<3; return; end
    w=.5-.5*cos(2*pi*(0:c.L-1)'/(c.L-1));
    indices=bsxfun(@plus,(0:c.L-1)',starts)+1;
    frames=bsxfun(@times,y(indices),w);
    tc=(indices-1)/c.scfFs-(n-1)/(2*c.scfFs); ts=starts/c.scfFs-(n-1)/(2*c.scfFs);
    T=n/c.scfFs; dF=c.dedupFrequencyFactor/T; dA=c.dedupRateFactor/T^2;
    band=-(c.smoothingBins-1)/2:(c.smoothingBins-1)/2;
    X=fft(frames,[],1); df=c.scfFs/c.L; nt=2^nextpow2(2*numel(ts));
    deltaAlpha=(-nt/2:nt/2-1)*(c.scfFs/c.hop)/nt;
    rates=unique([c.rateRange(1):c.rateStep:c.rateRange(2),c.rateRange(2)]);
    generated=zeros(0,4);
    for k=floor(c.fRange(1)/df):ceil(c.fRange(2)/df)
        fb=k*df;
        cross=X(mod(band+k,c.L)+1,:).*conj(X(mod(band-k,c.L)+1,:));
        v=mean(cross,1).*exp(-1j*4*pi*fb*ts);
        bank=bsxfun(@times,v,exp(-1j*2*pi*rates(:)*ts.^2));
        power=abs(fftshift(fft(bank,nt,2),2)).^2; f=fb+deltaAlpha/2;
        valid=deltaAlpha>=-df & deltaAlpha<df & f>=c.fRange(1) & f<=c.fRange(2);
        rows=branch_candidates(power(:,valid),f(valid),rates,c.branchPeaks,dF,dA,k);
        generated=[generated;rows]; %#ok<AGROW>
    end
    e=confirm_select(generated,frames,tc,w,band,c,T);
end
function e=confirm_select(G,frames,tc,w,band,c,T)
    e=empty_estimate(); e.generated=G; e.frameCount=size(frames,2);
    e.dedupScale=[c.dedupFrequencyFactor/T c.dedupRateFactor/T^2];
    if isempty(G); return; end
    [~,order]=sort(G(:,3),'descend'); sorted=G(order,:);
    keep=distinct_indices(sorted(:,1:2),e.dedupScale(1),e.dedupScale(2));
    distinct=sorted(keep,:); K=size(distinct,1); count=min(K,c.confirmCandidates);
    C=[distinct,NaN(K,1),(1:K)',NaN(K,1),zeros(K,1)];
    frameCache=cell(count,1);
    for j=1:count
        % 去重后不同候选各评分一次，无重复全谱FFT或后验真值探查。
        [C(j,5),frameCache{j}]=fsm_score(frames,tc,C(j,1),C(j,2),c,w,band); C(j,8)=1;
    end
    valid=find(C(:,8)>0 & isfinite(C(:,5)) & C(:,5)>=0);
    e.candidates=C; e.generatedSortOrder=order; e.distinctGeneratedIndices=order(keep);
    if isempty(valid); return; end
    [~,ranking]=sort(C(valid,5),'descend'); ranking=valid(ranking);
    C(ranking,7)=(1:numel(ranking))'; best=ranking(1);
    e.candidates=C; e.exactRanking=ranking; e.selectedCandidate=best;
    e.fLocal=C(best,1); e.rate=C(best,2); e.score=C(best,5);
    e.valid=all(isfinite([e.fLocal e.rate e.score])) && e.fLocal>=c.fRange(1) ...
        && e.fLocal<=c.fRange(2) && e.rate>=c.rateRange(1) && e.rate<=c.rateRange(2);
    e.boundary=min(abs(e.fLocal-c.fRange))<=min(e.dedupScale(1),.01*diff(c.fRange)) ...
        || min(abs(e.rate-c.rateRange))<=min(e.dedupScale(2),.01*diff(c.rateRange));
    e.frameSequence=frameCache{best};
    % q_j为第一层估计补偿后的第j帧复循环谱，先复数频率平滑，不逐帧取模再平均。
    % rho=|sum_j q_j|/sum_j|q_j|，0--1；重叠帧相关，不解释为独立概率/置信度。
    % 指标可能在错误峰也很高；本版不根据rho选择窗长或拒判，不预设其有效。
    denominator=sum(abs(e.frameSequence));
    if isfinite(denominator) && denominator>0 && all(isfinite(e.frameSequence))
        e.coherence=min(1,max(0,abs(sum(e.frameSequence))/denominator)); e.coherenceValid=true;
    end
end
function e=empty_estimate()
    e=struct('fLocal',NaN,'rate',NaN,'score',NaN,'valid',false,'boundary',false, ...
        'generated',zeros(0,4),'candidates',zeros(0,8),'generatedSortOrder',[], ...
        'distinctGeneratedIndices',[],'exactRanking',[],'selectedCandidate',0, ...
        'frameCount',0,'frameSequence',[],'coherence',NaN,'coherenceValid',false, ...
        'generationSeconds',0,'confirmationSeconds',0,'dedupScale',[NaN NaN]);
    e.generatedColumns={'Frequency_Hz','Rate_HzPerS','ApproxScore','FFTBranch'};
    e.candidateColumns={'Frequency_Hz','Rate_HzPerS','ApproxScore','FFTBranch', ...
        'ExactScore','ApproxRankAfterDedup','ExactRank','Confirmed'};
end

%% 八、保留原候选抽取/去重和精确频域FSM定义
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

function [spec,q]=fsm_spectrum(frames,tc,f,a,c,w)
    % alpha=2*f。相位psi=2*pi*f*t+pi*a*t^2；交叉乘积中的循环相位是2*psi。
    % a=0退化为普通移频FSM；a非零是多项式相位补偿的FSM，不声称原chirp平稳。
    phase=2*pi*f*tc+pi*a*tc.^2;
    Xp=fft(frames.*exp(-1j*phase),[],1);
    Xm=fft(frames.*exp( 1j*phase),[],1);
    cross=Xp.*conj(Xm);
    spec=fftshift(mean(cross,2))/(c.scfFs*sum(w.^2));
    bins=-(c.smoothingBins-1)/2:(c.smoothingBins-1)/2;
    q=mean(cross(mod(bins,c.L)+1,:),1)/(c.scfFs*sum(w.^2));
end

function [J,q]=fsm_score(frames,tc,f,a,c,w,band)
    [spec,q]=fsm_spectrum(frames,tc,f,a,c,w);
    % 先复数频率平滑，最后取模；17点平滑在f=0处的值作为统一精确指标。
    J=abs(mean(spec(c.L/2+1+band)))^2;
end

function [rmse,maxError]=track_error(df,da,jerk,n,fs)
    if any(~isfinite([df da jerk])) || n<2; rmse=Inf; maxError=Inf; return; end
    % u=t-t_ref，误差e(u)=df+da*u-jerk*u^2/2；参考时刻取整段采样中点。
    % 下式精确计算整段均匀采样点的均方误差，避免每次显式生成几万点轨迹。
    m2=(n^2-1)/(12*fs^2); m4=(n^2-1)*(3*n^2-7)/(240*fs^4);
    rmse=sqrt(max(0,df^2+da^2*m2+.25*jerk^2*m4-df*jerk*m2));
    bound=(n-1)/(2*fs); points=[-bound bound];
    if jerk~=0 && abs(da/jerk)<=bound; points(end+1)=da/jerk; end
    % 连续区间最大绝对误差（端点及二次函数顶点），不因稀疏画图漏掉峰值。
    maxError=max(abs(df+da*points-.5*jerk*points.^2));
end

function [rmse,bias,sd]=error_stats(v)
    if isempty(v); rmse=NaN; bias=NaN; sd=NaN;
    elseif any(~isfinite(v)); rmse=Inf; bias=NaN; sd=NaN;
    else; rmse=sqrt(mean(v.^2)); bias=mean(v); sd=std(v,0); end
    % 一旦存在非有限误差，完整样本RMSE报Inf、偏差/标准差报NaN；另有失效计数。
    % 不在删去失败样本后报告一个看似更好的有限RMSE。
end

function [lo,hi]=wilson(success,n)
    if n==0; lo=NaN; hi=NaN; return; end
    z=1.96; p=success/n; den=1+z^2/n;
    center=(p+z^2/(2*n))/den; half=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/den;
    lo=max(0,center-half); hi=min(1,center+half);
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


%% 九、统一时间参考的独立评价；频率主指标与率/联合指标分开
function m=evaluate_record(e,b,x,d,c,id,localRef,n)
    m=struct(); m.TrialID=id; m.BaseSeed=b.seed;
    m.EbN0_dB=d.EbN0; m.Duration_ms=1000*d.Duration; m.TrueRate_HzPerS=d.RateTruth;
    m.CommonReference_s=x.reference; m.LocalReference_s=localRef;
    m.TrueF0_Hz=c.fc+b.fd; % 共同参考为同一基础长记录原始采样中点
    m.LocalF_Hz=e.fLocal; m.EstimatedRate_HzPerS=e.rate;
    % 关键：估计坐标换算只使用估计率e.rate，绝不使用d.RateTruth。
    m.F0_Hz=e.fLocal+e.rate*(x.reference-localRef);
    m.F0Error_Hz=m.F0_Hz-m.TrueF0_Hz; m.RateError_HzPerS=e.rate-d.RateTruth;
    m.Valid=e.valid && isfinite(m.F0_Hz);
    m.F0Success=m.Valid && abs(m.F0Error_Hz)<=c.fTolerance; % 原20Hz门限，非联合成功
    m.RateSuccess=m.Valid && abs(m.RateError_HzPerS)<=c.rateTolerance; % 原200Hz/s，单列
    m.JointSuccess=m.F0Success && m.RateSuccess;
    % 真值仅在评价里计算对应局部中点频率，检查整段轨迹；不回写搜索输出。
    trueLocal=m.TrueF0_Hz+d.RateTruth*(localRef-x.reference);
    [m.TrackRMSE_Hz,m.TrackMaxError_Hz]=track_error(e.fLocal-trueLocal,m.RateError_HzPerS,0,n,d.ScfFs);
    m.BoundaryHit=e.boundary; m.ExactScore=e.score;
    m.FrameCoherence=e.coherence; m.CoherenceValid=e.coherenceValid; m.FrameCount=e.frameCount;
    m.GeneratedCount=size(e.generated,1); m.DistinctCount=size(e.candidates,1);
    m.ConfirmedCount=sum(e.candidates(:,8)>0);
    m.ConfirmedNonfiniteCount=sum(e.candidates(:,8)>0 & ~isfinite(e.candidates(:,5)));
    m.ContinuousRefinementCalls=0;
    m.Search_s=NaN; m.Generation_s=NaN; m.Confirmation_s=NaN;
    m.PreprocessAllocated_s=NaN; m.Evaluation_s=NaN; m.Compute_s=NaN;
    m.NoiseVariance=NaN; m.SampleSNR_dB=d.EbN0+10*log10(2*c.Rb/c.fs);
end
function T=trial_table(state)
    rows={};
    for kk=1:numel(state.conditions)
        ids=find(state.done(kk,:));
        for id=ids
            m=state.records{kk,id}.metrics; m.ConditionIndex=kk;
            m.Condition=string(state.conditions(kk).Name);
            m.ReceiverFs_Hz=state.conditions(kk).ScfFs;
            m.PNPhaseChip=state.conditions(kk).PNPhaseChip;
            m.LocalPNRatePPM=state.conditions(kk).LocalPNRatePPM;
            m.ConfirmLimit=state.conditions(kk).ConfirmCandidates;
            rows{end+1,1}=m; %#ok<AGROW>
        end
    end
    if isempty(rows); T=table(); else; T=struct2table(vertcat(rows{:})); end
end
function S=summarize_conditions(state)
    T=state.trials; rows=cell(numel(state.conditions),1);
    for kk=1:numel(state.conditions)
        d=state.conditions(kk); if isempty(T); q=table(); n=0;
        else; q=T(T.ConditionIndex==kk,:); n=height(q); end
        r=struct('Condition',string(d.Name),'EbN0_dB',d.EbN0,'Duration_ms',1000*d.Duration, ...
            'TrueRate_HzPerS',d.RateTruth,'ReceiverFs_Hz',d.ScfFs,'PNPhaseChip',d.PNPhaseChip, ...
            'LocalPNRatePPM',d.LocalPNRatePPM,'ConfirmLimit',d.ConfirmCandidates, ...
            'Mode',string(state.opts.mode),'Count',n,'PlannedCount',state.opts.trials,'Completed',all(state.done(kk,:)));
        if n>0
            f=q.F0Error_Hz; a=q.RateError_HzPerS; track=q.TrackRMSE_Hz;
            sf=sum(q.F0Success); sa=sum(q.RateSuccess); sj=sum(q.JointSuccess);
            r.InvalidCount=sum(~q.Valid); r.NonfiniteCount=sum(~isfinite(f)); r.F0FailureCount=n-sf;
            r.RateFailureCount=n-sa; r.JointFailureCount=n-sj;
            r.BoundaryCount=sum(q.BoundaryHit); r.ConfirmedNonfiniteCount=sum(q.ConfirmedNonfiniteCount);
        else
            f=[]; a=[]; track=[]; sf=0; sa=0; sj=0;
            r.InvalidCount=0; r.NonfiniteCount=0; r.F0FailureCount=0; r.RateFailureCount=0;
            r.JointFailureCount=0; r.BoundaryCount=0; r.ConfirmedNonfiniteCount=0;
        end
        r.F0SuccessRate=safe_ratio(sf,n); [r.F0Lo,r.F0Hi]=wilson(sf,n);
        r.RateSuccessRate=safe_ratio(sa,n); [r.RateLo,r.RateHi]=wilson(sa,n);
        r.JointSuccessRate=safe_ratio(sj,n); [r.JointLo,r.JointHi]=wilson(sj,n);
        [r.F0RMSE_Hz,r.F0Bias_Hz,r.F0Std_Hz]=error_stats(f);
        [r.RateRMSE_HzPerS,r.RateBias_HzPerS,r.RateStd_HzPerS]=error_stats(a);
        r.AbsF0Median_Hz=absolute_quantile(f,.5); r.AbsF0P90_Hz=absolute_quantile(f,.9);
        r.AbsF0P95_Hz=absolute_quantile(f,.95);
        [r.TrackRMSE_Hz,~,~]=error_stats(track);
        r.TrackMaxError_Hz=NaN; r.MeanCoherence=NaN; r.CoherenceInvalidCount=0;
        r.MeanSearch_s=NaN; r.MeanGeneration_s=NaN; r.MeanConfirmation_s=NaN;
        r.MeanPreprocessAllocated_s=NaN; r.MeanCompute_s=NaN;
        if n>0
            r.TrackMaxError_Hz=max(q.TrackMaxError_Hz);
            r.CoherenceInvalidCount=sum(~q.CoherenceValid);
            % 不静默剔除一致性无效值；均值保持NaN，诊断计数另列。
            r.MeanCoherence=mean(q.FrameCoherence);
            r.MeanSearch_s=mean(q.Search_s); r.MeanGeneration_s=mean(q.Generation_s);
            r.MeanConfirmation_s=mean(q.Confirmation_s);
            r.MeanPreprocessAllocated_s=mean(q.PreprocessAllocated_s); r.MeanCompute_s=mean(q.Compute_s);
        end
        r.QuantileDefinition=string('median standard; P90/P95 nearest rank; nonfinite absolute errors Inf');
        r.ContinuousRefinementCalls=0; rows{kk}=r;
    end
    S=struct2table(vertcat(rows{:}));
end
function y=safe_ratio(a,n)
    if n==0; y=NaN; else; y=a/n; end
end
function y=absolute_quantile(errors,p)
    % 全样本最近秩分位数Qp=x_(ceil(p*N))。非有限误差记为Inf，不删去失败。
    if isempty(errors); y=NaN; return; end
    values=abs(errors); values(~isfinite(values))=Inf; values=sort(values);
    if p==.5; y=median(values); else; y=values(max(1,ceil(p*numel(values)))); end
end

%% 十、少量无噪声点与首批优化核对；不要求离散网格零误差
function [T,benchmark]=noiseless_checks(c,o,x)
    fixtures=[100173.7 700;90625 700;90624.63 700];
    if o.rateComparison; fixtures(:,2)=[650;700;750]; end
    base=generate_base(c,x,o.seed); rows={}; C=make_conditions(c,o);
    durations=unique([C.Duration],'stable'); benchmark=struct('Kind','noiseless_reference_check','Verified',false);
    for duration=durations
        for j=1:size(fixtures,1)
            d=C(1); d.Duration=duration; d.RateTruth=fixtures(j,2); d.EbN0=Inf;
            d.ScfFs=c.scfFs; d.ConfirmCandidates=c.confirmCandidates; d.PNPhaseChip=0; d.LocalPNRatePPM=0;
            base.fd=fixtures(j,1)-c.fc;
            clean=base.spread.*cos(2*pi*(c.fc+base.fd)*x.t+pi*d.RateTruth*x.t.^2+base.phase);
            long=resample((clean.*x.pn)',1,31,c.filterN,c.filterBeta);
            n=round(duration*d.ScfFs); first=floor((numel(long)-n)/2)+1; y=long(first:first+n-1);
            r=receiver_config(c,d); clock=tic; estimate=stage1_receiver(y,r,true); elapsed=toc(clock);
            if j==1 && duration==durations(1) && o.verifyOptimization
                clock=tic; reference=stage1_reference(y,r); baselineElapsed=toc(clock);
                compare_estimates(estimate,reference);
                benchmark=struct('Kind','noiseless_reference_check','Verified',true, ...
                    'ReferenceSearch_s',baselineElapsed,'CachedSearch_s',elapsed);
            end
            localRef=(first-1)/d.ScfFs+(n-1)/(2*d.ScfFs);
            m=evaluate_record(estimate,base,x,d,c,0,localRef,n);
            m.Search_s=elapsed; m.NoiseVariance=0; m.Fixture=j;
            m.RateDistanceToGrid_HzPerS=min(abs(d.RateTruth-(c.rateRange(1):c.rateStep:c.rateRange(2))));
            rows{end+1,1}=m; %#ok<AGROW>
        end
    end
    T=struct2table(vertcat(rows{:}));
    disp(T(:,{'Fixture','Duration_ms','TrueRate_HzPerS','F0Error_Hz','RateError_HzPerS', ...
        'RateDistanceToGrid_HzPerS','F0Success','BoundaryHit'}));
    fprintf('无噪声点只报告离散网格误差/失败，不强制断言零误差或全部成功。\n');
end
function compare_estimates(a,b)
    assert(isequal(size(a.generated),size(b.generated)) && isequal(size(a.candidates),size(b.candidates)), ...
        '优化前后候选数量不一致。');
    assert_close(a.generated(:,1:2),b.generated(:,1:2),1e-12,1e-8);
    assert_close(a.generated(:,3),b.generated(:,3),1e-8,1e-25);
    assert(isequaln(a.generated(:,4),b.generated(:,4)),'分支编号不一致。');
    assert_close(a.candidates(:,1:2),b.candidates(:,1:2),1e-12,1e-8);
    assert_close(a.candidates(:,3),b.candidates(:,3),1e-8,1e-25);
    assert_close(a.candidates(:,5),b.candidates(:,5),1e-8,1e-25);
    assert(isequaln(a.candidates(:,[4 6 7 8]),b.candidates(:,[4 6 7 8])),'候选排名或确认状态不一致。');
    assert(a.valid==b.valid && a.selectedCandidate==b.selectedCandidate,'优化前后候选选择不一致。');
    assert_close([a.fLocal a.rate],[b.fLocal b.rate],1e-12,1e-8);
    assert_close(a.score,b.score,1e-8,1e-25);
    assert_close(a.frameSequence,b.frameSequence,1e-8,1e-18);
end
function assert_close(a,b,rtol,atol)
    assert(isequal(isnan(a),isnan(b)) && isequal(isinf(a),isinf(b)),'非有限状态不一致。');
    use=isfinite(a)&isfinite(b);
    assert(all(abs(a(use)-b(use))<=atol+rtol*max(abs(a(use)),abs(b(use)))),'优化前后数值不一致，请调查。');
end
function show_timing_estimate(state)
    T=trial_table(state); if isempty(T); return; end
    components=[mean(T.PreprocessAllocated_s),mean(T.Generation_s),mean(T.Confirmation_s)];
    names={'生成/预处理分摊','第一层候选生成','精确确认'}; [~,slow]=max(components);
    totalEstimate=mean(T.Compute_s)*numel(state.done);
    remaining=mean(T.Compute_s)*(numel(state.done)-sum(state.done(:)));
    fprintf('本机首批实测：主耗时=%s；平均计算%.3fs/估计。\n',names{slow},mean(T.Compute_s));
    fprintf('粗略估算本配置总计算%.1fs，剩余%.1fs；不含核对、无噪声、保存/绘图，非耗时承诺。\n',totalEstimate,remaining);
end

%% 十一、小型数值保存与单幅汇总图；不保存高采样率波形或FIG
function atomic_save(state,path)
    folder=fileparts(path); temp=[tempname(folder) '.mat'];
    cleanup=onCleanup(@() delete_if_exists(temp)); %#ok<NASGU>
    save(temp,'state','-v7.3'); [ok,msg]=movefile(temp,path,'f'); assert(ok,'%s',msg);
end
function delete_if_exists(path)
    if exist(path,'file'); delete(path); end
end
function report=export_tables(state,folder)
    names={'trials.csv','summary.csv','noiseless.csv'};
    values={state.trials,state.summary,state.noiseless}; report=struct('success',0,'failure',0,'messages',{{}});
    for j=1:numel(names)
        temp=[tempname(folder) '.csv'];
        try
            writetable(values{j},temp); [ok,msg]=movefile(temp,fullfile(folder,names{j}),'f'); assert(ok,'%s',msg);
            report.success=report.success+1;
        catch ME
            delete_if_exists(temp); report.failure=report.failure+1; report.messages{end+1}=ME.message;
            warning('FSM:CSV','CSV失败但数值MAT已保留：%s',ME.message);
        end
    end
end
function report=export_summary_plot(state,folder)
    report=struct('success',0,'failure',0,'message',''); h=[];
    try
        h=figure('Color','w','Name','FSM stage1 fixed-window baseline'); S=state.summary;
        % 非默认附加条件仅显示理想PN/800kHz/确认10，完整数据仍在CSV/MAT。
        use=S.PNPhaseChip==0 & S.LocalPNRatePPM==0 & S.ReceiverFs_Hz==state.cfg.scfFs ...
            & S.ConfirmLimit==state.cfg.confirmCandidates;
        Q=S(use,:); rates=unique(Q.TrueRate_HzPerS,'stable'); durations=unique(Q.Duration_ms,'stable');
        subplot(1,2,1); hold on;
        for a=rates(:)'
            for t=durations(:)'
                q=sortrows(Q(Q.TrueRate_HzPerS==a & Q.Duration_ms==t,:),'EbN0_dB');
                label=sprintf('%gms, a=%g',t,a);
                errorbar(q.EbN0_dB,q.F0SuccessRate,q.F0SuccessRate-q.F0Lo,q.F0Hi-q.F0SuccessRate,'o-', ...
                    'DisplayName',label);
                for j=1:height(q)
                    if isfinite(q.F0SuccessRate(j)); text(q.EbN0_dB(j),q.F0SuccessRate(j),sprintf(' N=%d',q.Count(j)),'FontSize',7); end
                end
            end
        end
        xlabel('信息比特Eb/N0(dB)'); ylabel('f0成功率及Wilson95%'); ylim([0 1.1]); grid on; legend('show','Location','best');
        title('第一层离散估计，20Hz门限；初步统计');
        subplot(1,2,2); hold on;
        for a=rates(:)'
            for t=durations(:)'
                q=sortrows(Q(Q.TrueRate_HzPerS==a & Q.Duration_ms==t,:),'EbN0_dB');
                semilogy(q.EbN0_dB,q.F0RMSE_Hz,'o-','DisplayName',sprintf('%gms, a=%g',t,a));
            end
        end
        set(gca,'YScale','log'); xlabel('信息比特Eb/N0(dB)'); ylabel('完整样本f0 RMSE(Hz)'); grid on; legend('show','Location','best');
        if ~isgraphics(h,'figure'); error('FSM:ClosedFigure','图窗已关闭，数值结果不受影响。'); end
        number=get(h,'Number'); filename=fullfile(folder,'summary.png'); exported=false;
        if exist('exportgraphics','file')~=0 || exist('exportgraphics','builtin')~=0
            try; exportgraphics(h,filename,'Resolution',150); exported=true;
            catch ME; report.message=ME.message; end
        end
        if ~exported
            assert(isgraphics(h,'figure'),'图窗已关闭。'); print(number,filename,'-dpng','-r150');
        end
        report.success=1;
    catch ME
        report.failure=1; report.message=ME.message;
        warning('FSM:PNG','单图导出失败，数值MAT已保存：%s',ME.message);
    end
end

%% 十二、无需Monte Carlo的原生自检入口（运行状态记录在state.checks）
function report=internal_checks(c)
    assert(isequal([c.Rb c.Rc c.fs c.fc c.scfFs c.L c.hop],[1e5 31e5 24.8e6 1e5 8e5 128 64]));
    assert(c.fTolerance==20 && c.rateTolerance==200 && c.rateStep==100 && c.branchPeaks==3 && c.confirmCandidates==10);
    o=struct('durations',.1,'rates',700,'ebN0',[-10 -7 -4], ...
        'windowComparison',false,'rateComparison',false,'syncComparison',false, ...
        'bandwidthComparison',false,'candidateComparison',false);
    C=make_conditions(c,o); assert(numel(C)==3);
    o.rateComparison=true; C2=make_conditions(c,o); assert(numel(C2)==9);
    assert(sum([C2.RateTruth]==700)==3); % 700条件不重复产生
    o.windowComparison=true; C3=make_conditions(c,o); assert(numel(C3)==36);
    receiverCfg=receiver_config(c,C(1)); assert(~isfield(receiverCfg,'RateTruth') && ~isfield(receiverCfg,'fd'));
    indices=distinct_indices([100000 0;100001 5;100010 0;100000 300],5,200);
    assert(isequal(indices,[1;3;4]));
    [r,b,s]=error_stats([1 NaN]); assert(isinf(r)&&isnan(b)&&isnan(s));
    assert(absolute_quantile([1 2 Inf],.95)==Inf);
    [lo,hi]=wilson(50,50); assert(lo<1 && abs(hi-1)<1e-12);
    common=(2480000-1)/(2*c.fs);
    for duration=[.001 .01 .05 .1]
        n=round(duration*c.scfFs); first=floor((80000-n)/2)+1;
        local=(first-1)/c.scfFs+(n-1)/(2*c.scfFs);
        f0=123456.7; a=700; localF=f0+a*(local-common);
        assert(abs(localF+a*(common-local)-f0)<1e-9);
        if duration==.1; assert(abs(a*(n-1)/(2*c.scfFs)-35)<.001); end
    end
    % 纯接收数据小例子：缓存与不缓存实现必须给出同一第一层结果。
    receiverCfg.confirmCandidates=3;
    y=cos(2*pi*100173.7*(0:799)'/c.scfFs);
    a=stage1_receiver(y,receiverCfg,true); b=stage1_reference(y,receiverCfg);
    compare_estimates(a,b);
    report=struct('Passed',true,'ContinuousRefinementCalls',0,'DefaultConditions',3, ...
        'RateComparisonNoDuplicate700',true,'CommonReference',true,'CandidateCacheEquality',true, ...
        'FailureStatistics',true,'Scope','Native only when this MATLAB entry actually runs');
end
