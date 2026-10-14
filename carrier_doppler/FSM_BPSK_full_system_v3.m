%% FSM频域循环谱：完整实通带DSSS-BPSK系统仿真（v3）
% 链路：随机比特 -> Gold码扩频 -> 矩形码片采样 -> 实通带BPSK -> 信道 -> 接收机。
% 接收核心：频移FFT -> Xp.*conj(Xm) -> 复数跨帧平均 -> 频率平滑 -> 最后取模。
% 不是信号平方估计器；默认以已知、准确同步PN为前提，不是盲码捕获或调制分类。
% MATLAB R2018b+；需要Signal Processing Toolbox的resample/freqz。
%
% 直接运行：默认quick，小规模验证，不自动运行数小时的full。
% 正式实验：fsm_options=struct('mode','full'); run('FSM_BPSK_full_system_v3.m');
% 恢复实验：fsm_options=struct('resumeFile','结果目录/checkpoint.mat'); run(...);
% 保存/候选/回退自检：fsm_options=struct('mode','selftest'); run(...);
% 原FSM_BPSK_full_system.m保持原样；本文件所有局部函数均在文件末尾。
clearvars -except fsm_options; clc;
if ~exist('fsm_options','var'); fsm_options=struct(); end
version='FSM_v3_paired_20261014_1';
sourceText=fileread([mfilename('fullpath') '.m']);

%% 一、参数：保留原信号、正式观测时间、先验、误差门限、1000次/条件
cfg.Rb=1e5; cfg.Rc=31e5; cfg.fs=24.8e6; cfg.fc=1e5;
cfg.T=.1; cfg.longT=.1; cfg.scfFs=8e5; cfg.L=128; cfg.hop=64;
cfg.smoothingBins=17; cfg.fRange=[50e3 150e3];
cfg.rateRange=[-1000 1000]; cfg.rateStep=100;
cfg.fTolerance=20; cfg.rateTolerance=200;
cfg.fdTest=45000; cfg.rateTest=900;
cfg.useKnownPN=true; cfg.rxPhaseChip=0; cfg.rxRatePPM=0; cfg.jerk=0;
cfg.filterN=10; cfg.filterBeta=5; % 显式Kaiser重采样滤波器，而非仅修改带宽标题
cfg.branchPeaks=3;              % 每个FFT分支保留3个二维分离候选
cfg.confirmCandidates=10;       % 精确评分确认数，与下行多起点数不是一个量
cfg.multiStarts=3;              % 默认从精确评分最高的3个不同候选启动
cfg.dedupFrequencyFactor=.5;    % dF=0.5/T Hz，边缘最大额外循环相位约pi
cfg.dedupRateFactor=2;          % dA=2/T^2 Hz/s，边缘额外循环相位约pi
cfg.maxRefineIter=200;
opts=struct('mode','quick','quickTrials',3,'seed',20261014,'makePlots',true, ...
    'ioSelfTest',true,'checkpointEveryTrials',50,'maxNewConditions',Inf);
fields=fieldnames(fsm_options);
for z=1:numel(fields); opts.(fields{z})=fsm_options.(fields{z}); end
assert(any(strcmp(opts.mode,{'quick','full','selftest'})),'Unknown mode.');
assert(opts.quickTrials>=1 && opts.quickTrials==floor(opts.quickTrials));
assert(opts.checkpointEveryTrials>=1 && opts.checkpointEveryTrials==floor(opts.checkpointEveryTrials));
assert(cfg.branchPeaks>=2 && cfg.multiStarts>=1 && cfg.multiStarts==floor(cfg.multiStarts));
assert(exist('resample','file')~=0 && exist('freqz','file')~=0,'需要Signal Processing Toolbox。');
assert(cfg.fs/cfg.Rc==round(cfg.fs/cfg.Rc));
assert(cfg.fc+cfg.fdTest+cfg.Rb<cfg.scfFs/2,'接收带宽不能覆盖第一主瓣。');
assert(mod(cfg.smoothingBins,2)==1 && cfg.smoothingBins<cfg.L);

%% 二、状态恢复：完成标记按条件和试验编号保存，禁止追加重复计数
if isfield(opts,'resumeFile') && ~isempty(opts.resumeFile)
    loaded=load(opts.resumeFile,'state'); state=loaded.state;
    assert(strcmp(state.version,version),'检查点版本不匹配。');
    assert(strcmp(state.sourceText,sourceText),'源文件改变；请使用原版本恢复或另建实验。');
    outputDir=fileparts(opts.resumeFile); if isempty(outputDir); outputDir=pwd; end
    cfg=state.cfg; savedOpts=state.opts;
    operational={'makePlots','ioSelfTest','checkpointEveryTrials'};
    for ii=1:numel(operational)
        if isfield(fsm_options,operational{ii}); savedOpts.(operational{ii})=opts.(operational{ii}); end
    end
    savedOpts.maxNewConditions=opts.maxNewConditions; opts=savedOpts; state.opts=opts;
    rng(state.rngState); % 同时保存/恢复全局状态；基础记录另按固定seed重建
    fprintf('恢复：%d/%d个条件完成。\n',sum(state.completed),numel(state.conditions));
else
    rng(opts.seed,'twister');
    outputDir=unique_directory(fullfile(pwd,['FSM_v3_' datestr(now,'yyyymmdd_HHMMSS')]));
    mkdir(outputDir);
    if strcmp(opts.mode,'full'); nTrials=1000; else; nTrials=opts.quickTrials; end
    conditions=make_conditions(cfg,nTrials,opts.mode);
    state=struct('version',version,'sourceText',sourceText,'cfg',cfg,'opts',opts, ...
        'conditions',conditions,'completed',false(numel(conditions),1), ...
        'records',{cell(numel(conditions),1)},'done',{cell(numel(conditions),1)}, ...
        'rngState',rng,'lastBaseRngState',[],'phase','initialized','lastException','');
    state.baseSeeds=opts.seed+(1:max(nTrials,1));
    for k=1:numel(conditions)
        state.records{k}=cell(conditions(k).nTrials,1);
        state.done{k}=false(conditions(k).nTrials,1);
    end
    atomic_state_save(state,fullfile(outputDir,'checkpoint.mat'));
end
fprintf('输出目录：%s\n',outputDir);
fprintf('模式%s；正式模式保持1000次/条件，当前不自动升级为full。\n',opts.mode);
fprintf('-7 dB信息比特Eb/N0 -> 原实采样SNR %.3f dB。\n',-7+10*log10(2*cfg.Rb/cfg.fs));

%% 三、内置小规模自检：不把自检或Python结果冒充正式实验
checks=internal_checks(cfg);
state.internalChecks=checks;
if opts.ioSelfTest
    state.ioChecks=io_selftest(outputDir);
end
atomic_state_save(state,fullfile(outputDir,'checkpoint.mat'));
if strcmp(opts.mode,'selftest')
    state.phase='selftest_complete';
    atomic_state_save(state,fullfile(outputDir,'results.mat'));
    fprintf('仅运行自检；未执行Monte Carlo。\n');
    return;
end
code=make_gold(); newConditions=0;

%% 四、配对实验：信噪比/时长/敏感性/候选数/真实带宽对照
% 条件循环便于逐条件保存。每个TrialID用同一个独立RandStream重建100ms长记录。
% 记录相同基础比特、初相、fD、a和标准高斯噪声。不同Eb/N0仅重标噪声幅度。
% duration先滤完整长记录，再截取同一中心段；不会重新生成短记录。
for k=1:numel(state.conditions)
    if state.completed(k); continue; end
    if newConditions>=opts.maxNewConditions; break; end
    currentCondition=state.conditions(k);
    fprintf('\n[%d/%d] %s (%d次)\n',k,numel(state.conditions),currentCondition.name,currentCondition.nTrials);
    if currentCondition.reuse>0
        assert(state.completed(currentCondition.reuse),'复用条件尚未完成。');
        state.records{k}=state.records{currentCondition.reuse};
        state.done{k}=state.done{currentCondition.reuse};
    else
        for trial=1:currentCondition.nTrials
            if state.done{k}(trial); continue; end
            seed=state.baseSeeds(trial);
            try
                base=make_base(cfg,code,seed);
                [record,demo]=run_trial(base,currentCondition,code,trial);
                state.records{k}{trial}=record;
                state.done{k}(trial)=true; % 仅在整条记录完整后提交这一试验编号
                state.lastBaseRngState=base.rngState;
                if strcmp(currentCondition.name,'EbN0_minus7') && trial==1; state.demo=demo; end
            catch ME
                % 编程/资源异常不是科学失效样本：不计为完成，下次重试同seed、同编号。
                % 数值无效估计则由run_trial正常返回并计为失败，不进入这里。
                state.lastException=sprintf('%s / trial %d: %s',currentCondition.name,trial,getReport(ME));
                state.rngState=rng;
                atomic_state_save(state,fullfile(outputDir,'checkpoint.mat'));
                rethrow(ME);
            end
            if mod(trial,opts.checkpointEveryTrials)==0
                state.rngState=rng;
                atomic_state_save(state,fullfile(outputDir,'checkpoint.mat'));
            end
            if mod(trial,25)==0 || trial==currentCondition.nTrials
                fprintf('  完成%d/%d\n',sum(state.done{k}),currentCondition.nTrials);
            end
        end
    end
    assert(all(state.done{k}));
    state.rngState=rng; state.phase='condition_trials_done';
    % 先提交逐次完成标记；若条件MAT导出失败，恢复时不重算试验，只重试保存。
    atomic_state_save(state,fullfile(outputDir,'checkpoint.mat'));
    bundle=struct('condition',currentCondition,'records',{state.records{k}}, ...
        'baseSeeds',state.baseSeeds,'rngState',state.rngState);
    atomic_bundle_save(bundle,fullfile(outputDir,[currentCondition.name '_numeric.mat']));
    state.completed(k)=true; state.phase='conditions_running';
    atomic_state_save(state,fullfile(outputDir,'checkpoint.mat'));
    newConditions=newConditions+1;
end

%% 五、先保存最终/当前数值结果，再进行任何正式绘图或导图
state.summary=make_summary(state);
state.groups=make_groups(state,cfg);
state.rngState=rng; state.numericComplete=all(state.completed);
state.phase='numerics_saved_before_plot';
atomic_state_save(state,fullfile(outputDir,'results.mat'));
atomic_state_save(state,fullfile(outputDir,'checkpoint.mat'));
state.numericExport=export_numeric_tables(state,outputDir);
atomic_state_save(state,fullfile(outputDir,'results.mat'));
disp(state.summary(:,{'Condition','Count','FinalJoint','FinalLo','FinalHi', ...
    'FinalFailureCount','InvalidFinalCount','FallbackCount','TrackRMSE_Hz'}));
fprintf('MAT数值结果已保存；CSV成功%d，失败%d。\n', ...
    state.numericExport.success,state.numericExport.failure);
if ~state.numericComplete
    fprintf('本次按maxNewConditions暂停；已完成%d/%d条件。用checkpoint.mat恢复，不重复计数。\n', ...
        sum(state.completed),numel(state.completed));
    return;
end

%% 六、绘图和图片导出失败不能破坏已经保存的数值结果
if opts.makePlots
    [figures,creationErrors]=create_figures(state);
    state.plotCreationErrors=creationErrors;
    state.imageExport=export_figures(figures,outputDir,false);
else
    state.imageExport=struct('pngSuccess',0,'pngFailure',0,'figSuccess',0,'figFailure',0,'details',[]);
end
state.phase='complete';
try
    atomic_state_save(state,fullfile(outputDir,'results.mat'));
    atomic_state_save(state,fullfile(outputDir,'checkpoint.mat'));
catch ME
    warning('FSM:FinalReport','导图报告保存失败；导图前数值MAT仍保留：%s',ME.message);
end
fprintf('PNG：成功%d/失败%d；FIG：成功%d/失败%d。\n',state.imageExport.pngSuccess, ...
    state.imageExport.pngFailure,state.imageExport.figSuccess,state.imageExport.figFailure);
fprintf('读取：S=load(''%s''); T=S.state.summary;\n',fullfile(outputDir,'results.mat'));

%% 七、条件表：保持原科学参数；quick缩小次数/无噪声扫描，不改正式定义
function plan=make_conditions(c,n,mode)
    plan=struct('name',{},'stage',{},'cfg',{},'eb',{},'nTrials',{},'reuse',{},'scanF',{},'scanA',{});
    df=c.scfFs/c.L;
    boundaries=(ceil(c.fRange(1)/df-.5):floor(c.fRange(2)/df-.5))*df+df/2;
    boundaries=boundaries(boundaries>c.fRange(1)&boundaries<c.fRange(2));
    if strcmp(mode,'full')
        f=unique([c.fRange 50037.7 100173.7 149963.3 boundaries-.37 boundaries boundaries+.37]);
        a=[-1000 -873.1 0 731.4 1000];
    else
        boundary=boundaries(ceil(numel(boundaries)/2));
        f=unique([c.fRange 50037.7 100173.7 boundary-.37 boundary boundary+.37]);
        a=[-1000 413.9 1000];
    end
    for i=1:numel(f)
        for j=1:numel(a)
            plan(end+1)=make_condition_entry(sprintf('Scan_%03d_%02d',i,j),'scan',c,Inf,1,0,f(i),a(j)); %#ok<AGROW>
        end
    end
    eb=[-10 -8 -7 -6 -4]; baseline=0;
    for x=eb
        name=sprintf('EbN0_minus%d',abs(x));
        plan(end+1)=make_condition_entry(name,'snr',c,x,n,0,NaN,NaN); %#ok<AGROW>
        if x==-7; baseline=numel(plan); end
    end
    for count=[5 10 Inf]
        d=c; d.confirmCandidates=count;
        reuse=0; if count==c.confirmCandidates; reuse=baseline; end
        if isinf(count); name='Candidates_all_distinct'; else; name=sprintf('Candidates_%d',count); end
        plan(end+1)=make_condition_entry(name,'candidates',d,-7,n,reuse,NaN,NaN); %#ok<AGROW>
    end
    for duration=[.001 .01 .05 .1]
        d=c; d.T=duration; reuse=0; if duration==c.T; reuse=baseline; end
        plan(end+1)=make_condition_entry(sprintf('Duration_%g_ms',duration*1000),'duration',d,-7,n,reuse,NaN,NaN); %#ok<AGROW>
    end
    names={'IdealSync','Phase_0p125_chip','Phase_0p5_chip','Phase_1_chip', ...
        'Rate_1_ppm','Rate_5_ppm','Rate_25_ppm','Jerk_1000','Jerk_5000','RawFSM_800k_BandlimitedNoPN'};
    phases=[0 .125 .5 1 0 0 0 0 0 0]; ppms=[0 0 0 0 1 5 25 0 0 0];
    jerks=[0 0 0 0 0 0 0 1000 5000 0];
    for i=1:numel(names)
        d=c; d.rxPhaseChip=phases(i); d.rxRatePPM=ppms(i); d.jerk=jerks(i);
        if i==numel(names); d.useKnownPN=false; end
        reuse=0; if i==1; reuse=baseline; end
        plan(end+1)=make_condition_entry(names{i},'sensitivity',d,-7,n,reuse,NaN,NaN); %#ok<AGROW>
    end
    for fsout=[8e5 16e5]
        d=c; d.scfFs=fsout; d.L=round(c.L*fsout/c.scfFs); d.hop=d.L/2;
        reuse=0; if fsout==c.scfFs; reuse=baseline; end
        plan(end+1)=make_condition_entry(sprintf('Bandwidth_%g_kSps',fsout/1000),'bandwidth',d,-7,n,reuse,NaN,NaN); %#ok<AGROW>
    end
end
function d=make_condition_entry(name,stage,c,eb,n,reuse,f,a)
    d=struct('name',name,'stage',stage,'cfg',c,'eb',eb,'nTrials',n,'reuse',reuse,'scanF',f,'scanA',a);
end

%% 八、可重建的基础长记录：seed与TrialID绑定，不依赖条件运行顺序
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

%% 九、发射/信道/前端与评价分离：真值不进入fsm_receiver
function [record,demo]=run_trial(b,d,code,id)
    c=d.cfg; fd=b.fd; a=b.rate;
    if isfinite(d.scanF); fd=d.scanF-c.fc; a=d.scanA; end
    spread=b.data.*b.pn;
    channelPhase=2*pi*fd*b.t+pi*a*b.t.^2+pi/3*c.jerk*b.t.^3+b.phase;
    clean=spread.*cos(2*pi*c.fc*b.t+channelPhase);
    variance=0;
    if ~isinf(d.eb); variance=mean(clean.^2)*c.fs/(2*c.Rb*10^(d.eb/10)); end
    % Eb/N0定义保持实AWGN双边N0/2约定。敏感性波形功率若改变，仅重标同一标准噪声。
    rx=clean+sqrt(variance)*b.noise;
    if c.useKnownPN
        Ns=c.fs/c.Rc;
        ci=floor((0:numel(rx)-1)*(1+c.rxRatePPM*1e-6)/Ns+c.rxPhaseChip);
        processed=rx.*code(mod(ci,31)+1);
    else; processed=rx; end
    [up,down]=rat(c.scfFs/c.fs,1e-12);
    [yLong,filterB]=resample(processed(:),up,down,c.filterN,c.filterBeta);
    % 先滤整条记录再截中心段；不同时长共用同样的滤波后噪声及数据。
    n=round(c.T*c.scfFs); first=floor((numel(yLong)-n)/2)+1;
    assert(first>=1 && first+n-1<=numel(yLong));
    y=yLong(first:first+n-1); startTime=(first-1)/c.scfFs;
    receiverCfg=receiver_config(c); % 只传搜索先验/FFT参数，不传fd/a/jerk/bits
    search=fsm_receiver(y,receiverCfg);
    reference=startTime+search.reference;
    delta=reference-b.reference;
    trueF=c.fc+fd+a*delta+.5*c.jerk*delta^2; trueA=a+c.jerk*delta;
    m=metric_template(); m.TrialID=id; m.BaseSeed=b.seed; m.ReferenceTime_s=reference;
    m.TrueF_Hz=trueF; m.TrueRate_HzPerS=trueA;
    m.CoarseF_Hz=search.fCoarse; m.CoarseRate_HzPerS=search.aCoarse;
    m.FineF_Hz=search.fFine; m.FineRate_HzPerS=search.aFine;
    m.CoarseF_Error=search.fCoarse-trueF; m.CoarseRate_Error=search.aCoarse-trueA;
    m.FineF_Error=search.fFine-trueF; m.FineRate_Error=search.aFine-trueA;
    [m.CoarseTrackRMSE_Hz,m.CoarseTrackMaxError_Hz]=track_error( ...
        m.CoarseF_Error,m.CoarseRate_Error,c.jerk,n,c.scfFs);
    [m.TrackRMSE_Hz,m.TrackMaxError_Hz]=track_error(m.FineF_Error,m.FineRate_Error,c.jerk,n,c.scfFs);
    m.CoarseValid=search.coarseValid; m.FinalValid=search.finalValid;
    m.CoarseSuccess=m.CoarseValid && abs(m.CoarseF_Error)<=c.fTolerance && abs(m.CoarseRate_Error)<=c.rateTolerance;
    m.FinalSuccess=m.FinalValid && abs(m.FineF_Error)<=c.fTolerance && abs(m.FineRate_Error)<=c.rateTolerance;
    m.OptimizerAnyConverged=any(search.refinements(:,6)>0);
    m.OptimizerValidCount=sum(search.refinements(:,7)>0);
    m.OptimizerAttemptCount=size(search.refinements,1);
    m.OptimizerNonConvergedCount=sum(search.refinements(:,6)<=0);
    m.OptimizerInvalidCount=sum(search.refinements(:,7)==0);
    m.ConfirmedNonfiniteScoreCount=sum(search.candidates(:,8)>0 & ~isfinite(search.candidates(:,5)));
    m.SelectedOptimizerConverged=search.selectedExitFlag>0 && ~search.fallback;
    m.OptimizerExitFlag=search.selectedExitFlag; m.RefineFallback=search.fallback;
    m.BoundaryHit=search.boundary; m.CoarseScore=search.coarseScore; m.FinalScore=search.finalScore;
    m.FrameCount=search.frameCount; m.GeneratedCandidates=size(search.generated,1);
    m.DistinctCandidates=size(search.candidates,1); m.ConfirmedCandidates=sum(search.candidates(:,8)>0);
    m.MultiStartsUsed=size(search.refinements,1); m.NoiseVariance=variance;
    m.Duration_s=c.T; m.ScfFs_Hz=c.scfFs;
    filterInfo=filter_metadata(filterB,c,up,down);
    record=struct('metrics',m,'search',search,'filter',filterInfo);
    demo=struct();
    if strcmp(d.name,'EbN0_minus7') && id==1
        ns=min(numel(rx),round(60e-6*c.fs));
        demo=struct('t_us',b.rawTime(1:ns)*1e6,'data',b.data(1:ns),'spread',spread(1:ns), ...
            'tx',spread(1:ns).*cos(2*pi*c.fc*b.t(1:ns)), ...
            'clean',clean(1:ns),'rx',rx(1:ns),'fd',fd,'rate',a,'phase',b.phase);
        demo.spectrum=display_spectrum(y,receiverCfg,search);
    end
end
function r=receiver_config(c)
    % 明确的接口白名单：任何真实发送参数都不在此结构中。
    names={'scfFs','L','hop','smoothingBins','fRange','rateRange','rateStep', ...
        'branchPeaks','confirmCandidates','multiStarts','dedupFrequencyFactor', ...
        'dedupRateFactor','maxRefineIter'};
    r=struct(); for i=1:numel(names); r.(names{i})=c.(names{i}); end
end

%% 十、FSM粗搜索：每分支多候选、基于观测时长的二维去重
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

    %% 十一、多起点精搜索：同一精确FSM指标，不以真值或误差筛选
    nStarts=min(c.multiStarts,numel(ranking)); R=NaN(nStarts,8); messages=cell(nStarts,1);
    normScore=max(e.coarseScore,realmin);
    for j=1:nStarts
        start=C(ranking(j),1:2);
        objective=@(p) refine_objective(p,frames,tc,start,T,c,w,band,normScore);
        flag=-999; point=[NaN NaN]; score=NaN;
        try
            [p,~,flag]=fminsearch(objective,[0 0],optimset('Display','off', ...
                'TolX',1e-5,'TolFun',1e-8,'MaxIter',c.maxRefineIter));
            point=start+[p(1)/T p(2)/T^2];
            if in_prior(point,c); score=fsm_score(frames,tc,point(1),point(2),c,w,band); end
        catch ME
            messages{j}=ME.message;
        end
        outputValid=in_prior(point,c) && isfinite(score) && score>=0;
        % ExitFlag<=0与数值无效是两回事。所有有效输出均按精确指标比较；
        % 未收敛但有效的输出也保留其状态，不冒充“收敛”。数值失败则不可选。
        R(j,:)=[start point score flag double(outputValid) double(flag>0)];
    end
    e.refinements=R; e.refineMessages=messages;
    e.refinementColumns={'StartF_Hz','StartRate_HzPerS','OutputF_Hz','OutputRate_HzPerS', ...
        'ExactScore','ExitFlag','OutputValid','Converged'};
    [point,score,selected,reason]=choose_final([e.fCoarse e.aCoarse e.coarseScore],R,c);
    e.fFine=point(1); e.aFine=point(2); e.finalScore=score;
    e.fallback=selected==0; e.fallbackReason=reason; e.selectedStart=selected;
    e.selectedExitFlag=NaN; if selected>0; e.selectedExitFlag=R(selected,6); end
    e.finalValid=in_prior(point,c) && isfinite(score) && score>=0;
    if e.finalValid
        assert(e.finalScore>=e.coarseScore,'最终指标不应低于最佳保留粗候选。');
        ftol=min(dF,.01*diff(c.fRange)); atol=min(dA,.01*diff(c.rateRange));
        e.boundary=min(abs(point(1)-c.fRange))<=ftol || min(abs(point(2)-c.rateRange))<=atol;
    end
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
function [point,score,selected,reason]=choose_final(coarse,R,c)
    point=coarse(1:2); score=coarse(3); selected=0; reason='NoValidRefinement';
    if ~in_prior(point,c) || ~isfinite(score); point=[NaN NaN]; score=NaN; reason='NoValidCoarse'; end
    for j=1:size(R,1)
        if R(j,7)>0 && in_prior(R(j,3:4),c) && isfinite(R(j,5)) && R(j,5)>=0
            if ~isfinite(score) || R(j,5)>score
                point=R(j,3:4); score=R(j,5); selected=j; reason='RefinementSelected';
            elseif selected==0; reason='NoStrictScoreImprovement'; end
        end
    end
    % 这里只保证指标不差，不保证真实频率/变化率误差或轨迹误差每次更小。
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

%% 十二、真正的频域循环谱及变化率补偿：不做时域信号平方
function spec=fsm_spectrum(frames,tc,f,a,c,w)
    % alpha=2*f。相位psi=2*pi*f*t+pi*a*t^2；交叉乘积中的循环相位是2*psi。
    % a=0退化为普通移频FSM；a非零是多项式相位补偿的FSM，不声称原chirp平稳。
    phase=2*pi*f*tc+pi*a*tc.^2;
    Xp=fft(frames.*exp(-1j*phase),[],1);
    Xm=fft(frames.*exp( 1j*phase),[],1);
    spec=fftshift(mean(Xp.*conj(Xm),2))/(c.scfFs*sum(w.^2));
end
function J=fsm_score(frames,tc,f,a,c,w,band)
    spec=fsm_spectrum(frames,tc,f,a,c,w);
    % 先复数频率平滑，最后取模；17点平滑在f=0处的值作为统一精确指标。
    J=abs(mean(spec(c.L/2+1+band)))^2;
end
function v=refine_objective(p,frames,tc,start,T,c,w,band,normalizer)
    point=start+[p(1)/T p(2)/T^2];
    if ~in_prior(point,c); v=1e12; return; end
    J=fsm_score(frames,tc,point(1),point(2),c,w,band);
    if ~isfinite(J); v=1e12; else; v=-J/normalizer; end
end
function d=display_spectrum(y,c,e)
    d=struct(); if ~e.finalValid; return; end
    [frames,tc,~,w]=prepare_frames(y,c);
    f=linspace(max(c.fRange(1),e.fFine-50),min(c.fRange(2),e.fFine+50),41);
    kernel=ones(c.smoothingBins,1)/c.smoothingBins; S=zeros(numel(f),c.L);
    for j=1:numel(f)
        spec=fsm_spectrum(frames,tc,f(j),e.aFine,c,w);
        S(j,:)=conv(spec,kernel,'same').';
    end
    d=struct('fAxis',(-c.L/2:c.L/2-1)*c.scfFs/c.L,'alphaAxis',2*f, ...
        'S',S,'fineF',e.fFine,'fineRate',e.aFine);
end
function info=filter_metadata(b,c,up,down)
    % 保存实际重采样FIR，而不是凭标题声称带宽改变。归一化到直流增益1。
    [h,f]=freqz(b,1,4096,c.fs*up); gain=abs(h)/max(abs(sum(b)),realmin);
    pass=find(gain>=10^(-1/20),1,'last'); stop=find(gain<=10^(-40/20),1,'first');
    info=struct('Numerator',b,'Up',up,'Down',down,'InputFs_Hz',c.fs, ...
        'OutputFs_Hz',c.scfFs,'Nyquist_Hz',c.scfFs/2,'FFTLength',c.L, ...
        'FFTSpacing_Hz',c.scfFs/c.L,'SmoothingWidth_Hz',c.smoothingBins*c.scfFs/c.L, ...
        'KaiserN',c.filterN,'KaiserBeta',c.filterBeta,'OneDBPassEdge_Hz',NaN,'FortyDBStopEdge_Hz',NaN);
    if ~isempty(pass); info.OneDBPassEdge_Hz=f(pass); end
    if ~isempty(stop); info.FortyDBStopEdge_Hz=f(stop); end
end

%% 十三、整段轨迹评价：与频率、斜率估计采用同一个参考时刻
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
function m=metric_template()
    names={'TrialID','BaseSeed','ReferenceTime_s','TrueF_Hz','TrueRate_HzPerS', ...
        'CoarseF_Hz','CoarseRate_HzPerS','FineF_Hz','FineRate_HzPerS', ...
        'CoarseF_Error','CoarseRate_Error','FineF_Error','FineRate_Error', ...
        'CoarseTrackRMSE_Hz','CoarseTrackMaxError_Hz','TrackRMSE_Hz','TrackMaxError_Hz', ...
        'CoarseValid','FinalValid','CoarseSuccess','FinalSuccess','OptimizerAnyConverged', ...
        'OptimizerValidCount','SelectedOptimizerConverged','OptimizerExitFlag','RefineFallback', ...
        'BoundaryHit','CoarseScore','FinalScore','FrameCount','GeneratedCandidates', ...
        'DistinctCandidates','ConfirmedCandidates','MultiStartsUsed','NoiseVariance','Duration_s','ScfFs_Hz', ...
        'OptimizerAttemptCount','OptimizerNonConvergedCount','OptimizerInvalidCount','ConfirmedNonfiniteScoreCount'};
    m=cell2struct(num2cell(NaN(1,numel(names))),names,2);
end

%% 十四、统计口径：所有完成试验进入分母，无效结果明确计数，不静默omitnan
function T=metrics_table(records,done)
    use=find(done); template=metric_template();
    if isempty(use); T=struct2table(template); T(1,:)=[]; return; end
    M=repmat(template,numel(use),1);
    for j=1:numel(use); M(j)=records{use(j)}.metrics; end
    T=struct2table(M);
    assert(width(T)==numel(fieldnames(template)) && height(T)==numel(use));
    assert(isequal(T.TrialID,use(:)) && numel(unique(T.TrialID))==height(T),'试验编号重复或错位。');
end
function T=make_summary(state)
    rows=cell(numel(state.conditions),1);
    for k=1:numel(state.conditions)
        d=state.conditions(k); m=metrics_table(state.records{k},state.done{k});
        r=summarize_metrics(m);
        r.Condition=string(d.name); r.Stage=string(d.stage); r.EbN0_dB=d.eb;
        r.Duration_ms=1000*d.cfg.T; r.PlannedCount=d.nTrials; r.Completed=state.completed(k);
        r.ConfirmLimit=d.cfg.confirmCandidates; r.RequestedMultiStarts=d.cfg.multiStarts;
        r.ReceiverFs_Hz=d.cfg.scfFs; r.FFTLength=d.cfg.L; r.ReusedConditionIndex=d.reuse;
        r.OneDBPassEdge_Hz=NaN; r.FortyDBStopEdge_Hz=NaN;
        id=find(state.done{k},1);
        if ~isempty(id)
            f=state.records{k}{id}.filter;
            r.OneDBPassEdge_Hz=f.OneDBPassEdge_Hz; r.FortyDBStopEdge_Hz=f.FortyDBStopEdge_Hz;
        end
        rows{k}=r;
    end
    T=struct2table(vertcat(rows{:}));
end
function r=summarize_metrics(m)
    n=height(m); r=struct(); r.Count=n;
    r.CoarseFailureCount=n-sum(m.CoarseSuccess); r.FinalFailureCount=n-sum(m.FinalSuccess);
    r.InvalidCoarseCount=sum(~logical(m.CoarseValid)); r.InvalidFinalCount=sum(~logical(m.FinalValid));
    r.NonfiniteFinalCount=sum(~isfinite(m.FineF_Error)|~isfinite(m.FineRate_Error));
    r.FallbackCount=sum(m.RefineFallback); r.OptimizerAnyConvergedCount=sum(m.OptimizerAnyConverged);
    r.SelectedOptimizerConvergedCount=sum(m.SelectedOptimizerConverged);
    r.ValidOptimizerOutputCount=sum(m.OptimizerValidCount); r.BoundaryHitCount=sum(m.BoundaryHit);
    r.OptimizerAttemptCount=sum(m.OptimizerAttemptCount);
    r.OptimizerNonConvergedCount=sum(m.OptimizerNonConvergedCount);
    r.OptimizerInvalidOutputCount=sum(m.OptimizerInvalidCount);
    r.ConfirmedNonfiniteScoreCount=sum(m.ConfirmedNonfiniteScoreCount);
    r.CoarseJoint=ratio(sum(m.CoarseSuccess),n); r.FinalJoint=ratio(sum(m.FinalSuccess),n);
    [r.CoarseLo,r.CoarseHi]=wilson(sum(m.CoarseSuccess),n);
    [r.FinalLo,r.FinalHi]=wilson(sum(m.FinalSuccess),n);
    [r.CoarseF_RMSE,r.CoarseF_Bias,r.CoarseF_Std]=error_stats(m.CoarseF_Error);
    [r.CoarseRate_RMSE,r.CoarseRate_Bias,r.CoarseRate_Std]=error_stats(m.CoarseRate_Error);
    [r.FineF_RMSE,r.FineF_Bias,r.FineF_Std]=error_stats(m.FineF_Error);
    [r.FineRate_RMSE,r.FineRate_Bias,r.FineRate_Std]=error_stats(m.FineRate_Error);
    [r.TrackRMSE_Hz,~,~]=error_stats(m.TrackRMSE_Hz);
    [r.CoarseTrackRMSE_Hz,~,~]=error_stats(m.CoarseTrackRMSE_Hz);
    r.TrackMaxError_Hz=max_or_nan(m.TrackMaxError_Hz);
    r.CoarseTrackMaxError_Hz=max_or_nan(m.CoarseTrackMaxError_Hz);
    r.StrictFrequencyWin=ratio(sum(abs(m.FineF_Error)<abs(m.CoarseF_Error)),n);
    r.StrictRateWin=ratio(sum(abs(m.FineRate_Error)<abs(m.CoarseRate_Error)),n);
end
function [rmse,bias,sd]=error_stats(v)
    if isempty(v); rmse=NaN; bias=NaN; sd=NaN;
    elseif any(~isfinite(v)); rmse=Inf; bias=NaN; sd=NaN;
    else; rmse=sqrt(mean(v.^2)); bias=mean(v); sd=std(v,0); end
    % 一旦存在非有限误差，完整样本RMSE报Inf、偏差/标准差报NaN；另有失效计数。
    % 不在删去失败样本后报告一个看似更好的有限RMSE。
end
function y=max_or_nan(v)
    if isempty(v); y=NaN; elseif any(~isfinite(v)); y=Inf; else; y=max(v); end
end
function y=ratio(a,b)
    if b==0; y=NaN; else; y=a/b; end
end
function [lo,hi]=wilson(success,n)
    if n==0; lo=NaN; hi=NaN; return; end
    z=1.96; p=success/n; den=1+z^2/n;
    center=(p+z^2/(2*n))/den; half=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/den;
    lo=max(0,center-half); hi=min(1,center+half);
end
function G=make_groups(state,c)
    G=struct('frequency',table(),'rate',table());
    k=find(strcmp({state.conditions.name},'EbN0_minus7'),1);
    if isempty(k); return; end
    m=metrics_table(state.records{k},state.done{k});
    G.frequency=one_group(m,m.TrueF_Hz-c.fc,[-45000 -30000 -15000 0 15000 30000 45000]);
    G.rate=one_group(m,m.TrueRate_HzPerS,[-900 -600 -300 0 300 600 900]);
end
function T=one_group(m,x,edges)
    rows=cell(numel(edges)-1,1);
    for j=1:numel(edges)-1
        valid=x>=edges(j)&x<edges(j+1);
        if j==numel(edges)-1; valid=x>=edges(j)&x<=edges(j+1); end
        r=summarize_metrics(m(valid,:)); r.Lower=edges(j); r.Upper=edges(j+1);
        r.Center=(edges(j)+edges(j+1))/2; rows{j}=r;
    end
    T=struct2table(vertcat(rows{:}));
end
function T=paired_table(state)
    % 每行保留配对编号，不把复用条件的数据再合并为独立样本。
    T=table(); b=find(strcmp({state.conditions.name},'EbN0_minus7'),1);
    if isempty(b) || ~state.completed(b); return; end
    baseline=metrics_table(state.records{b},state.done{b});
    for k=1:numel(state.conditions)
        if k==b || ~state.completed(k) || strcmp(state.conditions(k).stage,'scan'); continue; end
        m=metrics_table(state.records{k},state.done{k});
        assert(height(m)==height(baseline) && isequal(m.BaseSeed,baseline.BaseSeed));
        n=height(m);
        part=table(repmat(string(state.conditions(k).name),n,1),m.TrialID,m.BaseSeed, ...
            m.FinalSuccess-baseline.FinalSuccess,abs(m.FineF_Error)-abs(baseline.FineF_Error), ...
            abs(m.FineRate_Error)-abs(baseline.FineRate_Error),m.TrackRMSE_Hz-baseline.TrackRMSE_Hz, ...
            m.FinalValid & baseline.FinalValid,'VariableNames',{'ComparedCondition','TrialID', ...
            'BaseSeed','SuccessDifference','AbsoluteFErrorDifference','AbsoluteRateErrorDifference', ...
            'TrackRMSEDifference','BothValid'});
        if isempty(T); T=part; else; T=[T;part]; end %#ok<AGROW>
    end
end

%% 十五、数值保存：同目录临时文件写完后替换，旧检查点不被半写文件破坏
function atomic_state_save(state,path)
    folder=fileparts(path); temp=[tempname(folder) '.mat'];
    guard=onCleanup(@()delete_if_exists(temp)); %#ok<NASGU>
    save(temp,'state','-v7.3');
    [ok,msg]=movefile(temp,path,'f'); assert(ok,'原子替换失败：%s',msg);
end
function atomic_bundle_save(bundle,path)
    folder=fileparts(path); temp=[tempname(folder) '.mat'];
    guard=onCleanup(@()delete_if_exists(temp)); %#ok<NASGU>
    save(temp,'bundle','-v7.3');
    [ok,msg]=movefile(temp,path,'f'); assert(ok,'条件结果替换失败：%s',msg);
end
function delete_if_exists(path)
    if exist(path,'file'); delete(path); end
end
function report=export_numeric_tables(state,folder)
    names={'summary.csv','groups_frequency_minus7.csv','groups_rate_minus7.csv','paired_differences.csv'};
    tables={state.summary,state.groups.frequency,state.groups.rate,paired_table(state)};
    for k=1:numel(state.conditions)
        if any(state.done{k})
            names{end+1}=[state.conditions(k).name '_trials.csv']; %#ok<AGROW>
            tables{end+1}=metrics_table(state.records{k},state.done{k}); %#ok<AGROW>
        end
    end
    report=struct('success',0,'failure',0,'details',{{}});
    for j=1:numel(names)
        path=fullfile(folder,names{j}); temp=[tempname(folder) '.csv'];
        try
            writetable(tables{j},temp);
            [ok,msg]=movefile(temp,path,'f'); assert(ok,'%s',msg);
            report.success=report.success+1;
        catch ME
            delete_if_exists(temp); report.failure=report.failure+1;
            report.details{end+1}=sprintf('%s: %s',names{j},ME.message);
            warning('FSM:CSV','CSV失败但数值MAT已保存：%s',report.details{end});
        end
    end
end

%% 十六、绘图：各图独立try/catch；不改变任何原始误差或强制单调
function [figs,errors]=create_figures(state)
    figs=cell(9,1); errors=cell(9,1); S=state.summary;
    for number=1:9
        try
            h=figure('Color','w','Name',sprintf('FSM v3 figure %d',number)); figs{number}=h;
            switch number
                case 1
                    if ~isfield(state,'demo'); title('No demo record'); continue; end
                    d=state.demo;
                    values={d.data,d.spread,d.tx,d.clean,d.rx};
                    labels={'随机信息比特','Gold扩频和采样','发送实通带BPSK','多普勒信道无噪声输出','含噪接收：Eb/N0=-7 dB'};
                    for j=1:5
                        subplot(3,2,j); plot(d.t_us,values{j}); title(labels{j}); xlabel('时间(us)'); grid on;
                    end
                case 2
                    if ~isfield(state,'demo') || isempty(fieldnames(state.demo.spectrum)); title('No valid spectrum'); continue; end
                    d=state.demo.spectrum;
                    surf(d.fAxis/1e3,d.alphaAxis/1e3,abs(d.S),'EdgeColor','none');
                    xlabel('谱频率 f (kHz)'); ylabel('循环频率 alpha=2fc (kHz)');
                    zlabel('|S_x^alpha(f)| (power/Hz)'); title('真正FSM循环谱：复平均/频率平滑后取模'); view(45,35); colorbar;
                case 3
                    q=S(strcmp(S.Stage,'snr'),:);
                    subplot(1,2,1); errorbar(q.EbN0_dB,100*q.FinalJoint, ...
                        100*(q.FinalJoint-q.FinalLo),100*(q.FinalHi-q.FinalJoint),'o-');
                    hold on; plot(q.EbN0_dB,100*q.CoarseJoint,'s-'); yline(90,'k--');
                    ylim([0 105]); xlabel('Eb/N0(dB)'); ylabel('联合成功率(%)'); legend('最终及95%区间','粗估计','90%'); grid on;
                    subplot(1,2,2); plot(q.EbN0_dB,q.FallbackCount,'o-',q.EbN0_dB,q.InvalidFinalCount,'s-');
                    xlabel('Eb/N0(dB)'); ylabel('独立试验数'); legend('回退','最终无效'); grid on;
                case 4
                    q=S(strcmp(S.Stage,'snr'),:);
                    subplot(2,2,1); semilogy(q.EbN0_dB,q.CoarseF_RMSE,'o-',q.EbN0_dB,q.FineF_RMSE,'s-'); ylabel('频率RMSE(Hz)'); legend('粗','最终'); grid on;
                    subplot(2,2,2); semilogy(q.EbN0_dB,q.CoarseRate_RMSE,'o-',q.EbN0_dB,q.FineRate_RMSE,'s-'); ylabel('变化率RMSE(Hz/s)'); grid on;
                    subplot(2,2,3); errorbar(q.EbN0_dB,q.FineF_Bias,q.FineF_Std,'o-'); ylabel('频率偏差+/-标准差'); grid on;
                    subplot(2,2,4); errorbar(q.EbN0_dB,q.FineRate_Bias,q.FineRate_Std,'o-'); ylabel('变化率偏差+/-标准差'); grid on;
                case 5
                    q=S(strcmp(S.Stage,'duration'),:);
                    subplot(1,2,1); plot(q.Duration_ms,q.FinalJoint,'o-'); xlabel('观测时长(ms)'); ylabel('配对联合成功率'); grid on;
                    subplot(1,2,2); semilogy(q.Duration_ms,q.TrackRMSE_Hz,'o-'); xlabel('观测时长(ms)'); ylabel('整段轨迹RMSE(Hz)'); grid on;
                case 6
                    q=S(strcmp(S.Stage,'sensitivity'),:);
                    subplot(2,1,1); bar(q.FinalJoint); set(gca,'XTick',1:height(q),'XTickLabel',cellstr(q.Condition)); xtickangle(30); ylabel('联合成功率'); grid on;
                    subplot(2,1,2); bar([q.TrackRMSE_Hz q.TrackMaxError_Hz]);
                    set(gca,'XTick',1:height(q),'XTickLabel',cellstr(q.Condition)); xtickangle(30);
                    ylabel('轨迹误差(Hz)'); legend('整段RMSE','最大误差'); grid on;
                case 7
                    q=S(strcmp(S.Stage,'candidates'),:);
                    subplot(1,2,1); bar([q.CoarseJoint q.FinalJoint]); set(gca,'XTickLabel',cellstr(q.Condition)); legend('粗','最终'); ylabel('配对联合成功率'); grid on;
                    subplot(1,2,2); bar(q.TrackRMSE_Hz); set(gca,'XTickLabel',cellstr(q.Condition)); ylabel('轨迹RMSE(Hz)'); grid on;
                    title('确认数量5/10/全部，与多起点3严格区分');
                case 8
                    q=S(strcmp(S.Stage,'bandwidth'),:);
                    subplot(1,2,1); bar(q.FinalJoint); set(gca,'XTickLabel',cellstr(q.Condition)); ylabel('联合成功率'); grid on;
                    subplot(1,2,2); hold on;
                    for j=1:height(q)
                        index=find(strcmp({state.conditions.name},char(q.Condition(j))),1);
                        rec=state.records{index}{1}; f=rec.filter;
                        [H,F]=freqz(f.Numerator,1,4096,f.InputFs_Hz*f.Up);
                        plot(F/1e3,20*log10(max(abs(H)/abs(sum(f.Numerator)),1e-12)));
                    end
                    xlim([0 1000]); ylim([-100 5]); xlabel('实际物理频率(kHz)'); ylabel('实际FIR幅频响应(dB)');
                    legend(cellstr(q.Condition)); grid on;
                case 9
                    index=find(strcmp({state.conditions.name},'Jerk_5000'),1);
                    rec=state.records{index}{1}; m=rec.metrics; c=state.conditions(index).cfg;
                    half=(round(c.T*c.scfFs)-1)/(2*c.scfFs); u=linspace(-half,half,501);
                    actual=m.TrueF_Hz+m.TrueRate_HzPerS*u+.5*c.jerk*u.^2;
                    fitted=m.FineF_Hz+m.FineRate_HzPerS*u;
                    coarse=m.CoarseF_Hz+m.CoarseRate_HzPerS*u;
                    subplot(2,2,1); plot((u+half)*1000,actual-c.fc,'k', ...
                        (u+half)*1000,fitted-c.fc,'r--');
                    xlabel('段内时间(ms)'); ylabel('多普勒(Hz)'); legend('真实Jerk轨迹','最终线性模型'); grid on;
                    subplot(2,2,2); plot((u+half)*1000,coarse-actual,'b:', ...
                        (u+half)*1000,fitted-actual,'r');
                    xlabel('段内时间(ms)'); ylabel('轨迹误差(Hz)'); legend('粗','最终'); grid on;
                    title(sprintf('TrackRMSE=%.3g Hz; 连续最大误差=%.3g Hz',m.TrackRMSE_Hz,m.TrackMaxError_Hz));
                    base=find(strcmp({state.conditions.name},'EbN0_minus7'),1);
                    t=metrics_table(state.records{base},state.done{base});
                    subplot(2,2,3); histogram(t.CoarseF_Error); hold on; histogram(t.FineF_Error);
                    xlabel('频率误差(Hz)'); legend('粗','最终'); grid on;
                    subplot(2,2,4); histogram(t.CoarseRate_Error); hold on; histogram(t.FineRate_Error);
                    xlabel('变化率误差(Hz/s)'); legend('粗','最终'); grid on;
            end
        catch ME
            errors{number}=ME.message; warning('FSM:Figure','第%d图绘制失败：%s',number,ME.message);
        end
    end
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

%% 十七、MATLAB原生自检入口：真实运行状态保存到state.internalChecks/ioChecks
function report=internal_checks(c)
    assert(isequal([c.Rb c.Rc c.fs c.fc c.T c.scfFs c.L c.hop], ...
        [1e5 31e5 24.8e6 1e5 .1 8e5 128 64]));
    assert(isequal(c.fRange,[50e3 150e3]) && isequal(c.rateRange,[-1000 1000]));
    assert(c.fTolerance==20 && c.rateTolerance==200);
    assert(c.multiStarts>=1 && c.multiStarts==floor(c.multiStarts));
    r=receiver_config(c);
    assert(~isfield(r,'fd') && ~isfield(r,'bits') && ~isfield(r,'jerk') && ~isfield(r,'fdTest'));
    blank=empty_search(800,r);
    assert(numel(blank.candidateColumns)==size(blank.candidates,2));
    assert(numel(blank.refinementColumns)==size(blank.refinements,2));
    dF=.5/.1; dA=2/.1^2;
    ids=distinct_indices([100000 0;100001 5;100010 0;100000 300],dF,dA);
    assert(isequal(ids,[1;3;4]));
    P=[8 0 0;0 9 0;0 0 7];
    B=branch_candidates(P,[90000 100000 110000],[-1000 0 1000],3,dF,dA,1);
    assert(size(B,1)==3 && numel(unique(B(:,1)))==3);
    coarse=[100000 0 10];
    R=[100000 0 100010 0 12 0 1 0;100000 0 NaN 0 99 -999 0 0; ...
        100000 0 100020 0 11 1 1 1];
    [point,J,selected]=choose_final(coarse,R,r);
    assert(J==12 && selected==1 && point(1)==100010); % 有效但未收敛与“收敛”不混淆
    [~,J,selected]=choose_final(coarse,[100000 0 100010 0 9 1 1 1],r);
    assert(J==10 && selected==0); % 指标变差回退
    [~,J,selected]=choose_final(coarse,[100000 0 NaN 0 NaN -999 0 0],r);
    assert(J==10 && selected==0); % 数值失败回退
    [~,J,selected]=choose_final(coarse,[100000 0 200000 0 100 1 1 1],r);
    assert(J==10 && selected==0); % 越界不能因分数高而选中
    n=801; fs=8000; u=((0:n-1)-(n-1)/2)/fs; e=.2-20*u-.5*5000*u.^2;
    [tr,tm]=track_error(.2,-20,5000,n,fs);
    assert(abs(tr-sqrt(mean(e.^2)))<1e-10 && tm>=max(abs(e))-1e-12);
    [tr,tm]=track_error(NaN,0,0,n,fs); assert(isinf(tr)&&isinf(tm));
    m=metric_template(); names=fieldnames(m);
    for i=1:numel(names); m.(names{i})=0; end
    m.TrialID=1; m.FineF_Error=NaN; m.TrackRMSE_Hz=Inf;
    records={struct('metrics',m)}; t=metrics_table(records,true);
    result=summarize_metrics(t);
    assert(width(t)==numel(names) && result.Count==1 && result.FinalFailureCount==1);
    assert(result.NonfiniteFinalCount==1 && isinf(result.FineF_RMSE));
    y=cos(2*pi*100173.7*(0:4095)'/c.scfFs);
    [frames,tc,~,w]=prepare_frames(y,r);
    S=fsm_spectrum(frames,tc,0,0,r,w);
    psd=fftshift(mean(abs(fft(frames,[],1)).^2,2))/(c.scfFs*sum(w.^2));
    assert(max(abs(S-psd))<1e-12*max(psd));
    plus=fsm_spectrum(frames,tc,100173.7,413.9,r,w);
    minus=fsm_spectrum(frames,tc,-100173.7,-413.9,r,w);
    assert(max(abs(plus-conj(minus)))<1e-12*max(max(abs(plus)),realmin));
    report=struct('Passed',true,'ParameterPreservation',true,'CandidateDedup',true, ...
        'MultiOutputSelection',true,'Fallback',true,'TrajectoryFormula',true, ...
        'FailureDenominator',true,'TableDimensions',true,'FSM_PSDNormalization',true);
end
function report=io_selftest(folder)
    sub=unique_directory(fullfile(folder,'io_selftest')); mkdir(sub);
    % 先保存数值，然后故意关闭图窗及故意让导图失败，再重新读取数值。
    state=struct('sentinel',42,'done',[true false],'rngState',rng);
    path=fullfile(sub,'numerics_before_figures.mat'); atomic_state_save(state,path);
    restored=load(path,'state');
    assert(sum(restored.state.done)==1 && find(~restored.state.done,1)==2);
    report=struct('NumericsSavedFirst',true,'ReloadWithoutDuplicateCount',true, ...
        'GraphicsAvailable',false,'ClosedFigureFailureCount',NaN,'ForcedExportFailureCount',NaN, ...
        'NumericFileSurvived',false,'GraphicsError','');
    h=[]; h2=[];
    try
        h=figure('Visible','off'); h2=figure('Visible','off');
    catch ME
        close_if_valid(h); close_if_valid(h2);
        report.GraphicsError=ME.message; % 仅图窗创建失败可以标为后端不可用
    end
    if ~isempty(h) && ~isempty(h2) && isgraphics(h,'figure') && isgraphics(h2,'figure')
        guard=onCleanup(@()close_if_valid(h2)); %#ok<NASGU>
        close(h); closed=export_figures({h},sub,false);
        forced=export_figures({h2},sub,true);
        % 断言不放在“后端不可用”的catch里，真实测试失败必须抛出。
        assert(closed.pngFailure==1 && forced.pngFailure==1);
        report.GraphicsAvailable=true;
        report.ClosedFigureFailureCount=closed.pngFailure;
        report.ForcedExportFailureCount=forced.pngFailure;
    end
    loaded=load(path,'state'); assert(loaded.state.sentinel==42);
    report.NumericFileSurvived=true;
end
function close_if_valid(h)
    if ~isempty(h) && isgraphics(h,'figure'); close(h); end
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
