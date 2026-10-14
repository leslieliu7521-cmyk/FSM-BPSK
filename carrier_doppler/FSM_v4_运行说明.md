# FSM v4 运行与结果读取

完整单文件：`FSM_BPSK_full_system_v4.m`。需要MATLAB R2018b+和Signal Processing Toolbox（`resample`、`freqz`）。v3保留。请先小规模执行，查看数值与日志后再扩大。下列MATLAB命令是提供给用户的运行入口，本环境没有MATLAB，未原生执行。

## 1. 先自检，再少量短记录诊断

在文件所在目录执行：

```matlab
fsm_options=struct('mode','selftest');
run('FSM_BPSK_full_system_v4.m');
% 原生参数/候选/回退/轨迹/阶段分类/PN量化/I-O故障自检，不做Monte Carlo。

fsm_options=struct('mode','diagnostic','diagnosticTrials',1);
run('FSM_BPSK_full_system_v4.m');
```

默认smoke：1/10ms、固定离网格90643.7Hz与731.4Hz/s、理想同步/0.125码片/本地PN速率5ppm、无噪声及−7dB、5种方法，60条件，次数取`diagnosticTrials`（默认3）。所有方法都使用固定接收先验，案例真值只给发送端/评价。无噪声条件先执行完才执行含噪条件。可设置`maxNewConditions`分段保存、`makePlots=false`跳过正式图。

## 2. 30/100次配对定位：可只选关注条件

```matlab
fsm_options=struct('mode','diagnostic','diagnosticTrials',30,'makePlots',true);
fsm_options.diagnosticDurations=[.001 .01];
fsm_options.diagnosticScenes={'IdealSync','Rate_5_ppm'};
fsm_options.diagnosticCases={'OffGrid'};
fsm_options.diagnosticMethods={'V3','ConfirmAll','Branch9','DedupHalf','DirectExact'};
fsm_options.diagnosticEbN0=[Inf -7];
run('FSM_BPSK_full_system_v4.m');
% 需要100次时，在新实验中将diagnosticTrials设为100。
```

若只核对含噪失配、且已先验证无噪声，可设`diagnosticEbN0=-7`。不同方法精确重建同seed、同Case/Scene/时长/EbNo的接收记录，先验不依赖真值；TrialID与BaseSeed共同配对。每条件Count是独立基础试验数，不能将方法数乘入成功率分母。

## 3. 覆盖四种时长、分支交界和边界

```matlab
fsm_options=struct('mode','diagnostic','diagnosticProfile','coverage', ...
    'diagnosticTrials',1,'makePlots',false,'maxNewConditions',10);
run('FSM_BPSK_full_system_v4.m');
```

coverage默认时长1/10/50/100ms；案例OffGrid、BranchMinus、BranchPlus、PriorLow、PriorHigh、Random。默认总计720条件，先1次且分段执行，不能默认扩大为100次全部条件。案例由固定参数表注入，方法不接收这些参数。

可显式设置`diagnosticDurations`、`diagnosticScenes`、`diagnosticCases`缩小集合。场景还支持`Phase_0p0625_chip`、`Phase_0p1875_chip`、`Phase_0p5_chip`、`Rate_1_ppm`、`Rate_25_ppm`。

直接精确参考固定搜索90–91.25kHz，因此PriorLow/PriorHigh和多数Random案例范围外；仍保留记录并单独标记，不把范围外失败当作全先验漏峰。100ms直接参考约10521个精确网格点，可能非常耗时。先使用短记录；若要只看长记录生成链路，可将方法集合限制为`{'V3','ConfirmAll','Branch9','DedupHalf'}`。不能围绕每次真值缩小直接参考范围。

## 4. 原quick与正式full

```matlab
fsm_options=struct('mode','quick','quickTrials',3);
run('FSM_BPSK_full_system_v4.m');

% 仅在前面验证正常、明确准备好运行耗时实验后：
fsm_options=struct('mode','full','checkpointEveryTrials',50);
run('FSM_BPSK_full_system_v4.m');
```

quick/full保留原v3条件及V3算法；full含噪Monte Carlo仍每条件1000次（无噪声扫描每点1次）、100ms主实验、原Eb/N0/门限/先验，不自动替换成尚未证实的修复方法。所有模式的逐次MAT均保留独立阶段评价，diagnostic模式默认再对全部生成候选评分并汇总专用诊断表。quick的3/3不代表总体100%可靠，Wilson95%下限约43.8%。

## 5. 暂停与恢复

```matlab
fsm_options=struct('mode','diagnostic','diagnosticTrials',1, ...
    'maxNewConditions',2,'makePlots',false);
run('FSM_BPSK_full_system_v4.m');
% 记录输出目录，然后恢复同一源文件/同一科学配置：
fsm_options=struct('resumeFile','实际输出目录/checkpoint.mat', ...
    'maxNewConditions',2,'makePlots',false,'ioSelfTest',false);
run('FSM_BPSK_full_system_v4.m');
```

恢复时不用重填科学配置；源文本、配置、随机种子、逐次done和completed都在检查点中。只可改绘图/I-O自检/检查点间隔/本次条件预算。已完成TrialID不追加；条件MAT写出失败则重试保存而不重算已标记试验。v3检查点不能混用v4；修改源码后另建实验。数值`results.mat`/`checkpoint.mat`始终先于正式导图保存，单图失败不会中断数值保存。

## 6. 输出与读取

输出目录`FSM_v4_时间戳`：

- `results.mat`/`checkpoint.mat`：完整`state`，含`opts`实际模式/次数、`conditions`、逐次记录、done、completed、随机状态和源文本。
- `summary.csv`：模式、Count/PlannedCount、粗细RMSE/偏差/标准差、联合率/Wilson、失败/无效/收敛/回退/边界、轨迹误差。
- `diagnostic_trials.csv`：所有已完成诊断试验的标量指标、邻域尺度、各阶段频率/联合存在性、首次遗漏原因、计时与PN实际相移。
- `diagnostic_stages.csv`：每条件Count、频率/联合成功率与Wilson、互斥首次阶段计数、范围外参考及非有限计数。
- `diagnostic_pairs.csv`：相同PairKey/TrialID相对V3的原始误差、轨迹、成功与时间差；无效配对仍保留，`BothValid`另列。
- `D*_numeric.mat`及`*_trials.csv`：每条件完整记录/原41个标量指标。
- `figure_*.png/.fig`：逐图独立导出；导出成功/失败数在`state.imageExport`。

```matlab
L=load('实际输出目录/results.mat'); s=L.state;
disp(s.opts.mode); disp(s.summary(:,{'Condition','Mode','Count','PlannedCount','FinalJoint','FinalLo','FinalHi'}));
disp(s.diagnosticStages);
T=s.diagnosticTrials;
bad=T(~logical(T.FinalSuccess),:); % 保留全部失败，不把它们从统计中删除
bad(:,{'Condition','TrialID','FineF_Error','FineRate_Error','TrackRMSE_Hz', ...
    'FirstMissingFrequencyStage','SelectionExplanation','RefineFallback'})

k=find(strcmp({s.conditions.name},char(T.Condition(1))),1);
r=s.records{k}{T.TrialID(1)};
r.search.generated             % [f,rate,分数,FFT分支]，分数种类见generatedScoreKind
r.search.branchDiagnostics     % 各分支候选/频率和率栅格/近似二维表面
r.search.generatedApproxRank   % 每个生成候选的近似排名；DirectExact中实际是精确排名
r.search.generatedKept         % 是否经全局去重保留
r.search.dedupRepresentativeGeneratedIndex % 被哪个原始候选代表
r.search.distinctGeneratedIndices          % candidates行到generated行的映射
r.search.candidateColumns; r.search.candidates % f,a,近似分数,分支,精确分数,近似名次,精确名次,是否确认
r.search.refinementColumns; r.search.refinements % 起点/终点/分数/exitflag/有效/收敛
r.search.exactRanking; r.search.selectedStart; r.search.fallbackReason
r.diagnosis.GeneratedPosthocExact % 按generated原顺序：独立后验精确分数、精确排名
r.diagnosis.TruthProbe            % 事后真值邻域探查；绝不是接收机输出
r.diagnosis.RateProfileAtSelected % 最终频率处的率栅格及精确FSM分数
r.pnModel                        % 明确仅本地PN索引失配、实际采样量化分辨率
s.diagnosticPairs                % 差值<0表示该误差小于v3，不保证每次改善
```

`auditExactAll=false`可减少后验评分耗时，此时`GeneratedPosthocExact`为未评分NaN，`PosthocUnscoredCount`计数；不能将未评分误称数值失败。`saveBranchSurfaces=false`只减少表面存储，分支候选、排名和映射仍保留。设置这些选项要新建实验。

局部循环谱图以最终估计为中心，不能证明全先验捕获。指标高、有效、优化收敛、正式成功须分别读取，`FirstMissingFrequencyStage='None'`也不代表通过正式门限。
