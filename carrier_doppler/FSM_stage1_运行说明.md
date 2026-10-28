# 第一层700Hz/s基线运行说明

单文件包含全部MATLAB局部函数，需要MATLAB R2018b+和Signal Processing Toolbox。无需Parallel Computing Toolbox。原文件保留，先前不画图版本另存为`FSM_BPSK_stage1_700_metrics.m`。

## 一条直接运行指令

在文件目录执行：

```matlab
run('FSM_BPSK_stage1_700.m')
```

新工作区默认700Hz/s、100ms、Eb/N0=[−10 −7 −4]dB、每条件50次（150条含噪统计输出），先少量无噪声检查。默认恢复9张直观图并导出PNG，不存FIG/整条波形。若工作区保留以前的`stage1_options`，先清除它，或明确赋新结构，避免沿用上次设置。

## 更改SNR、次数和绘图开关

```matlab
stage1_options=struct('ebN0',[-10 -8 -7 -6 -4],'trials',50,'makePlot',true);
run('FSM_BPSK_stage1_700.m');
```

只做小样本排查：

```matlab
stage1_options=struct('trials',2,'makePlot',false);
run('FSM_BPSK_stage1_700.m');
```

原生自检、100次调试和1000次正式模式：

```matlab
stage1_options=struct('mode','selftest'); run('FSM_BPSK_stage1_700.m');
stage1_options=struct('mode','debug'); run('FSM_BPSK_stage1_700.m'); % 每条件100次
% 确认少量验证正常，再明确运行较大实验：
stage1_options=struct('mode','full'); run('FSM_BPSK_stage1_700.m');  % 每条件1000次
```

full仍仅对已选SNR/窗长/率集合运行，不会自动扩大到全部研究条件。不能在full设置较少次数来冒充正式模式；少量测试用compact（默认）或显式`trials`。

## 开启窗长或离网格变化率对照

```matlab
stage1_options=struct('trials',50,'windowComparison',true,'makePlot',true);
run('FSM_BPSK_stage1_700.m'); % 1/10/50/100ms，共600条含噪输出

stage1_options=struct('trials',50,'rateComparison',true);
run('FSM_BPSK_stage1_700.m'); % 650/700/750Hz/s，700不重复统计，共450条
```

同时开启两项会有1800条输出，勿默认扩大。也可以显式给小集合：

```matlab
stage1_options=struct('trials',30,'durations',[.01 .1],'rates',[650 700 750]);
run('FSM_BPSK_stage1_700.m');
```

可选开关`syncComparison`（理想、0.125码片、本地PN率5ppm）、`bandwidthComparison`（800k/1.6MHz实际重采样并调整FFT）、`candidateComparison`（确认5/10/全部不同生成候选）默认均false。它们是明确标注的方法/模型对照，不能把有限生成候选全部确认说成全局最优。PN率误差不是ADC时钟，0.125码片为整数输入采样点分辨率。

## 测速与优化核对

默认`verifyOptimization=true`：首个无噪声点与前两个基础试验的首个含噪条件，核对参考前端/逐次重建候选版本和缓存/线性滤波版本。允许合理浮点误差，排名或候选选择异常会抛出，先存检查点。

首批实际计算后控制台报告主要耗时项与本机总计算/剩余时间粗估。`Search_s`包括获取/建立帧计划、候选生成和精确确认；`Generation_s`与`Confirmation_s`分项记录。`PreprocessAllocated_s`包含该TrialID基础生成及共享滤波的分摊和本条件裁剪/缩放，`Compute_s`另加搜索和评价，**不含文件保存、绘图、无噪声及额外参考核对**。

`state.wallSeconds`是截至最终数值导出前的累计本次/恢复运行时间，包括初始化上下文、无噪声、核对和中途检查点，排除最终MAT/CSV/PNG导出；`ioSeconds`只累计中途检查点写入，最终导出未包括。`state.benchmarks`保存真实核对耗时。默认不会提前输出未经本机测量的“几分钟完成”。当前交付环境无MATLAB，不能给MATLAB150次真实耗时。

## 暂停与恢复

```matlab
stage1_options=struct('trials',50,'maxNewBaseTrials',5,'makePlot',false);
run('FSM_BPSK_stage1_700.m');
% 用控制台打印的绝对目录：
stage1_options=struct('resumeFile','实际输出目录/checkpoint.mat');
run('FSM_BPSK_stage1_700.m');
```

`maxNewBaseTrials`仅限制本次新完成基础TrialID数，不改计划50次。恢复按条件/TrialID的done标记继续，不重复计数；只允许改绘图开关、检查点间隔和本次预算，科学配置/种子/源文件变更必须新建实验。默认每10个基础TrialID保存小检查点；中途编程/资源异常保留已完成标记再抛出。数值无效估计正常记为失败，不走异常删除路径。

## 输出与读取

目录`FSM_stage1_700_时间戳`：

- `trials.csv`：全部逐次误差、f0/率/联合成功、轨迹、一致性、帧数/候选数/边界、耗时，保留失败。
- `summary.csv`：每条件实际Count/计划Count/模式、Wilson95%、频率与率RMSE/偏差/样本标准差、绝对误差中位数/P90/P95、失败/无效、轨迹和耗时。
- `noiseless.csv`：少量离网格/分支交界检查与率网格误差；默认率700，开启率对照时覆盖650/700/750。
- `results.mat`、`checkpoint.mat`：配置、源文本、条件、种子、完成标记、逐次候选和必要诊断；不含高采样率波形。
- `summary.png`：f0成功率/置信区间/Count及f0 RMSE。
- `simulation_01.png`至`simulation_08.png`：链路波形、局部循环谱、偏差/率/轨迹/耗时、逐次f0与分布/失败、率估计与轨迹、窗长对照、全先验候选、帧间一致性分析。全部在MAT数值保存后生成；默认不存FIG。

```matlab
L=load('实际输出目录/results.mat'); s=L.state;
disp(s.summary);
T=s.trials;
T(:,{'Condition','TrialID','F0_Hz','TrueF0_Hz','F0Error_Hz','F0Success', ...
    'EstimatedRate_HzPerS','RateError_HzPerS','RateSuccess','JointSuccess', ...
    'TrackRMSE_Hz','FrameCoherence','Compute_s'})

r=s.records{1,1};
r.search.generatedColumns; r.search.generated % 原始分支候选与近似分数
r.search.candidateColumns; r.search.candidates % 去重后、精确分数/排名/确认标记
r.search.generatedSortOrder; r.search.distinctGeneratedIndices
r.search.selectedCandidate; r.search.exactRanking
s.failureDiagnostics{1} % 默认每条件至多2条失败的附加每帧复谱诊断
s.benchmarks
assert(s.continuousRefinementCalls==0);
```

`search.generated`及`search.candidates`的Frequency_Hz均以局部窗中点为参考；CSV中的F0_Hz已换算到共同长记录中点，两者不要直接混作同一坐标。

每窗一致性摘要保存在所有逐次输出；需要每个试验的复谱序列时，新实验显式设置`saveFrameSequence=true`。默认仅少量失败诊断保留完整序列。

后续分析应直接使用逐次CSV计算分布、窗口配对变化、一致性与误差关系。没有真值辅助选窗结果或自适应性能曲线；一致性只作分析，不预设有效。50次/条件及其尾部分位数是初步观察，不能因样本成功率高就承诺总体超过90%。

## 恢复原仿真图（本次更新）

新文件默认`makePlot=true, simulationPlots=true`，统计仍是150条输出。明确开启：

```matlab
stage1_options=struct('makePlot',true,'simulationPlots',true);
run('FSM_BPSK_stage1_700.m');
```

已有旧版第一层results.mat时，只绘图、不重跑估计、不覆盖原统计：

```matlab
stage1_options=struct('plotFile','实际目录/results.mat');
run('FSM_BPSK_stage1_700.m');
```

图片存入原结果目录下新建的`simulation_plots_时间戳`，另存小型plot_report.mat。只重建一个预先按配置选择的试验（优先700Hz/s、100ms、−7dB、最早TrialID），不按真实误差挑选示例。局部循环谱以已保存的第一层估计为中心，不能当独立捕获证明。波形只显示约60us短片段，整条高采样率波形不写入MAT。每图独立计数，关闭任一图窗或单图导出失败不影响原数值。

只要原单幅汇总图：`makePlot=true, simulationPlots=false`；完全不画图：`makePlot=false`。默认窗长图仅100ms一个点，需要四窗长才开启windowComparison，不补造数据。新版本检查点不能与旧源码混合；旧检查点请用保留的metrics脚本恢复，得到results.mat后再调用只绘图入口。MATLAB原生图形仍未在本环境验证。
