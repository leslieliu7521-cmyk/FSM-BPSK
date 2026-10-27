# FSM正式实验交付

正式入口是`FSM_BPSK_formal_experiment.m`，调用未修改的完整v4引擎，保留v3。它将模式明确设为full，不会误跑默认3次诊断。需要MATLAB R2018b+和Signal Processing Toolbox。

当前云环境没有MATLAB/Octave，因此交付的是正式实验程序，**没有生成正式1000次实验结果，也不宣称已达到−7dB下90%以上**。原v4已完成的Python小样本不是本次正式实验结果。短记录误锁与PN失配下的指标歧义仍是已知限制；本包未偷偷替换成窄先验或未经验证的修复方法。

## 运行

解压后，在MATLAB进入解压目录，先做原生自检：

```matlab
fsm_options=struct('mode','selftest');
run('FSM_BPSK_full_system_v4.m');
```

然后明确启动完整正式实验（可能耗时数小时）：

```matlab
FSM_BPSK_formal_experiment
```

分段验证或恢复：

```matlab
FSM_BPSK_formal_experiment('',2) % 本次最多完成两个条件，不改正式样本数
FSM_BPSK_formal_experiment('实际输出目录/checkpoint.mat',2)
FSM_BPSK_formal_experiment('实际输出目录/checkpoint.mat') % 继续至完成
```

只能恢复本版本full检查点；不能把quick/diagnostic检查点转成正式结果。源文件更改会被引擎拒绝混合恢复。输出目录由控制台打印；保存配置、每条件done/completed、固定种子及随机状态。正式数值先保存，再绘图导出，单图导出失败不破坏已保存数值。

## 固定实验定义

| 项目 | 正式配置 |
|---|---|
| 链路 | 随机数据→Gold扩频→采样→实通带BPSK→多普勒/变化率信道→实AWGN→已知PN接收机 |
| 核心 | 频移FFT→谱交叉乘积→复数跨帧平均→频率平滑→最后取模 |
| 采样/载波/信息率 | 24.8MHz / 100kHz / 100kbit/s |
| 扩频码片率 | 3.1Mchip/s |
| FSM默认 | 800kHz、FFT128、hop64、17点复数平滑 |
| 主观测时间 | 100ms |
| 固定接收先验 | 频率50–150kHz；线性变化率±1000Hz/s |
| 随机信道 | 频偏±45kHz、变化率±900Hz/s；不提供给接收搜索 |
| 主Eb/N0 | −10、−8、−7、−6、−4dB，信息比特定义 |
| 正式联合成功 | 频率误差≤20Hz且变化率误差≤200Hz/s且输出有效 |
| 含噪试验次数 | 每条件1000次；失败和非有限输出不删除 |
| 基准搜索 | 原V3每分支3峰、确认10、多起点3，保留粗候选回退 |
| 统计 | 粗细RMSE/偏差/标准差、联合成功率/Wilson95%、轨迹误差、收敛/有效/回退分开 |

完整计划包括265个无噪声扫描点（每点1次）及24个含噪条件（每条件1000次）：5个Eb/N0、3个候选确认数、4个时长、10个同步/Jerk/Raw FSM条件、2个实际接收带宽。四个相同基线条件复用记录，实际需要约20000次含噪搜索；不能把复用或跨条件配对当作额外独立样本。

已知PN准确同步是理想基准；失配对照的rxRatePPM只改变本地PN索引，不是ADC时钟偏差。恒定相移分辨率0.125码片。Raw FSM为800kHz带限无PN基线。小范围DirectExact仅属于单独诊断模式，不参与本正式全先验实验。

## 读取与判定

```matlab
T=FSM_BPSK_formal_report('实际输出目录/results.mat');
L=load('实际输出目录/results.mat');
S=L.state.summary; % 所有扫描与含噪条件
k=find(strcmp({L.state.conditions.name},'EbN0_minus7'));
r=L.state.records{k}{1}; % 首次试验原始误差、所有候选/起点与独立评价
r.metrics
r.search
r.diagnosis
```

每个条件输出numeric.mat和trials.csv；总汇总在summary.csv，配对变化在paired_differences.csv，原始图在figure_*.png/.fig。读取函数不改变任何实验结果。

“1000次中观测成功率超过90%”与“Wilson95%下限超过90%”是不同要求。本报告同时显示两者，统计判据采用后者且必须完整1000次；不能把1000次中的901次成功说成已支持总体超过90%。即使达到该判据，也只适用于当前模型与误差门限，不意味着完整盲同步可靠性。

## 本次验证边界

核查正式入口显式full、恢复拒绝非full、含噪条件1000、报告不删除失败、v3/v4哈希不变；核对扫描265点、含噪24条件及4个复用条件。已做可用的源码检查。MATLAB原生运行、运行耗时、检查点恢复和正式统计结果本次均未执行，需在用户MATLAB环境实际运行。
