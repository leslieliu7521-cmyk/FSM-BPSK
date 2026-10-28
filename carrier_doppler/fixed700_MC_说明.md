# fixed700 第一层 Monte Carlo

新文件 `FSM_BPSK_fixed700_MC.m` 基于已交付的fixed700单记录文件另存，原文件保留。需要MATLAB R2018b+及Signal Processing Toolbox。

## 运行

先做3次流程检查：

```matlab
mc_options=struct('numTrials',3);
run('FSM_BPSK_fixed700_MC.m');
```

再运行默认50次（请明确重置旧选项）：

```matlab
clear mc_options;
run('FSM_BPSK_fixed700_MC.m');
```

或指定种子与次数：

```matlab
mc_options=struct('numTrials',50,'baseSeed',20261015);
run('FSM_BPSK_fixed700_MC.m');
```

不运行上述3次与50次两组，除非分别执行两条命令。每次运行创建新目录，不覆盖原结果，也不把不同运行自动混作一组。

## 改动及固定条件

只扩展统计流程，第一层接收函数、候选生成、去重、评分与原fixed700版本逐函数比对一致。每次以独立mt19937ar种子生成信息比特、初始相位及噪声，种子为 `baseSeed+(0:numTrials-1)`。每条记录独立完成接收机搜索后记录误差，绝不平均多次搜索指标找共同峰。

固定信息比特Eb/N0=-7dB、真实700Hz/s、100ms、中点真实载波100kHz，已知PN及同步不变。噪声仍为实AWGN双边N0/2约定：`mean(clean.^2)*fs/(2*Rb*10^(-7/10))`；目标原始采样SNR约-27.9342dB，逐次CSV包含实际原始采样SNR。生成与评价共用观测段中点0.049999375s。

接收机先验仍为50~150kHz与-1000~1000Hz/s，真实参数不进入接收机白名单。没有连续精搜索调用。没有SNR/窗长/离网格扫描、自适应或新估计器。

700恰在原变化率网格上；100kHz恰在原第一层分支中心频率网格上。单次噪声仍可能导致偏移或失败；本次不人为加扰动。结果仅评价该固定真值条件的重复性，不能外推覆盖所有频率、变化率。

## 输出与统计规则

新建 `FSM_fixed700_MC_时间戳` 目录：

- `trials.csv`：编号、种子、完成/有效标记、±20Hz成功标记、参考时刻、真实及估计f0/变化率及误差、实际原始采样SNR、单次耗时、连续搜索调用数、异常信息。`RuntimeException`单独标明运行异常。
- `summary.csv`：全部N条分母的成功率、Wilson双侧95%区间、频率RMSE/偏差/样本标准差、无效及失败数、运行异常数、变化率RMSE/偏差/样本标准差及耗时。
- `results.mat`：配置、源码、种子、逐次标量、汇总与完成状态。默认不保存完整波形或逐次候选数组。
- `checkpoint.mat`：每条结束后更新。未完成文件的 `complete=false`，未执行行保持 `Completed=false`。本版本不提供自动恢复入口，勿把未完成表当完整统计。
- `figure_01.png/.fig`：默认仅一幅逐次f0误差汇总图，含±20Hz门限，旁边给成功率、区间、RMSE等。无效试验在编号位置标红叉，红叉纵坐标0只是定位，不是零误差。
- `plot_report.mat`：图像导出报告；数值MAT先于CSV和图像保存。

所有已计划试验都执行一次并计入分母，不补抽样。成功为有效且 `abs(f0Error)<=20`；失败数=N−成功数，包含无效。无效/非有限估计的误差按NaN保存。任一完整样本误差非有限时，该指标完整样本RMSE=Inf，偏差/样本标准差=NaN，不删除失败样本重算。全部误差有限时按所有N条计算；样本标准差用N−1分母，N=1时标准差NaN。运行异常保留该行并计入无效/失败，单列异常数和原因，不能将有运行异常的结果视为已验证的估计器性能。

`TotalElapsed_s`是从试验循环前到循环结束的墙钟时间，含检查点保存和循环进度输出，不含最终CSV/图像导出及循环前初始化。`MeanElapsedPerTrial_s=TotalElapsed_s/N`。另有 `TotalCompute_s` 与 `MeanComputePerTrial_s`，只累计各试验计时段（含信号生成及估计，可选诊断计算，排除检查点导出）。

原单记录PSD/循环谱/切片图可选：

```matlab
mc_options=struct('numTrials',50,'singleRecordPlots',true);
run('FSM_BPSK_fixed700_MC.m');
```

只对预先指定的第1条记录计算这些额外数组并绘图，不按误差挑例。默认关闭，不逐次绘图。`makeSummaryPlot=false`可关闭汇总图；图窗有效性检查和兼容导出回退保留。

## 验证边界

当前环境无MATLAB/Octave，未原生执行3次流程检查，也未执行50次试验；没有提供虚构的成功率或估计结果。已通过MATLAB tree-sitter语法解析、第一层函数一致性与统计公式边界核对，共23项；详见 `fixed700_MC_validation.json`。语法解析与Python公式核对不能代替MATLAB原生运行和图形验证。文件中的实际统计值仅在用户运行后产生。
