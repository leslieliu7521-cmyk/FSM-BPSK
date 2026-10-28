# fixed700 最小修改说明

基于当前 `FSM_BPSK_spectrum_check.m` 另存完整脚本 `FSM_BPSK_spectrum_check_fixed700.m`。当前程序没有fixed700场景，本文件只运行这个新增场景；原诊断脚本及上传原文件不覆盖。

## 直接运行

MATLAB R2018b+，需要Signal Processing Toolbox。在文件目录执行：

```matlab
clear fixed700_options;
run('FSM_BPSK_spectrum_check_fixed700.m');
```

默认仅运行1条fixed700含噪记录并保留PSD、循环谱局部曲面、固定alpha及固定谱频率切片。可选只改种子、绘图开关或输出父目录，例如：

```matlab
fixed700_options=struct('seed',20261015,'makePlots',false);
run('FSM_BPSK_spectrum_check_fixed700.m');
```

## 改动

1. 固定信息比特Eb/N0=-7 dB、变化率700 Hz/s、标称载波100 kHz、观测段中点多普勒频偏0 Hz、100 ms。没有Monte Carlo、窗长扫描、自适应或第二条对照记录。信号位、相位、标准高斯噪声仍按原种子生成。
2. 解扩/滤波前先计算 `signalPower=mean(clean.^2)`，采用 `noiseVariance=signalPower*fs/(2*Rb*10^(-7/10))`，不对有限长度噪声另作归一化。目标原始采样SNR为 `-7+10*log10(2*Rb/fs)=-27.93421685 dB`。实际值用 `10*log10(signalPower/mean(noise.^2))` 测量，不预填实测结果，也不拿解扩滤波后的SNR冒充原始采样SNR。
3. 生成相位与接收评价统一参考为中心截取记录中点。本参数下该时刻为原始记录起点后的0.049999375 s；生成载波在此刻恰为100 kHz。原始高采样率网格中点与此刻约差0.6微秒，不能只在原始中点固定频偏却换另一个时刻评价。重采样后检查实际长度与参考时刻一致。
4. 第一层分支候选生成、去重、精确评分确认不变；输出其最佳确认候选。移除本文件的连续优化代码和优化辅助函数，不调用fminsearch。`continuousRefinementCalls=0`、`refinements`为空，并加入运行期断言。为兼容原诊断绘图，旧 `fFine/aFine/finalScore` 字段分别等于第一层 `fCoarse/aCoarse/coarseScore`，不表示运行过第二层。
5. 接收机白名单不变：载波先验50~150 kHz、变化率先验-1000~1000 Hz/s、步长100 Hz/s。真实700只用于信号生成及事后误差/标线，未用于缩小先验或选择候选。
6. 图标题标出信息比特Eb/N0、真实/估计载波、真实/估计变化率与参考时刻。曲面/切片仍由第一层估计变化率补偿，局部扫描只用于显示，不更新估计输出。PSD双边轴、坐标定义及单图容错导出保留。

## 输出与验证

输出到新的 `FSM_fixed700_时间戳` 目录。先保存 `spectrum_numeric.mat`，再生成 `spectrum_report.txt` 与各图PNG/FIG及 `plot_report.mat`。数值包含目标Eb/N0、目标和实际原始采样SNR、信号/噪声功率、同一参考时刻的f0及变化率误差、连续搜索调用数，以及原有谱峰/网格数据。单条记录不称为准确率。

当前环境无MATLAB/Octave，**未原生执行，未验证实测估计、峰值或图形结果**。已完成19项静态/代数核对，包含第一层代码一致性、无连续优化调用、单记录入口、噪声公式、参考时刻与目标SNR换算。详见 `fixed700_validation.json`。这不等同于MATLAB语法检查或原生运行通过。

本文件是仅第一层fixed700的单记录诊断，不代表已完成独立性能统计实验。
