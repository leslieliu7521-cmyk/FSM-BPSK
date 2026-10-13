# FSM-BPSK

完整实数 DSSS-BPSK 系统仿真：随机数据、Gold 码扩频、采样、BPSK、线性多普勒/变化率、按信息比特 Eb/N0 加入 AWGN，然后用 **频域 FSM 循环谱**估计接收载波与变化率。

- 主程序：[FSM_BPSK_full_system.m](carrier_doppler/FSM_BPSK_full_system.m)
- 说明：[FSM 版本使用与验证说明](carrier_doppler/FSM版本使用与验证说明.md)
- 打包下载：[FSM_BPSK_FSM_frequency_v3.zip](FSM_BPSK_FSM_frequency_v3.zip)
- 旧平方匹配算法对照：[FSM_BPSK_square_baseline.m](carrier_doppler/FSM_BPSK_square_baseline.m)

需要 MATLAB R2018b+ 和 Signal Processing Toolbox。主程序包含全部局部函数，无外部 mseq.m 依赖。先使用 `run_mode='quick'` 检查运行与图形；正式验证使用默认 `run_mode='full'`（1000 次/条件，完整套件可能耗时数小时）。

FSM 核心使用移频 FFT 交叉乘积、复数跨帧平均与谱频率平滑，非零变化率采用候选相位补偿；不是将信号平方的匹配搜索称为循环谱。默认已知 PN 同步前端，单独测试 PN 相位/码率失配、jerk 和取消 PN 的带限对照；并未验证完整盲捕获、调制分类或真实信道。

当前已完成 Python FSM 核心复核和 MATLAB 语法检查；MATLAB 正式试验未运行。100 次 Python −7 dB 结果及 25 点无噪声扫描保存在 `carrier_doppler/fsm_validation_*.json/csv` 和 `fsm_noiseless_scan.csv` 中，不能替代千次 MATLAB 及失配场景验证。所有 MATLAB 图与统计由实际运行生成。
