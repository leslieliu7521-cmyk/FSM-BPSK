# FSM-BPSK

第一层精简基线：[完整程序](carrier_doppler/FSM_BPSK_stage1_700.m)、[下载交付包](FSM_BPSK_stage1_700_delivery.zip)、[运行/恢复与选项](carrier_doppler/FSM_stage1_运行说明.md)、[验证报告](carrier_doppler/FSM_stage1_验证报告.md)。默认700Hz/s、100ms、三个Eb/N0、50基础试验/条件（150条输出），完全关闭连续精搜索，默认不画图。当前无MATLAB原生结果，Python检查不能代替正式统计。

固定变化率+700Hz/s：[下载完整实验包](FSM_BPSK_rate700_delivery.zip)，[完整MATLAB脚本](carrier_doppler/FSM_BPSK_rate700.m)，[运行说明](carrier_doppler/FSM_700Hz每秒实验说明.md)。接收先验仍为±1000Hz/s，主观测100ms，正式每个含噪条件1000次；当前无MATLAB原生统计结果。

正式实验：[下载正式实验包](FSM_BPSK_formal_experiment_delivery.zip)。解压后执行`FSM_BPSK_formal_experiment`，固定full、每个含噪条件1000次；读取/恢复见[FSM正式实验说明](carrier_doppler/FSM_正式实验说明.md)。当前没有MATLAB原生正式结果，不保证达到90%；v3/v4保持不变。

新增完整诊断版：[FSM_BPSK_full_system_v4.m](carrier_doppler/FSM_BPSK_full_system_v4.m)，保留v3与原文件。

- [下载v4完整交付包](FSM_BPSK_full_system_v4_delivery.zip)
- [v4运行、30/100次诊断、full与恢复](carrier_doppler/FSM_v4_运行说明.md)
- [v4修改与诊断尺度](carrier_doppler/FSM_v4_修改说明.md)
- [v4实际验证报告及未完成项目](carrier_doppler/FSM_v4_验证报告.md)

v4默认短记录/本地PN失配诊断：保存分支候选、去重映射、精确排名、多起点和独立真值评价，并对照原V3/全部确认/每分支9峰/半尺度去重/固定小范围直接FSM。尚未证明完整先验误锁已经修复；指标高不等于误差小。成功率图均带Count/Wilson；数值先保存，再逐图导出。

当前环境没有MATLAB/Octave，未原生运行；Python配对诊断不能代替MATLAB验证。先selftest与少量诊断，再决定是否扩大，不默认运行full。PN失配仅改变本地PN索引，非ADC时钟；相移分辨率0.125码片。

以下保留v3说明与下载：

v3完整程序：[FSM_BPSK_full_system_v3.m](carrier_doppler/FSM_BPSK_full_system_v3.m)。原`FSM_BPSK_full_system.m`保持不变，请运行新文件。

- [下载完整v3交付包](FSM_BPSK_full_system_v3_delivery.zip)
- [运行、恢复与读取结果](carrier_doppler/FSM_v3_运行说明.md)
- [详细修改说明](carrier_doppler/FSM_v3_修改说明.md)
- [实际验证报告与未验证项](carrier_doppler/FSM_v3_验证报告.md)

需要MATLAB R2018b+和Signal Processing Toolbox；脚本包含全部局部函数。默认quick小规模验证。正式模式保持原信号参数、100ms主观测、固定先验、20Hz/200Hz/s门限和每条件1000次试验，不默认启动可能耗时数小时的full。

核心保持频域FSM：频移FFT、谱交叉乘积、复数跨帧平均、频率平滑、最后取模。采用分支多候选、与观测时长相关的二维去重、精确FSM确认与默认3起点精搜索。最终指标不低于已确认的最佳粗候选，不代表每次真实误差更小。

新增可恢复检查点、数值先于导图保存、配对候选数/时长/失配/带宽实验和Jerk整段轨迹误差。已知PN同步仍是基线前提，Raw对照明确是800kHz带限无PN分支。

当前环境没有MATLAB/Octave，未完成原生运行或正式千次实验；已做小规模Python交叉检查和源审计，不能将其当作MATLAB完整运行结果。先按运行说明执行selftest，再quick，确认后再选择full。
