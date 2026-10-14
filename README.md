# FSM-BPSK

最新完整程序：[FSM_BPSK_full_system_v3.m](carrier_doppler/FSM_BPSK_full_system_v3.m)。原`FSM_BPSK_full_system.m`保持不变，请运行新文件。

- [下载完整v3交付包](FSM_BPSK_full_system_v3_delivery.zip)
- [运行、恢复与读取结果](carrier_doppler/FSM_v3_运行说明.md)
- [详细修改说明](carrier_doppler/FSM_v3_修改说明.md)
- [实际验证报告与未验证项](carrier_doppler/FSM_v3_验证报告.md)

需要MATLAB R2018b+和Signal Processing Toolbox；脚本包含全部局部函数。默认quick小规模验证。正式模式保持原信号参数、100ms主观测、固定先验、20Hz/200Hz/s门限和每条件1000次试验，不默认启动可能耗时数小时的full。

核心保持频域FSM：频移FFT、谱交叉乘积、复数跨帧平均、频率平滑、最后取模。采用分支多候选、与观测时长相关的二维去重、精确FSM确认与默认3起点精搜索。最终指标不低于已确认的最佳粗候选，不代表每次真实误差更小。

新增可恢复检查点、数值先于导图保存、配对候选数/时长/失配/带宽实验和Jerk整段轨迹误差。已知PN同步仍是基线前提，Raw对照明确是800kHz带限无PN分支。

当前环境没有MATLAB/Octave，未完成原生运行或正式千次实验；已做小规模Python交叉检查和源审计，不能将其当作MATLAB完整运行结果。先按运行说明执行selftest，再quick，确认后再选择full。
