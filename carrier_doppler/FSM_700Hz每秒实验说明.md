# 固定700 Hz/s的正式实验

完整独立脚本：`FSM_BPSK_rate700.m`。v3、v4及之前正式实验文件保留不变。需要MATLAB R2018b+及Signal Processing Toolbox。

## 运行

解压后进入该目录，在MATLAB执行：

```matlab
% 先执行原生自检，不做正式试验
fsm_options=struct('mode','selftest');
run('FSM_BPSK_rate700.m');
% 先做少量验证：无噪声频率扫描每点1次，其余每条件3次
fsm_options=struct('mode','quick','quickTrials',3);
run('FSM_BPSK_rate700.m');
% 确认验证正常后明确启动正式实验，可能耗时数小时
fsm_options=struct('mode','full');
run('FSM_BPSK_rate700.m');
```

直接运行脚本而不提供fsm_options时默认为full。重跑时请明确设置新的fsm_options，避免沿用上一次的selftest/quick选项。

分段运行与恢复：

```matlab
fsm_options=struct('mode','full','maxNewConditions',2,'makePlots',false);
run('FSM_BPSK_rate700.m');
fsm_options=struct('resumeFile','实际输出目录/checkpoint.mat');
run('FSM_BPSK_rate700.m');
```

输出目录为`FSM_rate700_时间戳`。不能混用v4或随机变化率版本检查点；恢复会核对版本与源文本，不重复统计已完成试验。

## 实验定义

- **真实线性变化率始终为+700Hz/s**，Jerk=0；100ms内频率漂移约70Hz。
- 每次仍随机生成数据、初相、频偏（±45kHz）和噪声；已知PN同步为模型前提。
- 接收机不知道700，只使用原50–150kHz频率先验和±1000Hz/s变化率先验；没有围绕700缩小搜索。
- 保留完整实通带链路及真正频域FSM核心、多起点与粗候选回退、轨迹误差和数值先存流程。
- 主观测100ms；Eb/N0为信息比特定义，取−10、−8、−7、−6、−4dB。
- 正式联合门限仍为频率误差≤20Hz且变化率误差≤200Hz/s且输出有效，不以收敛代替成功。
- 正式计划53个无噪声频率扫描点（每点1次）；5个Eb/N0条件和4个时长条件（每条件1000次）。100ms时长条件复用−7dB主实验，因此实际需要约8000次新的含噪搜索。跨条件配对不作为新增独立样本。
- 四种时长为1/10/50/100ms，使用相同长记录中心段。为了使整段率恒定700，本版不加入Jerk，也不额外扩展同步/带宽/候选消融研究。

## 图与结果

6张图依次为：完整链路波形、最终估计附近的局部循环谱、Eb/N0粗细成功率及Wilson95%区间与Count、频率/率/轨迹误差、时长对照、−7dB逐次估计率与误差。局部循环谱不能作为独立捕获证明。

```matlab
T=FSM_BPSK_formal_report('实际输出目录/results.mat');
L=load('实际输出目录/results.mat');
k=find(strcmp({L.state.conditions.name},'EbN0_minus7'));
m=L.state.records{k}{1}.metrics;
m.TrueRate_HzPerS        % 应为700
m.FineRate_HzPerS        % 估计值，未强制为700
m.FineRate_Error         % 估计值减700
m.TrackRMSE_Hz
```

formal_report仅用于full结果；quick请直接查看`L.state.summary`。输出保留每次失败/非有限结果、候选列表、精确排名、优化器exitflag、回退及边界诊断。数值MAT先于图保存，导图失败独立计数。

## 实际验证与限制

本次通过25项源码/公式检查：核心接收、评分、选择、评价与保存函数与v4一致；700真值未进入接收配置；先验/门限/正式次数保留；相位二次项拟合对应变化率700，100ms漂移70Hz；MATLAB语法解析通过。检查记录见`FSM_rate700_validation.json`。

当前环境没有MATLAB，**未原生运行或执行1000次正式统计，不能保证−7dB达到90%以上**。语法解析和公式检查不能替代原生运行。短记录误锁和原FSM指标的低信噪比限制仍可能存在，不通过删除失败或强制曲线单调掩盖。
