# FSM v3 运行、恢复和读取结果

主文件：`FSM_BPSK_full_system_v3.m`，所有局部函数已包含在一个文件内。不要误运行保留的旧版`FSM_BPSK_full_system.m`。需要MATLAB R2018b+及Signal Processing Toolbox；中文图形显示还需要系统有合适的中文字体。

## 1. 先运行快速自检（不做Monte Carlo）

在文件所在目录执行：

```matlab
fsm_options=struct('mode','selftest');
run('FSM_BPSK_full_system_v3.m');
```

检查参数、候选距离、多输出选择/回退、表列数、非有限样本分母、轨迹公式、FSM归一化；另先保存数字42，再故意关闭图窗/注入导图失败并重读数字，验证故障不覆盖数值。没有图形后端时在`state.ioChecks`中标为未验证，不把跳过算通过。自检断言失败会抛出，不能被“后端不可用”catch隐藏。

## 2. 快速配对验证（默认）

```matlab
fsm_options=struct('mode','quick','quickTrials',3);
run('FSM_BPSK_full_system_v3.m');
```

有21个无噪声点；其余条件每个3次，主要检查运行、诊断和保存，不能用于宣称90%统计性能。包括多起点，仍可能耗时数分钟，视机器而定。可设`makePlots=false`跳过正式绘图；原生I/O自检由`ioSelfTest`单独控制。

也可以分段验证恢复机制：

```matlab
fsm_options=struct('mode','quick','maxNewConditions',2,'makePlots',false);
run('FSM_BPSK_full_system_v3.m');
% 查看控制台输出的结果目录，再恢复：
fsm_options=struct('resumeFile','实际结果目录/checkpoint.mat', ...
    'maxNewConditions',2,'makePlots',false);
run('FSM_BPSK_full_system_v3.m');
```

已保存TrialID不会再追加；只有未完成条件/试验或未完成的条件保存会继续。源文本变化会阻止混合恢复；请保留同一文件。恢复可改绘图开关、I/O自检、检查点间隔和本次条件预算，不能改变科学配置、随机seed和试验次数。

## 3. 正式实验：明确选择full，1000次/条件

```matlab
fsm_options=struct('mode','full','checkpointEveryTrials',50);
run('FSM_BPSK_full_system_v3.m');
```

正式主实验仍100ms，5个Eb/N0点；加上配对候选、时长、敏感性和接收带宽对照，约20000次需要实际计算的含噪估计及265个无噪声点；另有4个配对条件直接复用基线，避免重复计算。可能耗时数小时，需要足够内存和磁盘保存完整候选诊断。默认程序不会自动运行full。

开始新实验请删除工作区旧的`fsm_options`或明确重新赋值；不要用quick检查点恢复成full。程序只清理工作区变量，保留`fsm_options`，不会关闭用户已有其他图窗。

可配置搜索参数在主文件第一节：`branchPeaks`、`confirmCandidates`、`multiStarts`、去重尺度因子及最大精搜索迭代数。默认每分支3候选、确认10个、多起点3个。改代码后请新建实验，不混合旧检查点。

## 4. 结果文件

新建唯一目录`FSM_v3_时间戳`，已有同名目录自动加后缀，不覆盖其他实验。

- `checkpoint.mat`：配置、条件、逐次完成标记、记录及随机状态；用于恢复。
- `results.mat`：先于正式图形保存的数值结果，之后补入CSV/图形导出报告。
- `条件名_numeric.mat`：该条件完整配置、逐次指标、候选排名、精确分数、所有起点的结果和FIR信息。
- `条件名_trials.csv`：逐次41列标量指标。
- `summary.csv`：每条件粗/最终统计、失败数、回退数、收敛/有效输出计数、轨迹误差、Wilson区间及Completed状态。
- `groups_*_minus7.csv`：按真实频偏/变化率分组，保留空组与样本量。
- `paired_differences.csv`：相对−7dB基线按TrialID配对的误差/成功/轨迹差值；没有删去无效配对。
- `figure_*.png/.fig`：正式图片；失败不会阻止数值结果保存。
- `io_selftest*/`：原生保存故障自检，不是科学仿真数据。

## 5. 读取方法

```matlab
A=load('实际结果目录/results.mat');
S=A.state;
disp(S.summary);
disp(S.numericExport); disp(S.imageExport);

k=find(strcmp({S.conditions.name},'EbN0_minus7'),1);
R=S.records{k}{1};
disp(R.metrics);                 % 参数、误差、轨迹、回退、边界等
C=array2table(R.search.candidates,'VariableNames',R.search.candidateColumns);
F=array2table(R.search.refinements,'VariableNames',R.search.refinementColumns);
disp(C); disp(F);                % 未确认候选的ExactScore为NaN，Confirmed为0
plot(R.filter.Numerator);        % 本次实际重采样FIR
```

自检模式读取`S.internalChecks`与`S.ioChecks`；此时没有Monte Carlo汇总表。

最终统计不以优化器exitflag替代估计成功；查看`FinalSuccess`、`FinalValid`、`SelectedOptimizerConverged`、`RefineFallback`等不同字段。回退时`OptimizerExitFlag=NaN`，每个起点的exitflag仍在F表。若完整样本含非有限估计，RMSE报Inf且有计数；不要删掉这些行再报告一个漂亮的RMSE。
