function FSM_BPSK_formal_experiment(resumeFile,maxNewConditions)
%% FSM频域循环谱正式实验入口：固定调用v4的full模式，不改变估计器
% 新建完整正式实验：FSM_BPSK_formal_experiment
% 分段运行（先完成两个条件）：FSM_BPSK_formal_experiment('',2)
% 恢复：FSM_BPSK_formal_experiment('实际输出目录/checkpoint.mat')
% 恢复且本次最多两个条件：FSM_BPSK_formal_experiment('实际目录/checkpoint.mat',2)
%
% 一、正式定义
% 24个含噪条件各1000次，另有265个无噪声参数扫描点，每点1次。
% 四个相同基线条件复用原记录，不作为新增独立试验；实际约20000次含噪搜索。
% 不运行诊断模式的固定小范围DirectExact，不将其称为全先验正式估计器。
% 接收先验50--150kHz、变化率+-1000Hz/s；误差门限20Hz、200Hz/s。
% 主观测100ms，另保留1/10/50/100ms配对对照；Eb/N0为信息比特定义。
% FSM核心、已知PN前提、随机种子、逐次失败、Wilson区间及保存逻辑均由原v4保持。
%
% 二、执行范围与限制
% 本入口明确启动可能耗时数小时的正式实验，请先完成v4 selftest/小样本验证。
% MATLAB R2018b+，Signal Processing Toolbox。v3、v4原文件均不修改。
% 此入口不保证达到90%；实际性能由正式输出评价，不强制单调或删除失败。
    if nargin<1 || isempty(resumeFile); resumeFile=''; end
    if nargin<2 || isempty(maxNewConditions); maxNewConditions=Inf; end
    assert(ischar(resumeFile) || (isstring(resumeFile) && isscalar(resumeFile)), ...
        'resumeFile必须是检查点路径或空字符串。');
    assert(isnumeric(maxNewConditions) && isscalar(maxNewConditions) && ...
        maxNewConditions>=1 && (isinf(maxNewConditions) || ...
        maxNewConditions==floor(maxNewConditions)),'条件预算必须为正整数或Inf。');
    engine=fullfile(fileparts(mfilename('fullpath')),'FSM_BPSK_full_system_v4.m');
    assert(exist(engine,'file')==2,'请将本入口与FSM_BPSK_full_system_v4.m放在同一目录。');

    %% 三、只设置运行选项；正式科学参数由未修改的v4配置提供
    fsm_options=struct('mode','full','seed',20261014,'makePlots',true, ...
        'ioSelfTest',true,'checkpointEveryTrials',50,'maxNewConditions',maxNewConditions);
    if ~isempty(resumeFile)
        assert(exist(char(resumeFile),'file')==2,'找不到检查点文件。');
        [pathOK,pathInfo]=fileattrib(char(resumeFile)); assert(pathOK,'无法解析检查点路径。');
        resumeFile=pathInfo.Name; % run会切换到引擎目录，预先转为绝对路径
        header=load(char(resumeFile),'state');
        assert(isfield(header,'state') && strcmp(header.state.opts.mode,'full'), ...
            '正式入口只能恢复full检查点，不能将quick/diagnostic结果冒充正式实验。');
        counts=[header.state.conditions.nTrials];
        scans=strcmp({header.state.conditions.stage},'scan');
        assert(all(counts(~scans)==1000),'检查点不是每个含噪条件1000次的正式配置。');
        clear header; % v4将再次读取完整状态，先释放临时检查副本
        fsm_options.resumeFile=char(resumeFile);
    end
    fprintf('正式FSM实验：full；每个含噪条件1000次；保留原全频率先验和全部失败。\n');
    fprintf('实验完成后调用 FSM_BPSK_formal_report(''输出目录/results.mat'') 检查统计。\n');
    %% 四、最后调用完整引擎；引擎clearvars后不再依赖本函数的局部变量
    run(engine);
end
