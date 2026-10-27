function T=FSM_BPSK_formal_report(resultsFile)
%% 读取正式实验；明确区分完成1000次、观测成功率和Wilson统计判据
% T=FSM_BPSK_formal_report('实际输出目录/results.mat');
% 不删除失败，不以优化器收敛替代真实误差成功，不把复用条件合并成独立样本。
% 本函数只显示/返回原结果及统计判据，不覆盖MAT，也不重新搜索或修改任何估计。
    assert(nargin==1 && exist(char(resultsFile),'file')==2,'请输入results.mat实际路径。');
    L=load(char(resultsFile),'state'); s=L.state;
    assert(strcmp(s.opts.mode,'full'),'该文件不是full正式实验。');
    assert(s.cfg.fTolerance==20 && s.cfg.rateTolerance==200,'误差门限与正式定义不符。');
    T=s.summary(strcmp(s.summary.Stage,'snr'),:);
    assert(height(T)==5 && all(T.PlannedCount==1000),'正式信噪比条件数/计划次数不符。');
    % 仍保留未完成条件，标记为未达到正式判据，不能将已完成的少量试验当1000次。
    T.ThousandTrialsComplete=T.Count==1000 & T.Completed;
    T.ObservedJointAbove90=T.FinalJoint>.90;
    T.WilsonLowerAbove90=T.FinalLo>.90;
    T.StatisticalCriterionMet=T.ThousandTrialsComplete & T.WilsonLowerAbove90;
    disp(T(:,{'Condition','EbN0_dB','Count','PlannedCount','Completed','FinalJoint', ...
        'FinalLo','FinalHi','FinalFailureCount','InvalidFinalCount','FallbackCount', ...
        'FineF_RMSE','FineRate_RMSE','TrackRMSE_Hz','StatisticalCriterionMet'}));
    target=T(T.EbN0_dB==-7,:); assert(height(target)==1);
    fprintf('-7dB：完成%d/1000，联合成功率%.2f%%，Wilson95%%=[%.2f%%, %.2f%%]。\n', ...
        target.Count,100*target.FinalJoint,100*target.FinalLo,100*target.FinalHi);
    if ~target.ThousandTrialsComplete
        fprintf('尚未完成1000次，不能给出正式验收结论，请恢复checkpoint继续。\n');
    elseif target.StatisticalCriterionMet
        fprintf('在本模型和20Hz/200Hz/s门限下，Wilson95%%下限超过90%%。\n');
    else
        fprintf('未达到Wilson95%%下限超过90%%的统计判据，失败样本仍完整计入。\n');
    end
    fprintf('优化收敛、有效输出、回退与真实成功分开统计；结果不外推到完整盲同步。\n');
end
