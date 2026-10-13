% 修改说明（2026-10-13）
% 保持Rb,Rc,fs,t1,Ns,fc,g1,g2、实cos调制、SNR、numMonte和粗细搜索参数。
% 保留用户原mseq函数依赖，不擅自更换PN码生成约定。
% 修正相位负号；时域准确频移FFT替代复FFT插值；窗能量归一化。
% 新增循环相关曲面/切片、循环谱曲面/切片、干净/含噪波形。
% V2：每次独立粗细估计后统计RMSE/Bias/Std；新增+50Hz网格对照与无噪声检查。
% 需要MATLAB R2018b或更高版本（xline）及原mseq.m；噪声由randn生成。
% 运算量较大：11个SNR，每个30次观测，含两套网格对照；重叠候选复用。

clear all; clc; close all;
rng(20261013,'twister'); % 固定种子便于重复核对

%% 接收数据参数
Rb = 1e5;        % 数据速率
Rc = 31e5;       % 伪码速率
fs = 2 * 124e5;    % 采样频率
ts = 1/fs;
length_gold = 31; 
t1 = 1e-3;
length_data = Rb*t1;
Ns = 8;
N = fs*t1;
t = 0:1/fs:(N-1)/fs;

fc = 1e5; % 载波频率

%% 产生平衡Gold码
g1 = [0 0 1 0 1]; 
g2 = [1 1 1 0 1]; 
m_1 = mseq(g1);
m_2 = mseq(g2);
gold_seq = zeros(33, 31); % 两条m序列 + 31条异或序列
for shift=0:30
    m_2shifted = circshift(m_2, [0, shift]);
    gold_seq(shift+3,:) = 2*(mod(m_1+m_2shifted, 2))-1;
end
m_1 = 2*m_1-1;
m_2 = 2*m_2-1;
gold_seq(1,:) = m_1;
gold_seq(2,:) = m_2;
gold1 = gold_seq(4,:);

%% 待发送数据
s1 = randi([0,1], 1, length_data);
s1(s1==0) = -1;

%% 扩频调制
ssp_data1 = zeros(1, length_data*length_gold);
for i = 1:length_data
    ssp_data1(((i-1)*length_gold+1):i*length_gold) = s1(i)*gold1;
end

%% 内插采样
chip_index = floor((0:N-1)/(fs/Rc)) + 1;
ssp_sample1 = ssp_data1(chip_index); % 整数采样位置，避免浮点边界错位

%% 实通带BPSK调制：保持原信号模型，2fc特征使用普通循环谱
% ssp_data_bpsk = ssp_sample1 .* exp(1j*2*pi*fc*t); 
ssp_data_bpsk = ssp_sample1 .* cos(2*pi*fc*t);

%% ==============FSM参数设置================
SNR = -10:1:0;
num_snr = length(SNR);

%% ================= FSM分段参数 =================

% 原来的 Nfft=248 会导致：
% df = fs/Nfft = 100 kHz
% 对1 Hz精搜索明显过粗
%
% 这里增大分段长度
frame_length = 4096;

% FFT点数
Nfft = 4096;

% 50%重叠
overrate = 0.5;

% 帧移
hop = round(frame_length * (1-overrate));

% 分段数
numFrame = floor((N-frame_length)/hop) + 1;

% 实际使用的数据长度
Nsample_use = (numFrame-1)*hop + frame_length;

fprintf('\n========== FSM参数 ==========\n');
fprintf('采样率 fs = %.3f MHz\n', fs/1e6);
fprintf('FFT长度 = %d\n', Nfft);
fprintf('FFT频率间隔 = %.3f Hz\n', fs/Nfft);
fprintf('帧长度 = %d\n', frame_length);
fprintf('Hop = %d\n', hop);
fprintf('帧数 = %d\n', numFrame);

%% ================= Monte Carlo =================

numMonte = 30;

%% ================= 循环频率搜索 =================

% 理论循环频率：
%
%       alpha0 = 2*fc
%
alpha_true = 2*fc;

% 粗搜索范围
alpha_coarse_min = max(0, alpha_true - 5000);
alpha_coarse_max = min(fs/2, alpha_true + 5000);

% 粗搜索步长
alpha_coarse_step = 100;

% 精搜索步长
alpha_fine_step = 1;

% 精搜索范围
fine_half_range = 500;

alpha_axis_coarse = ...
    alpha_coarse_min:alpha_coarse_step:alpha_coarse_max;

%% ================= 频率轴 =================

% fftshift以后：
%
%       -fs/2 ~ fs/2
%
f_axis = (-Nfft/2:Nfft/2-1) * (fs/Nfft);

% 只搜索载波附近的频率区域
%
% 对于alpha≈2fc，
% f+alpha/2 和 f-alpha/2
% 分别落在 ±fc附近
%
f_search_max = 500e3;

f_eval = f_axis(abs(f_axis) <= f_search_max).'; % 与FFT频点严格一致

%% ================= 窗函数 =================

window = hamming(frame_length);

snr_target = -7; % 保持原来的展示SNR
x_demo = [];

%% ================= 独立估计 + 配对Monte Carlo =================
% 不平均J后找峰：每次观测独立粗搜和细搜，最后统计误差。
% 主实验仍用原网格；对照仅偏移粗搜索网格50 Hz，信号参数不变。
alpha_grid_offset = alpha_axis_coarse + 50;
alpha_grid_offset = alpha_grid_offset(alpha_grid_offset <= alpha_coarse_max);
coarse_grids = {alpha_axis_coarse, alpha_grid_offset};
case_names = {'Original coarse grid','Coarse grid shifted +50 Hz'};
alpha_coarse_trials = zeros(num_snr,numMonte,2);
alpha_fine_trials = zeros(num_snr,numMonte,2);
% 各SNR复用同一份标准高斯噪声，只改变幅度，减少相邻点随机差异。
% 不同mc独立；粗细方法和网格对照始终使用完全相同的数据。
noise_bank = randn(numMonte,N);
signal_power = mean(abs(ssp_data_bpsk).^2);
frame_start = (0:numFrame-1).'*hop + 1;

clean_frames = make_frames(ssp_data_bpsk,frame_start,frame_length,window);
[c_clean,f_clean,diag_clean] = estimate_two_grids(clean_frames,coarse_grids, ...
    alpha_fine_step,fine_half_range,f_axis,f_eval,frame_start,fs);
fprintf('\n无噪声检查：真实载波 %.3f Hz\n',fc);
for c=1:2
    fprintf('%s: coarse %.3f Hz, fine %.3f Hz, fine bias %.3f Hz\n', ...
        case_names{c},c_clean(c)/2,f_clean(c)/2,f_clean(c)/2-fc);
end

for k=1:num_snr
    noise_scale=sqrt(signal_power/10^(SNR(k)/10));
    for mc=1:numMonte
        x0_awgn=ssp_data_bpsk + noise_scale*noise_bank(mc,:);
        frames=make_frames(x0_awgn,frame_start,frame_length,window);
        [c_hat,f_hat,diag_one]=estimate_two_grids(frames,coarse_grids, ...
            alpha_fine_step,fine_half_range,f_axis,f_eval,frame_start,fs);
        alpha_coarse_trials(k,mc,:)=reshape(c_hat,1,1,2);
        alpha_fine_trials(k,mc,:)=reshape(f_hat,1,1,2);
        if SNR(k)==snr_target && mc==1
            x_demo=x0_awgn;
            diag_demo=diag_one; % 同一次观测的粗细指标曲线
        end
        fprintf('SNR=%d dB, trial %d/%d\n',SNR(k),mc,numMonte);
    end
end

error_coarse=alpha_coarse_trials/2-fc;
error_fine=alpha_fine_trials/2-fc;
rmse_coarse=reshape(sqrt(mean(error_coarse.^2,2)),num_snr,2);
rmse_fine=reshape(sqrt(mean(error_fine.^2,2)),num_snr,2);
bias_coarse=reshape(mean(error_coarse,2),num_snr,2);
bias_fine=reshape(mean(error_fine,2),num_snr,2);
std_coarse=reshape(std(error_coarse,0,2),num_snr,2);
std_fine=reshape(std(error_fine,0,2),num_snr,2);
fine_win_rate=reshape(mean(abs(error_fine)<abs(error_coarse),2),num_snr,2);
% 标准差使用样本定义；RMSE直接对每次误差统计。
for c=1:2
    result_table=table(SNR(:),rmse_coarse(:,c),rmse_fine(:,c), ...
        bias_coarse(:,c),bias_fine(:,c),std_coarse(:,c),std_fine(:,c), ...
        fine_win_rate(:,c),'VariableNames',{'SNR_dB','Coarse_RMSE_Hz', ...
        'Fine_RMSE_Hz','Coarse_Bias_Hz','Fine_Bias_Hz','Coarse_Std_Hz', ...
        'Fine_Std_Hz','Fine_Win_Rate'});
    fprintf('\n%s（固定同一数据序列，仅噪声Monte Carlo）\n',case_names{c});
    disp(result_table);
end

%% ================= 统计图：不强制平滑、不强制单调 =================
figure('Name','Carrier frequency RMSE');
for c=1:2
    subplot(1,2,c);
    plot(SNR,rmse_coarse(:,c),'o-','LineWidth',1.3); hold on;
    plot(SNR,rmse_fine(:,c),'s-','LineWidth',1.3);
    xlabel('SNR (dB)'); ylabel('Carrier frequency RMSE (Hz)');
    title(case_names{c}); legend('Coarse','Fine','Location','best'); grid on;
end
figure('Name','Bias and standard deviation');
for c=1:2
    subplot(2,2,c);
    plot(SNR,bias_coarse(:,c),'o-',SNR,bias_fine(:,c),'s-');
    xlabel('SNR (dB)'); ylabel('Bias (Hz)'); title(case_names{c});
    legend('Coarse','Fine','Location','best'); grid on;
    subplot(2,2,c+2);
    plot(SNR,std_coarse(:,c),'o-',SNR,std_fine(:,c),'s-');
    xlabel('SNR (dB)'); ylabel('Standard deviation (Hz)');
    legend('Coarse','Fine','Location','best'); grid on;
end
figure('Name','Fine search strict win rate');
plot(SNR,fine_win_rate(:,1),'o-',SNR,fine_win_rate(:,2),'s-');
xlabel('SNR (dB)'); ylabel('P(|fine error| < |coarse error|)'); ylim([0 1]);
legend(case_names,'Location','best'); grid on;
title('Strict improvement rate; ties are not counted as wins');

% 同一个样本上比较干净和含噪指标，定位有限记录偏差。
for c=1:2
    figure('Name',sprintf('Search diagnostics case %d',c));
    subplot(2,1,1);
    plot(diag_clean(c).coarse_axis/1e3,diag_clean(c).coarse_J/max(diag_clean(c).coarse_J),'k-'); hold on;
    plot(diag_demo(c).coarse_axis/1e3,diag_demo(c).coarse_J/max(diag_demo(c).coarse_J),'b-');
    xline(alpha_true/1e3,'--r'); xlabel('alpha (kHz)'); ylabel('Normalized J');
    title([case_names{c},' - coarse']); legend('Clean','Noisy','True 2fc'); grid on;
    subplot(2,1,2);
    plot(diag_clean(c).fine_axis/1e3,diag_clean(c).fine_J/max(diag_clean(c).fine_J),'k-'); hold on;
    plot(diag_demo(c).fine_axis/1e3,diag_demo(c).fine_J/max(diag_demo(c).fine_J),'b-');
    xline(alpha_true/1e3,'--r'); xlabel('alpha (kHz)'); ylabel('Normalized J');
    title(sprintf('Fine - one observation at %d dB',snr_target));
    legend('Clean','Noisy','True 2fc'); grid on;
end
% 可直接查看工作区的error_*、rmse_*、bias_*、std_*和diag_clean。

%% ---------- 3. 类似cyclo_bpsk_demo的信号、循环相关和循环谱图 ----------
% 注意：信号参数不同，图的形式对应demo，峰形不应强行做成相同。
assert(~isempty(x_demo), '展示SNR必须包含在SNR列表中。');
alpha_plot = unique([0:1e4:5e5, Rb, 2*fc, Rc-Rb, Rc, Rc+Rb]);
max_lag = 64;  % 绘图范围：覆盖8个码片，属于分析参数
lags = -max_lag:max_lag;
R = cyclic_correlation(x_demo, lags, alpha_plot, fs);

figure('Name','DSSS-BPSK signal');
show_n = min(N, 80*Ns);
subplot(2,1,1);
plot(t(1:show_n)*1e6, ssp_data_bpsk(1:show_n)); grid on;
xlabel('Time (us)'); ylabel('Amplitude'); title('Clean real DSSS-BPSK');
subplot(2,1,2);
plot(t(1:show_n)*1e6, x_demo(1:show_n)); grid on;
xlabel('Time (us)'); ylabel('Amplitude'); title(sprintf('One observation: SNR = %d dB',snr_target));

figure('Name','Cyclic autocorrelation surface');
% 非均匀alpha轴使用真实坐标，不把Rc附近特征挤进低频曲面。
a_low = alpha_plot <= 5e5;
surf(lags/fs*1e6, alpha_plot(a_low)/1e3, abs(R(a_low,:)), 'EdgeColor','none');
xlabel('Lag (us)'); ylabel('Cyclic frequency alpha (kHz)'); zlabel('|R_x^alpha(tau)|');
title(sprintf('Cyclic autocorrelation, single observation, %d dB',snr_target));
view(45,35); colorbar; grid on;

figure('Name','Cyclic autocorrelation slices'); hold on;
alpha_slices = [0 Rb 2*fc Rc];
for a = alpha_slices
    row = find(alpha_plot==a,1);
    plot(lags/fs*1e6,abs(R(row,:)), 'LineWidth',1.2, ...
        'DisplayName',sprintf('alpha = %.1f kHz',a/1e3));
end
xlabel('Lag (us)'); ylabel('|R_x^alpha(tau)|');
title('Cyclic autocorrelation slices'); legend('Location','best'); grid on;

% 时间平滑(跨帧相干平均) + 频率平滑(复谱移动平均)。
% 先平均复谱，最后取模；不将PSD峰值冒充循环谱。
freq_smooth_bins = 9;
starts_demo = (0:numFrame-1)*hop + 1;
frames_demo = zeros(frame_length,numFrame);
for i=1:numFrame
    frames_demo(:,i) = x_demo(starts_demo(i):starts_demo(i)+frame_length-1).' .* window;
end
f_display_max = Rc + 2*fc;
mask_display = abs(f_axis) <= f_display_max;
f_display = f_axis(mask_display);
Smap = zeros(length(alpha_plot),sum(mask_display));
for ia=1:length(alpha_plot)
    spec = cyclic_spectrum_exact(frames_demo,alpha_plot(ia),starts_demo,fs,window);
    spec = conv(spec,ones(freq_smooth_bins,1)/freq_smooth_bins,'same');
    Smap(ia,:) = spec(mask_display).';
end

figure('Name','Cyclic spectrum surface');
surf(f_display/1e6,alpha_plot(a_low)/1e3,abs(Smap(a_low,:)), 'EdgeColor','none');
xlabel('Spectral frequency f (MHz)'); ylabel('Cyclic frequency alpha (kHz)');
zlabel('|S_x^alpha(f)| (power/Hz)');
title(sprintf('Smoothed cyclic spectrum, single observation, %d dB',snr_target));
view(45,35); colorbar; grid on;

figure('Name','Cyclic spectrum slices'); hold on;
for a=alpha_slices
    row=find(alpha_plot==a,1);
    plot(f_display/1e6,abs(Smap(row,:)), 'LineWidth',1.2, ...
        'DisplayName',sprintf('alpha = %.1f kHz',a/1e3));
end
xlabel('Spectral frequency f (MHz)'); ylabel('|S_x^alpha(f)| (power/Hz)');
title('Cyclic spectrum slices: alpha=0 is PSD'); legend('Location','best'); grid on;

fprintf('\n载波估计 = alpha_hat/2；1 Hz网格不是1 Hz精度保证。\n');
fprintf('保留原mseq依赖：请将原来的mseq.m放在MATLAB路径中。\n');

%% ================= 局部函数 =================
function J = cyclic_metric_exact(frames,f_axis,f_eval,alpha,frame_start,fs)
    % frames每列已经加窗；不再对复数FFT做线性插值。
    spec = cyclic_spectrum_exact(frames,alpha,frame_start,fs,[]);
    mask = f_axis(:)>=min(f_eval) & f_axis(:)<=max(f_eval);
    % 对复循环谱做频率平滑后再取模，降低单频点噪声尖峰。
    spec = conv(spec,ones(9,1)/9,'same');
    J = max(abs(spec(mask)));
end

function spec = cyclic_spectrum_exact(frames,alpha,frame_start,fs,window)
    nfft=size(frames,1);
    tau=(0:nfft-1).'/fs;
    Xp=fftshift(fft(bsxfun(@times,frames,exp(-1j*pi*alpha*tau)),nfft,1),1);
    Xm=fftshift(fft(bsxfun(@times,frames,exp( 1j*pi*alpha*tau)),nfft,1),1);
    phase=exp(-1j*2*pi*alpha*(frame_start(:).'-1)/fs);
    products=bsxfun(@times,Xp.*conj(Xm),phase);
    if isempty(window)
        window=0.54-0.46*cos(2*pi*(0:nfft-1).'/(nfft-1));
    end
    spec=mean(products,2)/(fs*sum(window.^2));
end

function R=cyclic_correlation(x,lags,alphas,fs)
    % 对称时延定义；奇数时延的中点允许半个采样间隔。
    x=x(:).'; N=length(x);
    R=zeros(length(alphas),length(lags));
    for il=1:length(lags)
        lag=lags(il);
        i=max(1,1-lag):min(N,N-lag);
        j=i+lag;
        midpoint=((i+j)/2-1)/fs;
        product=x(j).*conj(x(i));
        for ia=1:length(alphas)
            R(ia,il)=mean(product.*exp(-1j*2*pi*alphas(ia)*midpoint));
        end
    end
end

function frames=make_frames(x,starts,frame_length,window)
    frames=zeros(frame_length,length(starts));
    for i=1:length(starts)
        frames(:,i)=x(starts(i):starts(i)+frame_length-1).' .* window;
    end
end

function [coarse_hat,fine_hat,detail]=estimate_two_grids(frames,grids, ...
    fine_step,half_range,f_axis,f_eval,starts,fs)
    % 合并两套候选频率，每个频率仅计算一次；不改变搜索指标。
    union_coarse=unique([grids{1},grids{2}]);
    jc=zeros(size(union_coarse));
    for i=1:length(union_coarse)
        jc(i)=cyclic_metric_exact(frames,f_axis,f_eval,union_coarse(i),starts,fs);
    end
    coarse_hat=zeros(1,2); fine_hat=zeros(1,2);
    detail=struct('coarse_axis',{},'coarse_J',{},'fine_axis',{},'fine_J',{});
    fine_axes=cell(1,2);
    for c=1:2
        [~,loc]=ismember(grids{c},union_coarse);
        detail(c).coarse_axis=grids{c}; detail(c).coarse_J=jc(loc);
        [~,idx]=max(jc(loc)); coarse_hat(c)=grids{c}(idx);
        fine_axes{c}=max(0,coarse_hat(c)-half_range):fine_step: ...
            min(fs/2,coarse_hat(c)+half_range);
    end
    union_fine=unique([fine_axes{1},fine_axes{2}]);
    jf=zeros(size(union_fine));
    [reuse,loc]=ismember(union_fine,union_coarse);
    jf(reuse)=jc(loc(reuse));
    for i=find(~reuse)
        jf(i)=cyclic_metric_exact(frames,f_axis,f_eval,union_fine(i),starts,fs);
    end
    for c=1:2
        [~,loc]=ismember(fine_axes{c},union_fine);
        detail(c).fine_axis=fine_axes{c}; detail(c).fine_J=jf(loc);
        [~,idx]=max(jf(loc)); fine_hat(c)=fine_axes{c}(idx);
    end
end
