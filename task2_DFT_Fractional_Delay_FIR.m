%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 名称：教材Page179例6.1：小数时延FIR滤波器设计
% 时间：2014.04.08~2014.04.10
% 作者：詹飞
% 功能：假设信号没有确定解析表达式，设计长度为L的小数时延FIR滤波器设计
%       此处并未考虑阵元的个数，因此也不需要形成加权向量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%% part1: 变量初始化 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 线性调频信号LPM参数
f0 = 1000;      % f0
f1 = 0.5*f0;    
fu = f0;        % f1和fu为调频信号源的上下边界频率 
fs = 5*f0;      % 采样频率fs
Ts = 1/fs;      % 采样间隔Ts
N = 512;        % 采样点数512
T = N/fs;       % 信号持续时间，满足T*fs = 512个采样点


% FIR滤波器参数
L = 15;             % 滤波器长度L
tao = 0.12345*Ts;   % 期望延迟量τ
fc = 0.4*fs;        % 通带截止频率fc
f2 = 0.5*fs;        % 数字频带边缘值，无法到达
K = 100;            % 离散化频率点数K

STEP_d = f2/fs/K;    
F_PB = linspace(0,fc/fs,fc/f2*K+1); % 延迟滤波器的通带F_PB = [0,0.4]
fd_ideal = F_PB;
fd = 0:STEP_d:(f2/fs-STEP_d);       % 数字频带fd = [0,0.5),fd = f/fs 全频带

% 反向延迟D个节拍
if (tao>=-0.5*Ts && tao<0.5*Ts && mod(L,2)==1)  % L为奇数，τ∈[-0.5Ts,0.5Ts)
    D = (L-1)/2;
elseif (tao>=0 && tao<0.5*Ts && mod(L,2)==0)    % L为偶数，τ∈[0Ts,0.5Ts)
    D = L/2-1;
elseif (tao>=-0.5*Ts && tao<0 && mod(L,2)==0)   % L为偶数，τ∈[-0.5Ts,0)
	D = L/2+1;
end
   
%%%%%%%%%%%%%%%%%%%%% part2: FIR滤波器的期望频率响应 %%%%%%%%%%%%%%%%%%%%%%%
Hd = exp(-1i*2*pi*fd_ideal*(D+tao/Ts)).'; % fd_ideal = [0,0.4]
Amplitude_Hd = 10*log10(abs(Hd));
Phase_Hd = angle(Hd)/pi*180;

%%%%%%%%%%%%%%%%%% part3: FIR滤波器的设计频率响应 L∞范数准则 %%%%%%%%%%%%%%%
% 误差加权系数λk
% lambda_k = zeros(1,K);
% for i=1:length(F_PB)
%     lambda_k(i) = 1;
% end

% 频率响应向量 e(f) = [1,exp(-j2πf/fs),...,exp(-j(L-1)2πf/fs)].' 
%        |      1              ...        1                |
% e(f) = | exp(-j2πfd(1))     ...   exp(-j2πfd(81))      |  
%        |     ...             ...        ...              |
%        |exp(-j(L-1)2πfd(1)) ...   exp(-j(L-1)2πfd(81)) |
%       
e_f = exp(-1i*2*pi*kron((0:L-1).',fd_ideal));
e_f_full = exp(-1i*2*pi*kron((0:L-1).',fd));

% L∞-norm:optimal Chebyshev filter formulation
cvx_begin
    variable h_Inf(L,1)
    minimize(max(abs(e_f.'*h_Inf - Hd)))
cvx_end

% check if problem was successfully solved
disp(['Problem is ' cvx_status])
if ~strfind(cvx_status,'Solved')
  h_Inf = [];
end

% 使用L∞范数准则设计出的FIR滤波器
H_Inf = e_f_full.'*h_Inf;
Amplitude_H_Inf = 10*log10(abs(H_Inf));
Phase_H_Inf = angle(H_Inf)/pi*180;

%%%%%%%%%%%%%%%%%% part4: FIR滤波器的设计频率响应 L1范数准则 %%%%%%%%%%%%%%%
% L1-norm
cvx_begin
    variable h_1(L,1)
    minimize(norm(e_f.'*h_1 - Hd,1))
cvx_end

% check if problem was successfully solved
disp(['Problem is ' cvx_status])
if ~strfind(cvx_status,'Solved')
  h_1 = [];
end

% 使用L1范数准则设计出的FIR滤波器
H_1 = e_f_full.'*h_1;
Amplitude_H_1 = 10*log10(abs(H_1));
Phase_H_1 = angle(H_1)/pi*180;

%%%%%%%%%%%%%%%%%% part5: FIR滤波器的设计频率响应 L2范数准则 %%%%%%%%%%%%%%%
% L2-norm
cvx_begin
    variable h_2(L,1)
    minimize(norm(e_f.'*h_2 - Hd,2))
cvx_end

% check if problem was successfully solved
disp(['Problem is ' cvx_status])
if ~strfind(cvx_status,'Solved')
  h_2 = [];
end

% 使用L1范数准则设计出的FIR滤波器
H_2 = e_f_full.'*h_2;
Amplitude_H_2 = 10*log10(abs(H_2));
Phase_H_2 = angle(H_2)/pi*180;

%%%%%%%%%%%%%%%%%%%% part6: 理想延迟波形与FIR滤波器延迟波形 %%%%%%%%%%%%%%%%%
t = (0:1:N-1)./fs;
s_ideal_LPM = zeros(1,length(t));
s_delay_LPM = zeros(1,length(t));
for i=1:length(t)
    s_ideal_LPM(i) = sin(2*pi*(f1+(fu-f1)/(2*T)*t(i))*t(i));
	s_delay_LPM(i) = sin(2*pi*(f1+(fu-f1)/(2*T)*(t(i)-tao))*(t(i)-tao));
end

% L∞,L1,L2准则小数时延滤波器
y_FIR_delay = zeros(3,length(t)+L-1);       % 卷积后输出信号长度512+L-1
y_FIR_delay(1,:) = conv(s_ideal_LPM,h_Inf);	% L∞
ww=filter(h_Inf,1,s_ideal_LPM);
y_FIR_delay(2,:) = conv(s_ideal_LPM,h_1);   % L1
y_FIR_delay(3,:) = conv(s_ideal_LPM,h_2);   % L2

% 反向延迟D=7个节拍
y_FIR_delay_valid = zeros(3,length(t));
for i=1:length(t)
    y_FIR_delay_valid(1,i) = y_FIR_delay(1,D+i); % L∞
    y_FIR_delay_valid(2,i) = y_FIR_delay(2,D+i); % L1
    y_FIR_delay_valid(3,i) = y_FIR_delay(3,D+i); % L2
end

%%%%%%%%%%%%%%%%%%%%%% part7: plot all figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure1: FIR滤波器幅度响应 
subplot(2,1,1);
plot(fd_ideal,Amplitude_Hd,'or',...                     % FIR期望响应
                           'MarkerEdgeColor','r',...
                           'MarkerFaceColor','r',...
                           'MarkerSize',2);
hold on;
plot(fd,Amplitude_H_Inf,':k','LineWidth',2);    % L∞
hold on;
plot(fd,Amplitude_H_1,'--b','LineWidth',2);     % L1
hold on;
plot(fd,Amplitude_H_2,'-g','LineWidth',2);      % L2
legend('期望值','{\itL}_∞准则','{\itL}_1准则','{\itL}_2准则','Location','SouthWest');
ylabel('幅度/dB');
axis([0 0.5 -1.5 0.25]);

% figure2: FIR滤波器相位响应
subplot(2,1,2);
plot(fd_ideal,Phase_Hd,'or',...                     % FIR期望响应
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','r',...
                       'MarkerSize',2);
hold on;
plot(fd,Phase_H_Inf,':k','LineWidth',2);        % L∞  
hold on;
plot(fd,Phase_H_1,'--b','LineWidth',2);         % L1
hold on;
plot(fd,Phase_H_2,'-g','LineWidth',2);          % L2
xlabel('数字频率{\itf}_d');
ylabel('相位/(°)');
axis([0 0.5 -200 200]);

% figure3: 设计误差
figure;
plot(fd_ideal,abs(H_Inf(1:length(F_PB))-Hd),':k','LineWidth',2);    % L∞
hold on;
plot(fd_ideal,abs(H_1(1:length(F_PB))-Hd),'--b','LineWidth',2);     % L1
hold on;
plot(fd_ideal,abs(H_2(1:length(F_PB))-Hd),'-g','LineWidth',2);      % L2
legend('{\itL}_∞准则','{\itL}_1准则','{\itL}_2准则','Location','NorthWest');
xlabel('{\itf}_d');
ylabel('|{\itH}({\itf}_d)-{\itH}_d({\itf}_d)|');
axis([0 0.4 0 8*10^(-3)]);

% figure4: 输出波形对比
figure;
subplot(4,1,1);
plot(t*fs,s_delay_LPM);     % 理想延迟波形
ylabel('理想延迟波形');
axis([0 N-1 -1.2 1.2]);

subplot(4,1,2);
plot(t*fs,y_FIR_delay_valid(1,:));   % L∞ 准则输出
ylabel('{\itL}_∞准则');
axis([0 N-1 -1.2 1.2]);

subplot(4,1,3);
plot(t*fs,y_FIR_delay_valid(2,:));   % L1准则输出
ylabel('{\itL}_1准则');
axis([0 N-1 -1.2 1.2]);

subplot(4,1,4);
plot(t*fs,y_FIR_delay_valid(3,:));   % L2准则输出
xlabel('\iti');
ylabel('{\itL}_2准则');
axis([0 N-1 -1.2 1.2]);

% 各种准则下，波形误差
figure;
subplot(3,1,1);
plot(t*fs,y_FIR_delay_valid(1,:)-s_delay_LPM);   % L∞ 准则
axis([0 N-1 -5*10^(-3) 5*10^(-3)]);
ylabel('{\itL}_∞准则波形误差');

subplot(3,1,2);
plot(t*fs,y_FIR_delay_valid(2,:)-s_delay_LPM);   % L1准则
axis([0 N-1 -5*10^(-3) 5*10^(-3)]);
ylabel('{\itL}_1准则波形误差');

subplot(3,1,3);
plot(t*fs,y_FIR_delay_valid(3,:)-s_delay_LPM);   % L2准则
axis([0 N-1 -5*10^(-3) 5*10^(-3)]);
xlabel('\iti');
ylabel('{\itL}_2准则波形误差');
