%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 名称：优化阵列信号处理 P237 例8.4：
% 时间：2022 -10 -17
% 作者：Q！n
% 功能：设计长度为L，具有任意形式的频响的FIR滤波器（仿真中给出期望频响）
%        a unified framework for designing FIR filters with arbitrary
%        magnitude and phase response 
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
L = 64;             % 滤波器长度L


%% 通带和阻带设置 （归一化）
PB_fl=0.15;           % 延迟滤波器的通带F_PB = [0.15,0.275]
PB_fh=0.275;

SB_left_fl=0;         % 延迟滤波器的阻带F_SB = [0,0.0875] U [0.29375,0.49375]
SB_left_fh=0.0875;

SB_right_fl=0.29375;
SB_right_fh=0.49375;

STEP_d = 0.00625;     % 归一化频率的步长

F_PB=PB_fl:STEP_d:PB_fh;
F_SB_left=SB_left_fl:STEP_d:SB_left_fh;
F_SB_right=SB_right_fl:STEP_d:SB_right_fh;

F_SB=[F_SB_left,F_SB_right];    
fd=0:STEP_d:SB_right_fh;

%%%%%%%%%%%%%%%%%%%%% part2: FIR滤波器的期望频率响应 %%%%%%%%%%%%%%%%%%%%%%%
Hd_pass =[0.051627-1.882562i,-1.235489-0.515611i,-0.655313+0.698169i,0.311384+0.612522i,...\
    0.488205-0.066721i,0.062944-0.346711i,-0.222841-0.114074i,-0.117387+0.132750i,...
    0.078966+0.101126i,0.084503-0.052698i,-0.040880-0.076340i,-0.076726+0.032108i...
    0.019898+0.081157i,0.084205-0.002597i,0.018051-0.081697i,-0.071570-0.038550i,...
    -0.055156+0.054157i,0.031210+0.064871i,0.065498-0.005718i,0.019245-0.056405i,...
    -0.037965-0.038671i];  %来自参考文献
Hd_pass=Hd_pass.';
Amplitude_Hd_pass = 20*log10(abs(Hd_pass));
Phase_Hd_pass = angle(Hd_pass)/pi*180;
Hd_stop=zeros(1,length(F_SB));
%%%%%%%%%%%%%%%%%% part3: FIR滤波器的设计频率响应  不同约束准则 %%%%%%%%%%%%%%%

% 频率响应向量 e(f) = [1,exp(-j2πf/fs),...,exp(-j(L-1)2πf/fs)].' 
%        |      1              ...        1                |
% e(f) = | exp(-j2πfd(1))     ...   exp(-j2πfd(81))      |  
%        |     ...             ...        ...              |
%        |exp(-j(L-1)2πfd(1)) ...   exp(-j(L-1)2πfd(81)) |
%       
e_f_pass = exp(-1i*2*pi*(0:L-1).'*F_PB);
e_f_stop = exp(-1i*2*pi*(0:L-1).'*F_SB);
e_f_full = exp(-1i*2*pi*(0:L-1).'*fd);  


error_constraint=0.01;        
lamda_P=ones(1,length(F_PB));  %通带加权系数
lamda_K=ones(1,length(F_SB));  %阻带加权系数

%  阻带峰值误差约束、通带加权最小均方误差优化准则
cvx_begin
variable h_Inf(L,1)
minimize(norm(e_f_pass.'*h_Inf - Hd_pass,2))

subject to
max(lamda_K'.*abs(e_f_stop.'*h_Inf-Hd_stop.')) <= error_constraint;
cvx_end

% check if problem was successfully solved
disp(['Problem is ' cvx_status])
if ~strfind(cvx_status,'Solved')
  h_Inf = [];
end

% 使用L∞范数准则设计出的FIR滤波器
H_Inf = e_f_full.'*h_Inf;
Amplitude_H_Inf = 20*log10(abs(H_Inf));
Phase_H_Inf = angle(H_Inf)/pi*180;



%  通带加权均方误差约束、阻带最大误差最小
error_constraint_2=0.00001;        

cvx_begin
variable h_Inf_2(L,1)
minimize(max(lamda_K'.*abs(e_f_stop.'*h_Inf_2-Hd_stop.')))
subject to
sum(lamda_P'.*((e_f_pass.'*h_Inf_2-Hd_pass).*conj(e_f_pass.'*h_Inf_2-Hd_pass)))<= error_constraint_2;
cvx_end

disp(['Problem is ' cvx_status])
if ~strfind(cvx_status,'Solved')
  h_Inf_2 = [];
end


H_Inf_2 = e_f_full.'*h_Inf_2;
Amplitude_H_Inf_2 = 20*log10(abs(H_Inf_2));
Phase_H_Inf_2 = angle(H_Inf_2)/pi*180;


%%%%%%%%%%%%%%%%%%%%%% part7: plot all figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure1:  
figure();
subplot(2,1,1);
plot(F_PB,Amplitude_Hd_pass,'or',...                     % FIR期望响应
                           'MarkerEdgeColor','r',...
                           'MarkerFaceColor','r',...
                           'MarkerSize',2);
hold on;
plot(fd,Amplitude_H_Inf,':k','LineWidth',1);    % L∞

legend('Desired','Designed','Location','NorthEast');
ylabel('Magnitude(dB)');
xlabel('Normalized Frequency')
title('Least-squares passband subject to peak constrained stopbands')


% figure2: 

subplot(2,1,2);
plot(F_PB,Phase_Hd_pass,'or',...                     % FIR期望响应
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','r',...
                       'MarkerSize',2);
hold on;
plot(fd,Phase_H_Inf,':k','LineWidth',1);         

legend('Desired','Designed','Location','NorthEast');
ylabel('Phase(degree)');
xlabel('Normalized Frequency')


% figure3: 设计误差
figure();
subplot(2,1,1);
plot(F_PB,Amplitude_Hd_pass,'or',...                     % FIR期望响应
                           'MarkerEdgeColor','r',...
                           'MarkerFaceColor','r',...
                           'MarkerSize',2);
hold on;
plot(fd,Amplitude_H_Inf_2,':k','LineWidth',1);    % L∞

legend('Desired','Designed','Location','NorthEast');
ylabel('Magnitude(dB)');
xlabel('Normalized Frequency')
title('Minimax stopbands subject to mean-square constrained passband')


% figure2: FIR滤波器相位响应
subplot(2,1,2);
plot(F_PB,Phase_Hd_pass,'or',...                     % FIR期望响应
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','r',...
                       'MarkerSize',2);
hold on;
plot(fd,Phase_H_Inf_2,':k','LineWidth',1);        % L∞  

legend('Desired','Designed','Location','NorthEast');
ylabel('Phase(degree)');
xlabel('Normalized Frequency')

