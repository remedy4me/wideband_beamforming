%% =======================================================
% By Ning Jiangbo 
% 2018/05/03
% 小数延迟FIR滤波器
%% =======================================================
clc;
clear all;
close all;

%% ----------------- 基本参数
L = 80; % 滤波器长度
N_index = 0:L-1;
fs = 1e4; % 采样频率
Ts = 1./fs;
tao = 0.12345*Ts;
%D_group_delay = D_group_delay_function( L,tao,Ts);% 群时延的点数
if mod(L,2) ==1 && tao>=-0.5*Ts && tao<0.5*Ts
    D_group_delay = (L-1)/2;
elseif mod(L,2) ==0 && tao>=0 && tao<0.5*Ts
    D_group_delay = L/2 - 1;
elseif mod(L,2) ==0 && tao>=-0.5*Ts && tao<0
    D_group_delay = L/2 + 1;
end
%%
PB_scope = [0.15,0.3];% 通带范围
SB_scope = [0,0.1275;0.3225,0.5]; % 阻带范围，第一行：左阻部分和第二行：右阻部分
% SB_scope = [0,0.13;0.32,0.5];
f_pb_inter = 0.005; % 通带相邻间隔
f_sb_inter = 0.0025; % 阻带相邻间隔
f_PB = PB_scope(1):f_pb_inter:PB_scope(2);N_PB = length(f_PB);
f_SB_left = SB_scope(1,1):f_sb_inter:SB_scope(1,2);
f_SB_right = SB_scope(2,1):f_sb_inter:(SB_scope(2,2)-f_sb_inter);
f_SB = [f_SB_left,f_SB_right];N_SB = length(f_SB);
f_digital = [f_SB_left,f_PB,f_SB_right];
f_number = length(f_digital);
% costime = linspace(-pi/4,10*pi/6,N_PB);xxcos = abs(cos(costime))+0.5;
% Hf_pass = xxcos.*exp(-1i*2*pi*f_PB*(D_group_delay+tao/Ts)); % 期望滤波器，调制
Hf_pass = exp(-1i*2*pi*f_PB*(D_group_delay+tao/Ts)); % 期望滤波器(通带恒定‘1’)
Hf_desire = [zeros(size(f_SB_left)),Hf_pass,zeros(size(f_SB_right))];
Hf_stop = [zeros(size(f_SB_left)),zeros(size(f_SB_right))];
figure
subplot(211)
plot(f_PB,20*log10((abs(Hf_pass))),'ro-')
subplot(212)
plot(f_PB,angle(Hf_pass)*180/pi,'ro-')
E_fk = exp(-1i*2*pi*N_index'*f_PB);
E_fp = exp(-1i*2*pi*N_index'*f_SB);
lamda_k = ones(size(f_PB));lamda_p = ones(size(f_SB));
error_constraint = 0.01; % -40dB
%% 2 norm and infinite norm
% 阻带峰值误差约束，通带加权均方误差最小
cvx_begin
variable h2_inf(L)
minimize(lamda_k*((((E_fk.'*h2_inf-Hf_pass.')').').*...
    (E_fk.'*h2_inf-Hf_pass.')))
subject to
max(lamda_p'.*abs(E_fp.'*h2_inf-Hf_stop.')) <= error_constraint;
cvx_end
Hf_design2_inf = [E_fp(:,1:length(f_SB_left)),E_fk,...
    E_fp(:,length(f_SB_left)+1:end)].'*h2_inf;
subplot(211)
hold on;xlabel('fd');ylabel('幅度/dB');title('频率响应')
plot(f_digital,20*log10(abs(Hf_design2_inf)/max(abs(Hf_design2_inf))),'b-')
subplot(212)
hold on;
plot(f_digital,angle(Hf_design2_inf)*180/pi,'b-')
xlabel('fd');ylabel('相位/°');title('相位响应')
%% infinite norm and 2 norm
% 通带加权均方误差约束，阻带最大误差最小
% error_constraint = 0.0210;
cvx_begin
variable hinf_2(L)
minimize(max(lamda_p'.*abs(E_fp.'*hinf_2-Hf_stop.')))
subject to
lamda_k*((((E_fk.'*hinf_2-Hf_pass.')').').*...
    (E_fk.'*hinf_2-Hf_pass.')) <= error_constraint
cvx_end
Hf_designinf_2 = [E_fp(:,1:length(f_SB_left)),E_fk,...
    E_fp(:,length(f_SB_left)+1:end)].'*hinf_2;
subplot(211)
hold on;xlabel('fd');ylabel('幅度/dB');title('频率响应')
plot(f_digital,20*log10(abs(Hf_designinf_2)/max(abs(Hf_designinf_2))),'r-')
subplot(212)
hold on;
plot(f_digital,angle(Hf_designinf_2)*180/pi,'r-')
xlabel('fd');ylabel('相位/°');title('相位响应')
%% 2 norm and 2 norm
% 均方阻带约束，通带均方误差最小
% error_constraint = 1e-5;
cvx_begin
variable h2_2(L)
minimize(lamda_k*((((E_fk.'*h2_2-Hf_pass.')').').*...
    (E_fk.'*h2_2-Hf_pass.')))
subject to
lamda_p*((((E_fp.'*h2_2-Hf_stop.')').').*...
    (E_fp.'*h2_2-Hf_stop.')) <= error_constraint
cvx_end
Hf_design2_2 = [E_fp(:,1:length(f_SB_left)),E_fk,...
    E_fp(:,length(f_SB_left)+1:end)].'*h2_2;
subplot(211)
hold on;xlabel('fd');ylabel('幅度/dB');title('频率响应')
plot(f_digital,20*log10(abs(Hf_design2_2)/max(abs(Hf_design2_2))),'g-')
grid on;
subplot(212)
hold on;
plot(f_digital,angle(Hf_design2_2)*180/pi,'g-')
grid on;
xlabel('fd');ylabel('相位/°');title('相位响应')
legend('期望值','阻带最大误差，通带均方最小','通带均方，阻带最大误差最小',...
    '阻带均方，通带均方最小')
%%
figure
plot(f_digital,abs(Hf_design2_inf-Hf_desire.'),'r-')
hold on;
plot(f_digital,abs(Hf_designinf_2-Hf_desire.'),'c')
plot(f_digital,abs(Hf_design2_2-Hf_desire.'),'g.-')
xlabel('fd');ylabel('设计误差');
legend('阻最大误差约束，通均最小','通均约束，阻最大误差最小',...
    '阻均约束，通均最小')
%% ================LFM调频信号================
N = 512;
t = (0:N-1)/fs;
T = t(end);
fu = 2.7e3;fl = 1.7e3;
sig_LFM = sin(2*pi*(fl+(fu-fl)*t/(2*T)).*t);
figure
t_pie = t - tao;
sig_LFM_delay = sin(2*pi*(fl+(fu-fl)*t_pie/(2*T)).*t_pie);
subplot(311)
plot(t,sig_LFM_delay,'b-');grid on;
xlim([0 inf]);
xlabel('time/s');ylabel('Amplitude');title('理想延迟波形')
y_out = conv(h2_2,sig_LFM);% h2_inf,hinf_2,h2_2
y_out_D_N = y_out(1:N); % 取出滤波器的前N个点
y_out_D = [y_out_D_N(D_group_delay+1:end),zeros(1,D_group_delay)];% 延迟
subplot(312)
plot(t,y_out_D,'r-');grid on;
xlim([0 inf]);
xlabel('time/s');ylabel('Amplitude');title('滤波器延迟波形')
subplot(313)
plot(t,y_out_D-sig_LFM_delay,'r-');grid on;
ylim([-5e-2 5e-2])
xlim([0 inf]);
xlabel('time/s');ylabel('Amplitude');title('波形误差')


