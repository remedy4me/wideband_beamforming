%% =======================================================
% By Ning Jiangbo 
% 2018/05/01
% 小数延迟FIR滤波器
%% =======================================================
clc;
clear all;
close all;
%%
L = 16; % 滤波器长度
N_index = 0:L-1;
fs = 1e4; % 采样频率
Ts = 1./fs;
tao = 0.12345*Ts;
f_number = 100; % 频点的个数
% D_group_delay = D_group_delay_function( L,tao,Ts);% 群时延的点数
if mod(L,2) ==1 && tao>=-0.5*Ts && tao<0.5*Ts
    D_group_delay = (L-1)/2;
elseif mod(L,2) ==0 && tao>=0 && tao<0.5*Ts
    D_group_delay = L/2 - 1;%
elseif mod(L,2) ==0 && tao>=-0.5*Ts && tao<0
    D_group_delay = L/2 + 1;%
end
f_tem = linspace(0,0.5,f_number+1);
f_digital =f_tem(1:f_number); % 数字频率
f_d_inter = (f_digital(2)-f_digital(1));
f_PB = 0:f_d_inter:0.4;
f_SB =  f_PB(end)+f_d_inter:f_d_inter:f_digital(end);
Hf_pass = exp(-1i*2*pi*f_PB*(D_group_delay+tao/Ts));
% Hf_stop = exp(-1i*2*pi*f_SB*(D_group_delay+0/Ts));
Hf_desire = [Hf_pass,zeros(1,f_number-length(Hf_pass))];
figure
subplot(211)
plot(f_PB,20*log10(round(abs(Hf_pass))),'ro-')
subplot(212)
plot(f_PB,angle(Hf_pass)*180/pi,'ro-')
E_f = exp(-1i*2*pi*N_index'*f_digital);
lamda_k = [ones(1,length(Hf_pass)),zeros(1,f_number-length(Hf_pass))];
%% -------------- norm 1 ,L1范数
cvx_begin
variable h1(L)
minimize(lamda_k*(abs(E_f.'*h1-Hf_desire.')))
cvx_end
Hf_passign1 = E_f.'*h1;
subplot(211)
hold on;
plot(f_digital,20*log10(abs(Hf_passign1)/max(abs(Hf_passign1))),'g-')
subplot(212)
hold on;
plot(f_digital,angle(Hf_passign1)*180/pi,'g-')
%% --------------- norm 2 ,L2范数
cvx_begin
cvx_q = cvx_quiet(true)
variable h2(L)
minimize(lamda_k*((((E_f.'*h2-Hf_desire.')').').*(E_f.'*h2-Hf_desire.')))
cvx_end
Hf_passign2 = E_f.'*h2;
subplot(211)
hold on;xlabel('fd');ylabel('幅度/dB');title('频率响应')
plot(f_digital,20*log10(abs(Hf_passign2)/max(abs(Hf_passign2))),'b-')
subplot(212)
hold on;
plot(f_digital,angle(Hf_passign2)*180/pi,'b-')
xlabel('fd');ylabel('相位/°');title('相位响应')
%% ------------------ norm infinite ,无穷范数
cvx_begin
cvx_q = cvx_quiet(true)
variable h_inf(L)
minimize(max(lamda_k'.*abs(E_f.'*h_inf-Hf_desire.')))
cvx_end
Hf_passign_inf = E_f.'*h_inf;
subplot(211)
hold on;
plot(f_digital,20*log10(abs(Hf_passign_inf)/max(abs(Hf_passign_inf))),'c-')
legend('期望值','L1范数','L2范数','无穷范数')
subplot(212)
hold on;
plot(f_digital,angle(Hf_passign_inf)*180/pi,'c-')
%%
figure
plot(f_digital,abs(Hf_passign1-Hf_desire.'),'r--')
hold on;
plot(f_digital,abs(Hf_passign2-Hf_desire.'),'c')
plot(f_digital,abs(Hf_passign_inf-Hf_desire.'),'g.-')
xlabel('fd');ylabel('设计误差');
legend('L1范数','L2范数','无穷范数')
xlim([0,0.4])
%% ================LFM调频信号================
N = 512;
t = (0:N-1)/fs;
T = t(end);
fu = 2e3;fl = 1e3;
sig_LFM = sin(2*pi*(fl+(fu-fl)*t/(2*T)).*t);
figure
t_pie = t - tao;
sig_LFM_delay = sin(2*pi*(fl+(fu-fl)*t_pie/(2*T)).*t_pie);
subplot(311)
plot(t,sig_LFM_delay,'b-');grid on;
xlim([0 inf]);
xlabel('time/s');ylabel('Amplitude');title('理想延迟波形')
y_out = conv(h1,sig_LFM); % 可选 h1,h2,h_inf
y_out_D_N = y_out(1:N); % 取出滤波器的前N个点
y_out_D = [y_out_D_N(D_group_delay+1:end),zeros(1,D_group_delay)];% 延迟
subplot(312)
plot(t,y_out_D,'r-');grid on;hold on;
plot(t(end-D_group_delay+1:end),zeros(1,D_group_delay),'bp-')
xlim([0 inf]);
xlabel('time/s');ylabel('Amplitude');title('滤波器延迟波形')
subplot(313)
plot(t,y_out_D-sig_LFM_delay,'r-');grid on;
xlim([0 inf]);
ylim([-3e-3 3e-3])
xlabel('time/s');ylabel('Amplitude');title('波形误差')


