clc 
clear all 
close all 

%%  常数定义
deg2rad = pi/180;        % deg -> rad
rad2deg = 180/pi;        % rad -> deg
jay = sqrt(-1);          % 为了避免和循环变量i,j混淆

%% 阵列及 信号参数设置
N_array = 12;           % 阵元数
N_snap = 512;           % 快拍数 
f0=125;
fl = f0/2;              % 入射LFM信号最低频率
fu = f0;               % 入射LFM信号最高频率
fs = 5*f0;              % 采样频率(说明：这里设置的fu和fs恰好能保证信号在相邻阵元间的延迟为一个快拍点)
T = 1/fs;               % 采样时间间隔
t = 0:T:(N_snap-1)*T;   % 采样时间点 
TT = (N_snap-1)*T;      % 信号持续时间
s = sin(2*pi*(fl+(fu-fl)/(2*TT)*t).*t); % 按照书上169页公式构造LFM信号


c = 1500;                      % 声速 
lamda = c/fu;                  % 波长 
d = lamda/2;                   % 阵元间距 
ang_deg = [30];                % 信源方位角(角度)  
ang_rad = ang_deg * deg2rad;   % 信源方位角(弧度)
tau = d/c*sin(ang_rad);        % 相邻阵元间延迟时间，为-0.002，即-T

%% 各个阵元接收信号
for i=1:N_array                
    tt(i,:)=t-(i-1)*tau;
    ttt=tt(i,:);
    ttt(find(ttt<0))=0;
    ttt(find(ttt>TT))=0;
    ss(i,:)= sin(2*pi*(fl+(fu-fl)/(2*TT)*ttt).*ttt);
end

%%  FIR滤波器参数
L = 15;             % 滤波器长度L
f2=0.5*fs;
tao=mod(tau,T);


% 反向延迟D个节拍
if (tao>=-0.5*T && tao<0.5*T && mod(L,2)==1)  % L为奇数，τ∈[-0.5Ts,0.5Ts)
    D = (L-1)/2;
elseif (tao>=0 && tao<0.5*T && mod(L,2)==0)    % L为偶数，τ∈[0Ts,0.5Ts)
    D = L/2-1;
elseif (tao>=-0.5*T && tao<0 && mod(L,2)==0)   % L为偶数，τ∈[-0.5Ts,0)
	D = L/2+1;
end

%% 各子带 归一化频点划分
PB_f1=fl/fs;
PB_fh=fu/fs;
K = 160;                         % 离散化频率点数K
STEP_d = f2/fs/K;    
F_PB = PB_f1:STEP_d:PB_fh; 
fd = 0:STEP_d:f2/fs-STEP_d;      

%% %%%%%%%%%%%%%%%%%%% part2: FIR滤波器的期望频率响应 %%%%%%%%%%%%%%%%%%%%%%%

for ii=1:length(F_PB)
        
     weight_all(:,ii)= exp(-jay*2*pi*F_PB(ii)*fs*d/c*sin(ang_rad).*[0:N_array-1]')./N_array; % 常规波束形成，加权向量即为每个频点处阵列响应向量

end

for ii=1:N_array    
    
    Tm(ii)=-round(tau*(ii-1)/T+D)*T;
    Hd_pass(:,ii)=(conj(weight_all(ii,:)).*exp(jay*2*pi*F_PB*fs*Tm(ii))).';
    
end

%% %%%%%%%%%%%%%%%% part3: FIR滤波器的设计频率响应  2范数约束准则 %%%%%%%%%%%%%%%

% 频率响应向量 e(f) = [1,exp(-j2πf/fs),...,exp(-j(L-1)2πf/fs)].' 
%        |      1              ...        1                |
% e(f) = | exp(-j2πfd(1))     ...   exp(-j2πfd(81))      |  
%        |     ...             ...        ...              |
%        |exp(-j(L-1)2πfd(1)) ...   exp(-j(L-1)2πfd(81)) |
%       
e_f_pass = exp(-1i*2*pi*(0:L-1).'*F_PB);
e_f_full = exp(-1i*2*pi*(0:L-1).'*fd);  


for ii=1:N_array
    Hd_pass_m=Hd_pass(:,ii);
    cvx_begin
        variable h_2(L,1)
        minimize(norm(e_f_pass.'*h_2 - Hd_pass_m,2))
    cvx_end

    % check if problem was successfully solved
    disp(['Problem is ' cvx_status])
    if ~strfind(cvx_status,'Solved')
      h_2 = [];
    end
    
    h_m(ii,:)=h_2;                    % 不同阵元的滤波器系数， 行数为阵元 ，列数为 h(m)
    H_norm_2(:,ii) = e_f_full.'*h_2;  %不同阵元的FIR 响应  列为阵元  行为不同频点响应
    
end

%% filtering 

y_out=zeros(1,N_snap);

for ii=1:N_array
    
   T_m=round(tau*(ii-1)/T+D) ;
   ss_delay =[ ss(ii,1:end) ];
   y_FIR_out=filter(h_m(ii,:),1,ss_delay);
   temp=length(y_FIR_out(T_m+1:end));
   y_add=[y_FIR_out(T_m+1:end) zeros(1,N_snap-temp)];
   y_out=y_out+ y_add;
   
end

%% 画图 

figure()
plot(1:length(t),s); % LFM信号时域波形
xlim([1 N_snap]);
xlabel('(fs = 500Hz,512个样点)'); 
ylabel(''); 
title('LFM信号时域波形') 

ffts = fft(s);
% ffts = ffts/fs;%由于时域抽样会有一个 1/Ts的衰减，所以必须乘以Ts也即除以fs  
df = fs/N_snap;%频率分辨率  
f = [0:df:df*(N_snap-1)] - fs/2;%频率倒转  

%绘制频谱图
figure() 
% plot(f,abs(ffts));%直接画出FFT结果
plot(f,fftshift(abs(ffts)));%先把数据右半部分和左半部分对调成为真实频谱
title('LFM信号幅频特性');  
xlabel('f/Hz')

% plot all channel signal
figure();
max_value=max(ss(1,:));
for ii=1:size(ss,1)
    
    plot(ss(ii,:)+(ii-1)*max_value*2);
    hold on
    
end
 yticks([0:max_value*2:size(ss,1)*max_value*2])
for i=1:size(ss,1)
    ticklabels{i} = num2str(i);
end
yticklabels(ticklabels)
xlabel('时间序号i');
ylabel('阵元序号m');
ylim([-2 24])
xlim([0 N_snap])
title('各阵元接收信号波形')



Amplitude_Hd_pass = 20*log10(abs(Hd_pass(:,2)));  % 第二个阵元对应的FIR滤波器的 幅度 和相位 响应
Phase_Hd_pass = angle(Hd_pass(:,2))/pi*180;

Amplitude_H_norm_2 = 20*log10(abs(H_norm_2(:,2)));
Phase_H_norm_2 = angle(H_norm_2(:,2))/pi*180;

figure()
subplot(2,1,1);
plot(F_PB,Amplitude_Hd_pass,'or',...                     % FIR期望响应
                           'MarkerEdgeColor','r',...
                           'MarkerFaceColor','r',...
                           'MarkerSize',2);
hold on;
plot(fd,Amplitude_H_norm_2,':k','LineWidth',1);         

legend('期望值','设计值','Location','SouthWest');
ylabel('幅度/dB');
xlabel('归一化频率');
ylim([-26 -20])
title('第2号阵元对应FIR滤波器期望与设计频率响应');

subplot(2,1,2);
plot(F_PB,Phase_Hd_pass,'or',...                     % FIR期望响应
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','r',...
                       'MarkerSize',2);
hold on;
plot(fd,Phase_H_norm_2,':k','LineWidth',1);        % L∞  

ylabel('相位/(°)');
xlabel('归一化频率');


figure()
subplot(2,1,1);
plot(y_out);
xlim([0 N_snap])
xlabel('时间序号i')
ylabel('y(i)')
title('FIR波束输出序列及信号源波形失真大小')


subplot(2,1,2);
plot(y_out-s);
xlim([0 N_snap])
ylabel('y(i)-s(i)')
xlabel('时间序号i')
ylim([-2*1e-8 2*1e-8])