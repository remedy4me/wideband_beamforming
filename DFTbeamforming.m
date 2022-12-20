close all;    %关闭所有正在运行的窗口
clear ;    %清空缓存
clc;          %清屏 命令窗口


%%  parameter
M=12;          %假设空间 阵列的 阵元数目 M=8
d_=0.5;       %因为距离d一般设置为波长λ的一半，所以此处直接令：d/λ=0.5
theta=-30*pi/180;
W=2*pi*d_.*sin(theta);    % 导向矢量中，三个信号源的空间相位
f0=100;       %线性调频信号
f1=f0/2;
fu=f0;
fs=5*f0;
N=512;
T=N/fs;   %总时间长度
i=0:N;
t=i/fs;   %时间点数

s=sin(2*pi*(f1+((fu-f1)/2*T)*t).*t);  %LFM
noise=wgn(M,N+1,0);                   %加噪声



%%  
for m=1:M                       % 大循环为 天线阵元M
     for q=1:N+1                   % 小循环 为每个天线阵元上接收到的1000个采样信号点的值，q从1 到1000，默认步长为1
         Y=[s(q)];  % 构造矩阵Y，其中元素为第q时刻，三个信号源的接收样值。
                                 % m代表第几个阵元，q代表接收信号的时刻点
         f=f1+(fu-f1)*q/(fs*T);
         X(m,q)=Y(1)*exp(-j*(m-1)*W(1)*f/1500);
                                 %X(m，q)代表的是第m个阵元，第q时刻收到的信号
                                 %其中，W（i）是第i个信号源的加权值，
     end
end

 X=X;
 syms XK
 for h=1:M
     eval(['XK',num2str(h),'=fft(X(h,:))']);
 end
 
 XKS=[];
 
for h=1:M
    XKS= eval(['[XKS; XK' num2str(h) ']']);
end

v=1:M;
 theta0=-30;                     % 允许50度，方向的信号通过，其他方向滤除
 theta0=theta0.*pi/180;          % 化为弧度形式
 fy0=2*pi*d_.*sin(theta0);
 
 for q=1:N+1
      f=f1+(fu-f1)*q/(fs*T);
 a0=exp(j*(v-1)*fy0*f/1500);           % 将导向矢量 a(θ)中的元素表示出来，并赋值给a0
 W0(q,1:12)=a0/M;
 end
 
 K=zeros(1,N+1);
 
 for k=1:N+1
    K(k)=W0(k,:)*XKS(:,k);
 end
KI=ifft(K);
 
%%  PLOT 
figure(1)
plot(real(KI));
xlabel('点数');
ylabel('幅度 ');
title('时域波形波束输出');
grid on;                    % 给 当前坐标轴添加主要的格点


deta=abs(s)-abs(KI);
figure(2)
plot(deta);
xlabel('点数');
ylabel('误差 ');
title('DFT波束输出误差情况');
grid on;                    % 给 当前坐标轴添加主要的格点