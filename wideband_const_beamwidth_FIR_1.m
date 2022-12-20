clc;
clear all;
close all;

%%
M=12;
f0=1000;
fl=f0/2;
fu=f0;
fs=3.125*f0;

Ts=1/fs;
c=1500;
d=c/fu/2*[0:M-1];

Fpb=(0.16:0.005:0.32);                                            %通带
Fsb=[(0:0.005:0.13) (0.35:0.005:0.50)];                           %阻带
% Ftb=[(0.14:0.005:0.15) (0.33:0.005:0.34)];                      %过渡带
Fd=0:0.005:0.5;  

fpb=Fpb*fs;
fsb=Fsb*fs;
% ftb=Ftb*fs;

theta=(-90:2:90);
thetaML=(-8:2:28);
thetaSL_left=(-90:2:-12);
thetaSL_right=(32:2:90);
thetaSL=[thetaSL_left thetaSL_right];

thetaML_index_first=find(theta==thetaML(1));
thetaML_index_end=find(theta==thetaML(end));
thetaML_index=[thetaML_index_first:thetaML_index_end];


thetas=10;
tau = d/c*sin(thetas/180*pi);
%%
%%%求期望波束%%%
SL=-25;
wd=exp(-1i*2*pi*fl*d'*sin(thetas/180*pi)/c)/M;
a=exp(-1i*2*pi*fl*d'*sin(theta/180*pi)/c);   
cbf_p=wd'*a;
energy_cbf_P=20*log10(abs(cbf_p));
     
energy_cbf_PML=energy_cbf_P(:,thetaML_index);
energy_cbf_PSL_right=SL*ones(1,length(thetaSL_right));
energy_cbf_PSL_left=SL*ones(1,length(thetaSL_left));
    

% wd=exp(-1i*2*pi*100*d'*sin(thetas/180*pi)/c)/M;
% a=exp(-1i*2*pi*100*d'*sin(theta/180*pi)/c);   
% cbf_p=wd'*a;
% energy_cbf_P2=20*log10(abs(cbf_p));

figure();
hold on
plot(theta,energy_cbf_P,'k-');
% plot(theta,energy_cbf_P2,'k-');
hold on
scatter(thetaML,energy_cbf_PML,'*');
hold on
plot(thetaSL_right,energy_cbf_PSL_right,'r');
hold on
plot(thetaSL_left,energy_cbf_PSL_left,'r');
legend('常规波束','期望主瓣','期望旁瓣')
title('期望波束')
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-60 3])
xlim([-90 90])
grid on
 %% 频域加权值优化设计
p_pb_ideal=cbf_p(:,thetaML_index);
L=64;
N_fk=length(Fpb);
D=L/2;
for ii=1:M    
   Tm(ii)=-round(tau(ii)/Ts+D)*Ts;    
end  
%% 恒定主瓣响应 -频域各子带加权值优化设计
for ii=1:length(fpb)
    
    p_theta_ML=exp(-1i*2*pi*fpb(ii)*d'*sin(thetaML/180*pi)/c);  
    p_theta_SL=exp(-1i*2*pi*fpb(ii)*d'*sin(thetaSL/180*pi)/c);
 
    cvx_begin
    variable w_f(M,1) complex
    minimize (norm(w_f'*p_theta_ML-p_pb_ideal,2))
       
    subject to
    
    Beam_SL=w_f'*p_theta_SL;
    
    abs(Beam_SL)<=10^(SL/20)
    norm(w_f,2)<=0.4217
    cvx_end
    w_f_all(:,ii)=w_f;
end
 
 %%  FIR 滤波器优化设计
error_constraint=0.01;
e_f_pass = exp(-1i*2*pi*(0:L-1).'*Fpb);
e_f_stop = exp(-1i*2*pi*(0:L-1).'*Fsb);
e_f_full = exp(-1i*2*pi*(0:L-1).'*Fd);  
    
delay_pb_matrix=exp(1i*2*pi*Tm'*Fpb*fs); % 预延迟通带相位矩阵
delay_fd_matrix=exp(1i*2*pi*Tm'*Fd*fs);  %预延迟全频带相位矩阵
Hd_pass=conj(w_f_all).*delay_pb_matrix;
    
for ii=1:M

       Hd_pass_m=Hd_pass(ii,:);
       cvx_begin
       variable h(1,L)
       minimize(norm(h*e_f_pass - Hd_pass_m,2))
       
       subject to
       max(abs(h*e_f_stop)) <= error_constraint
       cvx_end
       h_m(ii,:)=h;
end
    
W_FK=h_m*e_f_pass ; %FIRl滤波器通带频响（未预延迟）
H_m=h_m*e_f_full;  %FIR滤波器全频带响应 （未预延迟）

%% PLOT
figure();
for ii=1:length(fpb)    
    w_ff=w_f_all(:,ii);
    a=exp(-1i*2*pi*fpb(ii)*d'*sin(theta/180*pi)/c);   
    Beam_fre_domain(ii,:)=w_ff'*a;
    plot(theta,20*log10(abs(Beam_fre_domain(ii,:))));
    hold on
end
xlim([-90 90])  
ylim([-60 5])
title('恒定主瓣响应宽带波束图')
ylabel('波束/dB')
xlabel('方位/(°)')

figure();
[degree,normalized_freq]=meshgrid(theta,fpb/fs);
surf(degree,normalized_freq,20*log10(abs(Beam_fre_domain)));
zlim([-60 0])
xlim([-90 90])
xlabel('方位/(^o)')
ylabel('归一化频率')
zlabel('波束/dB')
title('恒定主瓣响应宽带波束图');

figure();
for ii=1:length(fpb)
    w_ff=W_FK(:,ii).*conj(delay_pb_matrix(:,ii));
    a=exp(-1i*2*pi*fpb(ii)*d'*sin(theta/180*pi)/c);   
    Beam_time_domain(ii,:)=w_ff.'*a;
    plot(theta,20*log10(abs(Beam_time_domain(ii,:))));
    hold on 
end
xlim([-90 90])
ylabel('波束/dB')
xlabel('方位/(°)')
title('恒定主瓣响应FIR波束形成器波束图');

figure();
surf(degree,normalized_freq,20*log10(abs(Beam_time_domain)));
zlim([-60 0])
xlim([-90 90])
xlabel('方位/(^o)')
ylabel('归一化频率')
zlabel('波束/dB')
title('恒定主瓣响应FIR波束形成器波束图');

% figure();
% for ii=1:length(fpb)    
%     
%     wd=exp(-1i*2*pi*fpb(ii)*d'*sin(thetas/180*pi)/c)/M;
%     a=exp(-1i*2*pi*fpb(ii)*d'*sin(theta/180*pi)/c);   
%     cbf_p=wd'*a;
%     energy_cbf_P=20*log10(abs(cbf_p));
%     
%     plot(theta,energy_cbf_P);
%     hold on
% end

% xlim([-90 90])  
% ylim([-80 5])
% title('常规宽带波束图')
% ylabel('波束/dB')
% xlabel('方位/(°)')

%% 主瓣逼近均方根误差

fpb_ML_Beam_error_frq_domian=Beam_fre_domain(:,thetaML_index)-repmat(p_pb_ideal,length(fpb),1);
fpb_ML_Beam_error_time_domian=Beam_time_domain(:,thetaML_index)-repmat(p_pb_ideal,length(fpb),1);

for ii=1:length(fpb)
    error_square_fre(ii)=sqrt(fpb_ML_Beam_error_frq_domian(ii,:)*fpb_ML_Beam_error_frq_domian(ii,:)'/length(fpb));
    error_square_time(ii)=sqrt(fpb_ML_Beam_error_time_domian(ii,:)*fpb_ML_Beam_error_time_domian(ii,:)'/length(fpb));    
end

figure();
plot(1:length(fpb),error_square_fre,'k-o');
hold on
plot(1:length(fpb),error_square_time,'b-s');
title('主瓣均方根误差')
xlabel('频率序号')
ylabel('主瓣逼近均方根误差')

 %%  滤波器 设计 PLOT
 
Amplitude_Hd_pass = 20*log10(abs(Hd_pass(2,:)));  % 第二个阵元对应的FIR滤波器的 幅度 和相位 响应
Phase_Hd_pass = angle(Hd_pass(2,:))/pi*180;

Amplitude_H_norm_2 = 20*log10(abs(H_m(2,:)));
Phase_H_norm_2 = angle(H_m(2,:))/pi*180;

figure();
subplot(2,1,1);
plot(Fpb,Amplitude_Hd_pass,'or',...                     % FIR期望响应
                           'MarkerEdgeColor','r',...
                           'MarkerFaceColor','r',...
                           'MarkerSize',2);
hold on;
plot(Fd,Amplitude_H_norm_2,':k','LineWidth',1);         
legend('期望值','设计值','Location','SouthWest');
ylabel('幅度/dB');
xlabel('归一化频率');
ylim([-60 0])
xlim([0 0.5])
title('第2号阵元对应FIR滤波器期望与设计频率响应');

subplot(2,1,2);
plot(Fpb,Phase_Hd_pass,'or',...                     % FIR期望响应
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','r',...
                       'MarkerSize',2);
hold on;
plot(Fd,Phase_H_norm_2,':k','LineWidth',1);        % L∞  
ylabel('相位/(°)');
xlabel('归一化频率');