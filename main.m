clc;
clear all;
close all;
%% 绘制目标圆锥螺旋线轨迹
pi=3.14;
t=0:0.2:500;%总的轨迹时间
t_size=size(t);
w=0.08;%圆锥螺旋线的角速度
v_r=-0.05;%z轴轨迹的速度
v_zr=-0.5*ones(1,t_size(2));
v_xr=0.025*cos(w*t);%x轴轨迹的速度
v_yr=0.025*sin(w*t);%y轴轨迹的速度
x_r=t.*v_xr;%x轴轨迹的位置
y_r=t.*v_yr;%y轴轨迹的位置
z_r=t.*v_r;%z轴轨迹的位置
%% MPC轨迹跟踪
K=20;
dt=0.2;
w_p=100;
w_v=1;
w_a=1;
w_j=1;
p_0=[0,8,0];%x,y,z
v_0=[0,0,0];
a_0=[0,0,0];
log=[0 p_0 v_0 a_0 0 0 0];
i=1;
for t=0:0.2:400
    [Tp, Tv, Ta, Bp, Bv, Ba] =getPredictionMatrix(K, dt, p_0, v_0, a_0);
    Pr=[x_r(i:i+K-1)';y_r(i:i+K-1)';z_r(i:i+K-1)'];
    Vr=[v_xr(i:i+K-1)';v_yr(i:i+K-1)';v_zr(i:i+K-1)'];
    H=w_j*eye(3*K)+w_p*(Tp'*Tp)+w_v*(Tv'*Tv)+w_a*(Ta'*Ta);
    F=w_p*(Bp'-Pr')*Tp+w_v*(Bv'-Vr')*Tv+w_a*Ba'*Ta;
    A=[Tv;-Tv;Ta;-Ta;eye(3*K);-eye(3*K)];
    
    b=[ones(3*K,1)*6-Bv;ones(2*K,1)*6+Bv(1:2*K);ones(K,1)+Bv(2*K+1:3*K);...
        ones(3*K,1)*3-Ba;ones(2*K,1)*3+Ba(1:2*K);ones(K,1)+Ba(2*K+1:3*K);...
        ones(2*K,1)*3;ones(K,1)*2;ones(2*K,1)*3;ones(K,1)*2];
    
    J=quadprog(H,F,A,b);
%     J=quadprog(H,F,[],[]);
    
    jx=J(1);
    jy=J(1+K);
    jz=J(1+2*K);
    
    p_0(1)=p_0(1)+v_0(1)*dt+0.5*a_0(1)*dt^2+1/6*jx*dt^3;
    v_0(1)=v_0(1)+a_0(1)*dt+0.5*jx*dt^2;
    a_0(1)=a_0(1)+jx*dt;
    
    p_0(2)=p_0(2)+v_0(2)*dt+0.5*a_0(2)*dt^2+1/6*jy*dt^3;
    v_0(2)=v_0(2)+a_0(2)*dt+0.5*jy*dt^2;
    a_0(2)=a_0(2)+jy*dt;

    p_0(3)=p_0(3)+v_0(3)*dt+0.5*a_0(3)*dt^2+1/6*jz*dt^3;
    v_0(3)=v_0(3)+a_0(3)*dt+0.5*jz*dt^2;
    a_0(3)=a_0(3)+jz*dt;
    
    log=[log;t p_0 v_0 a_0 jx jy jz];
    i=i+1;
end
figure(1)
subplot(2,2,1)
p_size=size(log);
plot3(x_r(1:p_size(1)),y_r(1:p_size(1)),z_r(1:p_size(1)),'r','linewidth',2)
hold on
plot3(log(:,2),log(:,3),log(:,4),'b','linewidth',2);
legend('ref','traj');
title('position');

subplot(2,2,2)
plot(log(:,1),log(:,5),log(:,1),log(:,6),log(:,1),log(:,7));
legend('x','y','z');
% hold on
% plot(log(1:p_size(1),1),v_xr(1:p_size(1)),log(1:p_size(1),1),v_yr(1:p_size(1)),log(1:p_size(1),1),v_zr(1:p_size(1)));
% legend('x','y','z','x_r','y_r','z_r');
title('v');
subplot(2,2,3)
plot(log(:,1),log(:,8),log(:,1),log(:,9),log(:,1),log(:,10));
legend('x','y','z');
title('a');
subplot(2,2,4)
plot(log(:,1),log(:,11),log(:,1),log(:,12),log(:,1),log(:,13));
legend('x','y','z');
title('jerk');

