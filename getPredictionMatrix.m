%% 计算三阶积分器模型的预测矩阵
function [Tp, Tv, Ta, Bp, Bv, Ba] =getPredictionMatrix(K, dt, p_0, v_0, a_0)
Ta=zeros(3*K);
Tv=zeros(3*K);
Tp=zeros(3*K);

for i=1:K
    Ta(i,1:i)=ones(1,i)*dt;
    Ta(i+K,1+K:i+K)=ones(1,i)*dt;
    Ta(i+2*K,1+2*K:i+2*K)=ones(1,i)*dt;
end

for i=1:K
    for j=1:i
        Tv(i,j)=(i-j+0.5)*dt^2;
        Tv(i+K,j+K)=(i-j+0.5)*dt^2;
        Tv(i+2*K,j+2*K)=(i-j+0.5)*dt^2;
    end
end

for i=1:K
    for j=1:i
        Tp(i,j)=((i-j+1)*(i-j)/2+1/6)*dt^3;
        Tp(i+K,j+K)=((i-j+1)*(i-j)/2+1/6)*dt^3;
        Tp(i+2*K,j+2*K)=((i-j+1)*(i-j)/2+1/6)*dt^3;
    end
end

Ba= [ones(K,1)*a_0(1);ones(K,1)*a_0(2);ones(K,1)*a_0(3)];
Bv= [ones(K,1)*v_0(1);ones(K,1)*v_0(2);ones(K,1)*v_0(3)];
Bp= [ones(K,1)*p_0(1);ones(K,1)*p_0(2);ones(K,1)*p_0(3)];

for i=1:K
   Bv(i)=Bv(i)+i*dt*a_0(1); 
   Bp(i)=Bp(i)+i*dt*v_0(1)+i^2/2*a_0(1)*dt^2;
   
   Bv(i+K)=Bv(i+K)+i*dt*a_0(2); 
   Bp(i+K)=Bp(i+K)+i*dt*v_0(2)+i^2/2*a_0(2)*dt^2;
   
   Bv(i+2*K)=Bv(i+2*K)+i*dt*a_0(3); 
   Bp(i+2*K)=Bp(i+2*K)+i*dt*v_0(3)+i^2/2*a_0(3)*dt^2;
end





