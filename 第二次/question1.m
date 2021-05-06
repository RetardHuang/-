clear all;clc;close all;
%%
M=2;%m为几,就是几抽头
expectedpower=2;%期望信号功率
%%
% 第（1）问
N=4;%总信号采样数量
noisepower=.5;
figure(1)
% Jw=noisepwoer-2p'*w+w'*R*w%代价函数（误差性能面）
w0=-20:1:20;
w1=-20:1:20;
[w0,w1]=meshgrid(w0,w1);
Jw=(0.5+noisepower)*(w0.^2+w1.^2)+cos(2*pi/N)*w0.*w1+2*w1*sin(2*pi/N)+2;
surf(w0,w1,Jw)
title("第一问的性能曲面")
%%
% 第（2）问
N=4;
noisepower=2;
figure(2)
% Jw=noisepwoer-2p'*w+w'*R*w%代价函数（误差性能面）
w0=-20:1:20;
w1=-20:1:20;
[w0,w1]=meshgrid(w0,w1);
Jw=(0.5+noisepower)*(w0.^2+w1.^2)+cos(2*pi/N)*w0.*w1+2*w1*sin(2*pi/N)+2;
surf(w0,w1,Jw)
title("第二问的性能曲面")

rm=zeros(1,M);%自相关函数,索引1对应着n=0的情况
rm(1)=0.5+noisepower;
rm(2:end)=cos(2*pi*(1:length(rm)-1)/N);
%这里表示比较蠢 比如r(0)=rm（1），r(1)=rm（2），r(-1)=conj（rm（2）)
%这里的自相关函数显然是实数，因此我们可以把上面的conj去掉

pnegm=-sin(2*pi*(0:M-1)/N);%这里是维纳滤波里面的p
p=pnegm';

R=toeplitz(rm);%协方差矩阵，此东西是N*N维的！

wmin2=inv(R)*p;
Jwmin2=expectedpower-2*p'*wmin2+wmin2'*R*wmin2;

%%
% 第（3）问
N=16;
noisepower=0.5;
figure(3)
% Jw=noisepwoer-2p'*w+w'*R*w%代价函数（误差性能面）
[w0,w1]=meshgrid(-20:.1:20,-20:.1:20);
Jw=(0.5+noisepower)*(w0.^2+w1.^2)+cos(2*pi/N)*w0.*w1+2*w1*sin(2*pi/N)+2;
surf(w0,w1,Jw)
title("第三问的性能曲面")

rm=zeros(1,M);%自相关函数,索引1对应着n=0的情况
rm(1)=0.5+noisepower;
rm(2:end)=cos(2*pi*(1:length(rm)-1)/N);
%这里表示比较蠢 比如r(0)=rm（1），r(1)=rm（2），r(-1)=conj（rm（2）)
%这里的自相关函数显然是实数，因此我们可以把上面的conj去掉

pnegm=-sin(2*pi*(0:M-1)/N);%这里是维纳滤波里面的p
p=pnegm';

R=toeplitz(rm);%协方差矩阵，此东西是N*N维的！

wmin3=R\p;
Jwmin3=expectedpower-p'*wmin3+wmin3'*R*wmin3;hold on;

%%
% 第（4）问
n_axis=-N/4:N-1;%n的向量
n=n_axis(1:N);

phi=2*pi*rand(1);%均匀分布的随机相位
sn_ori=sin(2*pi*n_axis/N+phi);
sn=sn_ori(N/4+1:end);
dn=-2*sn_ori(1:N);

vn=noisepower*randn(1,N/4+N);
un=sn+vn(N/4+1:end);

EIGENVALUE=eig(R)%这里显示出 协方差矩阵的特征值 
MAXMU=2/max(EIGENVALUE);
mu=0.2;
maxerror=.0004;%最大允许误差

w=[0 20]';%确定w0
wnew=zeros(M,1); 
Jwnew=0;
Jw=0;
deltaw=zeros(M,1);
for iternumber=1:100

    Jw=expectedpower-2*p'*w+w'*R*w;
%     Jw=(0.5+noisepower)*(w(1).^2+w(2).^2)+cos(2*pi/N)*w(1).*w(2)+2*w(2)*sin(2*pi/N)+2;
    error=Jw-Jwmin3;

    if(error<maxerror)
        break
    else
        gradient=-2*p+2*R*w;
        deltaw=-mu*gradient/2;
        wnew=w+deltaw;
        Jwnew=expectedpower-2*p'*wnew+wnew'*R*wnew;
        figure(3)
        quiver3(w(1),w(2),Jw,deltaw(1),deltaw(2),Jwnew-Jw);
        figure(4)
        title("学习曲线")
        plot(iternumber,error,"o")
        hold on
        w=wnew;
        Jw=Jwnew;
    end
end
%%

figure(5)
subplot(2,2,1)
plot(n,sn) 
xlabel("n")
ylabel("s(n)")

subplot(2,2,2)
plot(n,dn) 
xlabel("n")
ylabel("d(n)")

subplot(2,2,3)
plot(n,un) 
xlabel("n")
ylabel("u(n)")