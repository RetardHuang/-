noisepower=0.4;%噪声功率
m=19;%m为几,就是2*几抽头
N=4*m;%总抽头数量
n_axis=-N/4:N-1;%n的向量
n=n_axis(1:N);

phi=2*pi*rand(1);%均匀分布的随机相位
sn_ori=sin(2*pi*n_axis/N+phi);
sn=sn_ori(N/4+1:end);
dn=-2*sn_ori(1:N);

rm=zeros(1,N/4+N);%自相关函数,索引1对应着n=0的情况
rm(1)=0.5+noisepower;
rm(2:end)=cos(2*pi*(1:length(rm)-1)/N);


vn=noisepower*randn(1,N/4+N);
un=sn+vn(N/4+1:end);

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