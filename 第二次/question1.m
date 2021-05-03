noisepower=0.4;%��������
m=19;%mΪ��,����2*����ͷ
N=4*m;%�ܳ�ͷ����
n_axis=-N/4:N-1;%n������
n=n_axis(1:N);

phi=2*pi*rand(1);%���ȷֲ��������λ
sn_ori=sin(2*pi*n_axis/N+phi);
sn=sn_ori(N/4+1:end);
dn=-2*sn_ori(1:N);

rm=zeros(1,N/4+N);%����غ���,����1��Ӧ��n=0�����
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