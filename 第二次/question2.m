%%
resov=[1 -0.975 0.95;-0.975 1.95 0;0.95 -0.975 1];
inv(resov)*[0.0731;0;0]
%%
a1=0.975;
a2=-0.95;
sigmav=0.0731;
%%
%我们将层级滤波器转化为 传递函数,可得:
%H(z)=1/(1-a1/z-a2/z^2)
% 第二问
N=512;%长度
vn=sqrt(sigmav)*randn(1,N);
NUM=1;DEN=[1 -a1 -a2];
un=filter(NUM,DEN,vn);
plot(un);
variationofunis=var(un)
%%
% 第三问
