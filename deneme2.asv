d=4;
N=10;
x=randn(d,1);
x=x/norm(x);
a=randn;
A=a*x*x';

for i=1:N-1;
    x=randn(d,1);
    x=x/norm(x);
    a=randn;
    B=a*x*x';;
    A=[A B];
end

%%
N=40;
%m=5;
%s=5;
d=[10:2:50];
D=size(d,2);
L=zeros(D,1);

M=zeros(N,1);

% we fix the dimension of each subspace to 1
for i=1:N
    M(i)=1;
end
Mtot=sum(M);


for j=1:D
    Frame=zeros(d(j),Mtot);
    for k=1:N
        r=randn(d(j),1);
        r=r/norm(r);
        Frame(:,k)=r;
    end
    L(j)=cos(min_angle(Frame,M,N,d(j))); 
end

%%

N=15;
m=8;
s=3;
d=10;

M=zeros(N,1);

% we fix the dimension of each subspace to 1
for i=1:N
    M(i)=1;
end
Mtot=sum(M);
Frame=zeros(d,Mtot);
for k=1:N
    r=randn(d,1);
    r=r/norm(r);
    Frame(:,k)=r;
end
X0=gen_sparse_vec(Frame,M,N,d,s);
A=bernoulli(m,N);
Y=A*X0;
opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
X = spg_mmv(A,Y,0,opts);



