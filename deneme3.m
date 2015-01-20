x=randn(20,6);
y=randn(20,7);

x=orth(x);
y=orth(y);

p=x*x';
q=y*y';



s=svd(x'*y);

%%
% calculates lambda for various, d, N, w values.. 

N=8;
d=5;
m=4;

Frame=zeros(d,N);
M=zeros(N,1);

% we fix the dimension of each subspace to 1
for i=1:N
    M(i)=2;
end

Mtot=sum(M);
[Frame,Proj]=gen_frame(M,N,d);

S=min_angle(Frame,M,N);

for j=1:N
plot([Frame(1,j) -Frame(1,j)],[Frame(2,j) -Frame(2,j)]);
hold on;
end

x=gen_sparse_vec(Frame,M,N,d,5);

A=randn(m,N)/sqrt(m);
AP=gen_AP_matrix(Proj,A,m,N,d);

%%
% try group sparsity of SPGL
n=60;
nGroups=30; 
groups=sort(ceil(rand(n,1) * nGroups));
p=randperm(nGroups); p=p(1:3);
idx=ismember(groups,p);

x0=zeros(n,1); x0(idx)= randn(sum(idx),1);


 