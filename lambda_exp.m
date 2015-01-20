% m=6;
% 
% A = rand(m);
% A = triu(A) + triu(A,1)';
% for i=1:m
%     A(i,i)=0;
% end
% 
% S=0;
% 
% for i=1:m
%     S=S+max(A(i,:));
% end
% 
% S
% mean=sum(A(:))/(m-1)

N=200;
s=20;
d=10;


M=zeros(N,1);
% we fix the dimension of each subspace to 1
for i=1:N
    M(i)=1;
end

clear AP;
clear Proj;

[Frame,Proj]=gen_frame(M,N,d);

p=randperm(N);
ind=p(1:s);
ind=sort(ind);

PCOL=zeros(N*d,d);

for j=1:N
    b=bernoulli(1,1);
    %b=randn;
    PCOL((j-1)*d+1:j*d,:)= b * block(Proj,d,j,j);
end

[alfa,ALFA]=min_angle(Frame,M,N);
lambda=max(alfa)
mean=sum(alfa)/length(alfa)

BETA=triu(ALFA); %upper triangular part of ALFA


x=norm(PCOL)

c1=sqrt(1+lambda*N)
c2=sqrt(1+mean*N)
c3=sqrt(1+norm(ALFA))
c4=sqrt(1+norm(ALFA,1))

e=1+norm(ALFA(:,ind),'inf')


for i=1:s  % we set diangonals to 0 in submatrix "S x S"
    ALFA(ind(i),ind(i))=0;
end

B=ALFA(:,ind); % submatrix of set [N] x S
C=ALFA(ind,ind); % submatrix of set S x S
lamb=max(B(:));
ave=sum(B(:))/length(B(:));
d1=1+norm(C,'inf')
d2=1+norm(B,'inf')
d3=1+mean*s   % that was calculated during the frame experiments
d4=1+lambda*s % that was theoretical bound before
s


%%

N=100;
d=8;

A=zeros(N,N);

for i=1:N-1
    for j=i+1:N
        a=rand;
        A(i,j)=a;
        A(j,i)=a;
    end
end

tot=N*(N-1);
mean=sum(A(:))/tot;

d1=mean*N
d2=norm(A)
d3=norm(A,1)
d4=max(A(:))*N



