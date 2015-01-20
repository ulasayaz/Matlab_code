function [Frame,Proj]=gen_frame(M,N,d); 
% this function generates a random frame for R^d
% with dimensions in vector M, each
% subframe is an orthonormal basis for W_i

Mtot=sum(M);
Frame=zeros(d,Mtot);
%Proj=[];
Proj=zeros(d*N,d*N);

index=1;
for j=1:N
    x=randn(d,M(j));
    x=orth(x);
    Frame(:,index:index+M(j)-1)=x;
    %Proj=blkdiag(Proj,x*x');
    Proj((j-1)*d+1:j*d,(j-1)*d+1:j*d)=x*x';
    index=index+M(j);
end
    