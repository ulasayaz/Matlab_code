function [X]=gen_compress_vec(Frame,M,N,d,q); 
% this function generates a compressible block vector X with norms
% according to ||x_i||_2 = i^(-1/q)

lambda=0;                                            

X=zeros(N*d,1);

%choosing the sparse index set randomly
% p=randperm(N);
% ind=p(1:s);
% ind=sort(ind);
%ind
%%%

for i=1:N
    X((i-1)*d+1:i*d)=subfr(Frame,M,i)*randn(M(i),1);
    a=norm(X((i-1)*d+1:i*d))*i^(1/q);
    X((i-1)*d+1:i*d)=X((i-1)*d+1:i*d)/a;
end

% for j=1:s-1
%     for k=j+1:s
%         y=X((ind(j)-1)*d+1:ind(j)*d);
%         z=X((ind(k)-1)*d+1:ind(k)*d);
%         y=y/norm(y);
%         z=z/norm(z);
%         a=abs(y'*z);
%         if a > lambda
%             lambda=a;
%         end
%     end
% end

