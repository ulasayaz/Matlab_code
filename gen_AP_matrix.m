function [X]=gen_AP_matrix(Proj,A,m,N,d); 
% this function generates the block matrix AP with respect
% to the Proj matrix and
% coefficient matrix A

X=zeros(m*d,N*d);

for i=1:m
    for j=1:N
        X((i-1)*d+1:i*d,(j-1)*d+1:j*d)=A(i,j)*Proj((j-1)*d+1:j*d,(j-1)*d+1:j*d);
    end
end
      
    
