function [B]=bernoulli(m,N);

B=zeros(m,N);
for i=1:m
    for j=1:N
        r=rand;
        if r <= 0.5
            B(i,j)=-1;
        else
            B(i,j)=1;
        end
    end
end
            
