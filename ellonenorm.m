function [a]=ellonenorm(X,N,d);

a=0;

for i=1:N
    a=a+norm(X((i-1)*d+1:i*d));
end
