function [y]=trans(x,k);

d=length(x);
y=zeros(d,1);

for t=1:d
    ind=mod(t-k,d);
    if ind==0
        ind=d;
    end
    y(t)=x(ind);
end
