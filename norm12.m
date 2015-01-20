function [t]=norm12(x);

t=0;
N=size(x,1);

for i=1:N
    t=t+norm(x(i,:));
end
