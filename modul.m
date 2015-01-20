function [y]=modul(x,l);

d=length(x);
y=zeros(d,1);

for k=1:d
    y(k)=exp(1i*2*pi*l*k/d)*x(k);
end
