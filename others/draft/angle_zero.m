d=10;
N=20;

A=rand(d,N);
B=normc(A);

D=zeros(d,d);

for i=1:d 
    D=D+ B(:,i) * B(:,i)';
end
norm(D)