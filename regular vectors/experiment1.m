count=500;

C=zeros(count,1);



N=10;

for t=1:count
    a=randn(N,1);
    b=randn(N,1);
    
    a=(2/3)*a/norm(a);
    b=0.48*b/norm(b);
    
    a=abs(a);
    b=abs(b);
    
    C(t)=dot(a,b);
end

c=min(C)