count=500;
C=zeros(count,1);
L=[2:10:302];

n=length(L)
ERR=zeros(n,1);

for j=1:n
    j
    N=L(j);
    
    for t=1:count
        
        x=randn(N,1);
        dx=zeros(N,1);
        
        x=x/norm(x);
        x=abs(x);
        
        for i=1:N
            dx(i)=discrete(x(i),N);
        end
        
        %         norm(dx)
        %         norm(x-dx)
        
        %error=norm(x-dx)*(1+norm(dx));
        error=norm(dx);
        
        if norm(x-dx) > 0.3
            C(t)=norm(dx);
        else
            C(t)=0;
        end
        
        %C(t)=error;
    end
    
    ERR(j)=max(C);
    
end

ERR
plot(L,ERR)