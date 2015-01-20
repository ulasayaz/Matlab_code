count=500;

C=zeros(count,1);
D=zeros(count,1);

N=12;
    
    for t=1:count
        
        x=randn(N,1);
        dx=zeros(N,1);
        
        x=x/norm(x);
        x=abs(x);
        
        for i=1:N
            dx(i)=discrete(x(i),N);
        end
        
%         norm(dx)
        D(t)=norm(dx);
        
        %error=norm(x-dx)*(1+norm(dx));
        error=norm(x-dx);
        C(t)=error;
    end
    
    c=max(C)
    d=max(D)
    


