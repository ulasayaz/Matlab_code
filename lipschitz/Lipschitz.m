m=20;
N=60;
d=20;

trial=100;
ctr=0;

for i=1:trial;
    i
    % we fix the dimension of each subspace to 1
    for i=1:N
        M(i)=10;
    end
    
    [Frame,Proj]=gen_frame(M,N,d);
    
    A=randn(m,N).^2-randn(m,N);
    AP=gen_AP_matrix(Proj,A,m,N,d);
    
    a=norm(A);
    b=norm(AP);
    
    if a < b
        fprintf('Alert');
        ctr=ctr+1;
    end
   
end
ctr