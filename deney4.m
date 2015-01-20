% in this experiment we fix N,s and d, so lambda is fixed
% then we plot success rate vs number of measurements 'm' for two cases of
% block sparsity and fusion frame..
% it is the same setup in deney1..

N=200;
d=3;
s=20;

%sub_lambda=zeros(1,L);

try_m=[2:2:100];
num=length(try_m)

Prob=zeros(2,num);
trials=100;

M=zeros(N,1);
% we fix the dimension of each subspace to 1
for i=1:N
    M(i)=1;
end

clear AP;
clear Proj;

[Frame,Proj]=gen_frame(M,N,d);

[alfa,ALFA]=min_angle(Frame,M,N);
lambda=max(alfa)
mean=sum(alfa)/length(alfa)

groups=zeros(N*d,1);
for i=1:N
    groups((i-1)*d+1:i*d)=i;
end

[X0,index,x_lambda]=gen_sparse_vec(Frame,M,N,d,s); % FF sparse vector
X0=X0*10;

lamb_eff=1+norm(ALFA(:,index),'inf');

Xb0=zeros(N,d); %block sparsity vector
for i=1:N
    Xb0(i,:)=X0((i-1)*d+1:i*d);
end

for ind=1:num
    ind   
    m=try_m(ind);
    
    % find 'success rate' for FF setup
    count=0;
    for t=1:trials
        %A=randn(m,N)/sqrt(m);
        
        A=randn(m,N);
        AP=gen_AP_matrix(Proj,A,m,N,d);
        
        Y=AP*X0;
        
        opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
        X = spg_group(AP,Y,groups,0,opts);
        
        %clear AP;
        %norm(X-X0)
        if norm(X-X0) < 1e-4
            count=count+1;
        end
    end
    'finished FF'
    Prob(1,ind)=count/trials; %write 'success rate' for FF
    %     if count == trials
    %         break;
    %     end
    
    count=0;
    % find 'success rate' for block sparse setup
    for t=1:trials
        %A=randn(m,N)/sqrt(m);
        
        A=randn(m,N);
        Yb=A*Xb0;
        
        opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
        Xb = spg_mmv(A,Yb,0,opts);
        
        if norm(Xb-Xb0,'fro') < 1e-4
            count=count+1;
        end
    end
    'finished block'
    Prob(2,ind)=count/trials; %write 'success rate' for block
end

plot(try_m,Prob(1,:));
hold
plot(try_m,Prob(2,:));


