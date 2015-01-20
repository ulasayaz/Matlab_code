% in this experiment we fix N and d, so lambda is fixed
% then we vary 's' and find the number of measurements 'm'
% that guarantees recovery with %96 for both FF and block sparsity methods.

% we write measurement numbers into a 2-row matrix called 'Meas'

N=200;
d=7;
Sparse=[5:5:35];
L=length(Sparse)

%sub_lambda=zeros(1,L);
lamb_eff=zeros(1,L);

try_m=[2:2:130];
num=length(try_m)

Meas=zeros(2,L);
trials=50;

M=zeros(N,1);
% we fix the dimension of each subspace to 1
for i=1:N
    M(i)=2;
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

for ind=1:L
    ind
    s=Sparse(ind);
    [X0,index,x_lambda]=gen_sparse_vec(Frame,M,N,d,s); % FF sparse vector
    X0=X0*10;
    lamb_eff(ind)=1+norm(ALFA(:,index),'inf');
    
    %sub_lambda(ind)=sub_angle(Frame,index,M);
    
    Xb0=zeros(N,d); %block sparsity vector
    for i=1:N
        Xb0(i,:)=X0((i-1)*d+1:i*d);
    end
    
    % find 'm' for FF setup
    for j=1:num
        j
        m=try_m(j);
        count=0;
        for t=1:trials
            %A=randn(m,N)/sqrt(m);
            if count < t-3
                break;
            end
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
        if count/trials >= 0.95 % we want at least 48 counts over 50
            Meas(1,ind)=m; %write m for FF
            break;
        end
    end %found 'm' for that 's'
    'finished FF'
    
    % find 'm' for block sparse setup
    for k=j:num
        k
        m=try_m(k);
        count=0;
        for t=1:trials
            %A=randn(m,N)/sqrt(m);
            if count < t-3
                break;
            end
            A=randn(m,N);
            Yb=A*Xb0;
            
            opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
            Xb = spg_mmv(A,Yb,0,opts);
            
            if norm(Xb-Xb0,'fro') < 1e-4
                count=count+1;
            end
        end
        if count/trials >= 0.95
            Meas(2,ind)=m; % write m for block sparse
            break;
        end
    end %found 'm' for that 's'
    'finished block'
end





