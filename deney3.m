% in this experiment we investigate the effect of the change in lambda on
% the number of measurements needed. we fix N and s. we vary 'omega' in order to
% get different lambda.

N=180;
s=18;
d=35;

omega=[1:10]; % dimesion of the subspaces
L=length(omega)

lambda=zeros(1,L);
%sub_lambda=zeros(1,L);
mean=zeros(1,L);
lamb_eff=zeros(1,L);

try_m=[2:2:150];
num=length(try_m)

Meas=zeros(1,L);
trials=50; % be careful about this setting

groups=zeros(N*d,1); 
for i=1:N
    groups((i-1)*d+1:i*d)=i;
end

clear AP;
clear Proj;

thres=1;
for ind=1:L
    ind
    clear Proj;
    M=zeros(N,1);
    for i=1:N
        M(i)=omega(ind);
    end
    
    [Frame,Proj]=gen_frame(M,N,d);
    
    [alfa,ALFA]=min_angle(Frame,M,N);
    lambda(ind)=max(alfa);
    mean(ind)=sum(alfa)/length(alfa);
    
    [X0,index,x_lambda]=gen_sparse_vec(Frame,M,N,d,s); % FF sparse vector
    X0=X0*10;
    lamb_eff(ind)=1+norm(ALFA(:,index),'inf');
    
    % sub_lambda(ind)=sub_angle(Frame,index,M);
    
    % find 'm' for FF setup
    for j=thres:num
        j
        m=try_m(j);
        count=0;
        for t=1:trials
            %A=randn(m,N)/sqrt(m);
            if count < t-2
                break;
            end
            A=randn(m,N);
            AP=gen_AP_matrix(Proj,A,m,N,d);
            
            Y=AP*X0;
            
            opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
            X = spg_group(AP,Y,groups,0,opts);
            
            clear AP;
            %norm(X-X0)
            if norm(X-X0) < 1e-4
                count=count+1;
            end
        end
        if count/trials >= 0.97
            Meas(ind)=m; %write m for FF
            thres=j; % avoids extra calculations for 'm'
            m
            break;
        end
    end %found 'm' for that 'omega'
    'finished FF'
end