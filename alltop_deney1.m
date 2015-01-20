% in this experiment for ALLTOP frames we fix N and d, so lambda is fixed
% then we vary 's' and find the number of measurements 'm'
% that guarantees recovery with %98 for both FF and block sparsity methods.

% we write measurement numbers into a 2-row matrix called 'Meas'

d=13;
N=d^2;
Sparse=[5:5:35];
L=length(Sparse)

%sub_lambda=zeros(1,L);

try_m=[2:2:150];
num=length(try_m)

Meas=zeros(2,L);
trials=50;

M=zeros(N,1);
% we fix the dimension of each subspace to 1
for k=1:N
    M(k)=1;
end

clear AP;
clear Proj;

%%%%% create alltop vector
x=zeros(d,1);

for k=1:d
    x(k)=exp(1i*2*pi*k^3/d)/sqrt(d);
end
%%%%%%  alltop end

%%%%%%%%%  generate frame
Frame=zeros(d,N);
Proj=zeros(d*N,d*N);
for l=1:d
    for k=1:d
        loc=(l-1)*d+k;
        y=modul(trans(x,k),l);
        Frame(:,loc)=y;
        Proj((loc-1)*d+1:loc*d,(loc-1)*d+1:loc*d)=y*y';
    end
end
%%%%%%%%%  end frame

alfa=min_angle(Frame,M,N);
lambda=max(alfa)
mean=sum(alfa)/length(alfa)

groups=zeros(N*d,1); %%%%  define 'groups' vector
for a=1:N
    groups((a-1)*d+1:a*d)=a;
end

for ind=1:L
    ind
    s=Sparse(ind);    
    %%%%%%  generate sparse vector
    X0=zeros(N*d,1);
    p=randperm(N);
    index=p(1:s);
    index=sort(index);
    r=complex(randn(s,1),randn(s,1));
    for a=1:s       
        X0((index(a)-1)*d+1:index(a)*d)=subfr(Frame,M,index(a))*r(a);
    end
    X0=X0*10;
    %%%%%  end sparse vector
    
    %sub_lambda(ind)=sub_angle(Frame,index,M);
    
    Xb0=zeros(N,d); %block sparsity vector
    for a=1:N
        Xb0(a,:)=X0((a-1)*d+1:a*d);
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
            %%%%%  generate random AP matrix
            A=randn(m,N);
            AP=zeros(m*d,N*d);
            for a=1:m
                for b=1:N
                    AP((a-1)*d+1:a*d,(b-1)*d+1:b*d)=A(a,b)*Proj((b-1)*d+1:b*d,(b-1)*d+1:b*d);
                end
            end
            %%%%%%  end AP matrix
            
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
            m
            break;
        end
    end %found 'm' for that 's'
    'finished FF'
    
    % find 'm' for block sparse setup
%     for k=j:num
%         k
%         m=try_m(k);
%         count=0;
%         for t=1:trials
%             %A=randn(m,N)/sqrt(m);
%             if count < t-3
%                 break;
%             end
%             A=randn(m,N);
%             Yb=A*Xb0;
%             
%             opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
%             Xb = spg_mmv(A,Yb,0,opts);
%             
%             if norm(Xb-Xb0,'fro') < 1e-4
%                 count=count+1;
%             end
%         end
%         if count/trials >= 0.95
%             Meas(2,ind)=m; % write m for block sparse
%             break;
%             m
%         end
%     end %found 'm' for that 's'
%     'finished block'
end





