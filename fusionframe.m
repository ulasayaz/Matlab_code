N=150;
%m=5;
%s=5;
d=10;

M=zeros(N,1);

% we fix the dimension of each subspace to 1
for i=1:N
    M(i)=1;
end
Mtot=sum(M);
Frame=zeros(d,Mtot);

% X=randn(d,m(1));
% X=GramS(X);
% 
% Frame(:,1:m(1))=X;
% 
% subfr(Frame,1,m)

teta=5; % fix the minimum angle between subspaces
rad=teta*pi/180;
lambda=cos(rad);

%%% create the subspaces with angles iteratively
% e1=[1,0,0]';
% e2=[0,1,0]';
% Frame(:,1)=e1;
% Frame(:,2)=e2;
% 
% %rotation matrix about e1 vector
% R=eye(3)*cos(rad)+sin(rad)*[0,0,0;0,0,-1;0,1,0] + (1-cos(rad))*[1,0,0;0,0,0;0,0,0];R^2;
% 
% count=1;
% for j=1:N-1
%     c1=count+M(j);
%     count=c1;
%     Frame(:,c1)=e1;
%     Frame(:,c1+1)=R^j*e2;
% end

% for i=1:N
%     r=randn(d,1);
%     r=r/norm(r);
%     Frame(:,i)=r;
% end

%%%%%%
%subspace(subfr(Frame,M,1),subfr(Frame,M,2))
 
%we try now 
% A=bernoulli(m,N);
% X0=gen_sparse_vec(Frame,M,N,d,s);

% Y=A*X0;
% opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
% X = spg_mmv(A,Y,0,opts);            % Choose sigma = 0
    
    
% Test for fixed s and N, fixed random vector X, with varying m, 
% for each m, try 100 times, and record the success rate. 

%X0=gen_sparse_vec(Frame,M,N,d,s);

S=[5:1:12];
sizeS=size(S,2);
L=zeros(sizeS,1);
Meas=zeros(sizeS,1);

try_m=[5:2:90];
num=size(try_m,2);

trial=50;



for i=1:sizeS
    i
    s=S(i);
    %create a random frame
    for k=1:N
        r=randn(d,1);
        r=r/norm(r);
        Frame(:,k)=r;
    end  
    L(i)=cos(min_angle(Frame,M,N,d)); %calculate its lambda
    X0=gen_sparse_vec(Frame,M,N,d,s); %generate random sparse vector
    %find number of measurements with 90% success
    for l=1:num
        count=0;
        for j=1:trial
            A=bernoulli(try_m(l),N);
            Y=A*X0;
            opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
            X = spg_mmv(A,Y,0,opts);            % Choose sigma = 0
            if norm(X-X0,'fro') < 1e-3
                count=count+1;
            end
        end
        if count/trial > 0.9
            Meas(i)=try_m(l);
            break;
        end
    end %finished calculating measurement number
end

        
        



    
    
    