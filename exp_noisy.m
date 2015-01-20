N=200;
m=50;
s=30;
d=35;

M=zeros(N,1);

% we fix the dimension of each subspace to 1
for i=1:N
    M(i)=5;
end

clear AP;
clear Proj;

[Frame,Proj]=gen_frame(M,N,d);

[alfa,ALFA]=min_angle(Frame,M,N);
lambda=max(alfa)
mean=sum(alfa)/length(alfa)

A=randn(m,N)/sqrt(m);
AP=gen_AP_matrix(Proj,A,m,N,d);

[X0,index,x_lambda]=gen_sparse_vec(Frame,M,N,d,s);
X0=X0/norm(X0);
lamb_eff=1+norm(ALFA(:,index),'inf')

erb=randn(m,d);  % 
erb=erb/norm(erb,'fro');
f=erb';
err=f(:);
sigma=0;

q=0.4;
X0=gen_compress_vec(Frame,M,N,d,q);
X0=X0/norm(X0);

% x_lambda
% sub_lambda=sub_angle(Frame,index,M)

groups=zeros(N*d,1);
for i=1:N;
    groups((i-1)*d+1:i*d)=i;
end

Y=AP*X0+sigma*err;

opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
X = spg_group(AP,Y,groups,sigma,opts);
r=abs(X-X0);
norm(X-X0)

%%% block sparse experiment without frame structure
% p=randperm(N); p=p(1:s);
% Xb0=zeros(N,d);
% Xb0(p,:)=10*randn(s,d);
% 
% Yb=A*Xb0;

%make vector X0 from FF into a block form
Xb0=zeros(N,d);
for i=1:N
    Xb0(i,:)=X0((i-1)*d+1:i*d);
end
Yb=A*Xb0+sigma*erb;

opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
Xb = spg_mmv(A,Yb,sigma,opts);
norm(Xb-Xb0,'fro')

% cvx_begin quiet
% variable XbC(N,d) 
% %minimize( sum( sqrt(sum(Xb.*Xb.'',2))) )
% minimize(norm12(XbC))
% subject to
% %norm(A*Xb-Yb,'fro') <= 1e-4
% A*XbC==Yb;
% %       norm12(Q)
% for a=1:N
%     % norm(subfr(Frame,M,a)*subfr(Frame,M,a)'*XbC(a,:).'-XbC(a,:).') <= 1e-3;
%     subfr(Frame,M,a)*subfr(Frame,M,a)'*(XbC(a,:).') == XbC(a,:).';
% end
% cvx_end
% % CVX program END
% 
% norm(XbC-Xb0,'fro')


%t=0:0.01:1;
%hist(alfa,t);

clear AP;
clear Proj;



