N=200;
m=90;
s=20;
d=35;

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

%A=randn(m,N)/sqrt(m);
A=randn(m,N);
AP=gen_AP_matrix(Proj,A,m,N,d);
rank(AP)

[X0,index,x_lambda]=gen_sparse_vec(Frame,M,N,d,s);
X0=X0*10;

lamb_eff=1+norm(ALFA(:,index),'inf')
eff=lamb_eff/s

% x_lambda
% sub_lambda=sub_angle(Frame,index,M)

groups=zeros(N*d,1);
for i=1:N;
    groups((i-1)*d+1:i*d)=i;
end

Y=AP*X0;

opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
X = spg_group(AP,Y,groups,0,opts);
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
Yb=A*Xb0;

opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
Xb = spg_mmv(A,Yb,0,opts);
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


t=0:0.01:1;
%hist(alfa,t);

clear AP;
clear Proj;



