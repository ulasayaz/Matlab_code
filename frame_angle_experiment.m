N=100;

d=[5:5:70];

S=size(d,2);
lambda=zeros(S,1);

M=zeros(N,1);

% we fix the dimension of each subspace to 1
for i=1:N
    M(i)=1;
end
Mtot=sum(M);

trial=10;
for i=1:S
    for t=1:trial
        Frame=gen_frame(d(i),M,N);
        lambda(i)=lambda(i)+min_angle(Frame,M,N);
    end
    lambda(i)=lambda(i)/trial;
end

