function [b]=discrete(a,N);

%N=100;

k=-log2(a);

kc=ceil(k);
kf=floor(k);

b=2^(-kc);

% if (a-2^(-kc)) >= (2^(-kf)-a)
%     b=2^(-kf);
% else
%     b=2^(-kc);
% end

if k > ceil(log2(3*N))
    b= 2^(-ceil(log2(3*N)));    
end