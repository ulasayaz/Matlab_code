function [B]=normc(B);

d=size(B,2);
for i=1:d
    B(:,i)=B(:,i)/norm(B(:,i));
end
