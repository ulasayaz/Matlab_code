function [B]=block(A,d,it,jt);
% returns the i, j th block entry of the
%block matrix of A with block size d OR
% i-th block of a vector X if j==0

if jt==0
    B=A((it-1)*d+1:it*d);
else
    B=A((it-1)*d+1:it*d,(jt-1)*d+1:jt*d);
end
