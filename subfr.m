function [X]=subfr(f,M,ind); 
%subfr returns the orthonormal matrix consisting of a basis for the 
%subspace W_i

s=sum(M(1:ind-1));
X=f(:,(s+1):(s+M(ind)));