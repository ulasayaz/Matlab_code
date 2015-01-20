function [lambda]=sub_angle(Frame,index,M);
%%% finds the maximum mutual lambda between subspaces restrited to an index
%%% set of size 's'

s=length(index);

lambda=0;

for i=1:s-1
    for j=i+1:s
        sing=svd(subfr(Frame,M,index(i))'*subfr(Frame,M,index(j)));
        a = max(sing);
        if a > lambda
            lambda=a;
        end
    end
end



