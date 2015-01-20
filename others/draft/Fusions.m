function [weightedF, localDualF, weightedLocalDualF, Sw, invSw, globalDual, fusionDual] = FusionFrame(K,F,v)

%
%
%   General: Given a collection of weighted subspaces defined by a collection of local frames 
%            and associated weights, this program computes various data important for fusion 
%            frame applications.
%
%   Syntax:  [weightedF, localDualF, weightedLocalDualF, Sw, invSw, globalDual, fusionDual] 
%             = FusionFrame(K,F,v).
%  
%   Input:   K=[K1,K2,...,KN]   is the collection of numbers of frame elements Ki spanning Wi.
%            F=[F1 F2... FN]    is the collection of all frame sequences in Column-by-Column format.   
%                              Fi consists of the frame elements spanning subspace Wi. 
%            v=[v1, v2,..., vN] is the collection of weights for each subspace.
%
%   Output:  weightedF          is the weighted collection of frame sequences.
%            localDualF         is the collection of local dual frame sequences.
%            weightedLocalDualF is the weighted collection of local dual frame sequences.
%            Sw                 is the fusion frame operator.
%            invSw              is the (pseudo) inverse of the fusion frame operator.
%            globalDual         is the global dual frame to $\{v_i f_{ij}\}_{i,j}$.
%            fusionDual         is the fused global dual frame.
%          
%   Version: 1.0
%   Date:    12/9/2006
% 


N=length(K); % extracting the number of subspaces.
[M,L]=size(F);
weightedF=zeros(M,L);
localDualF=[];
weightedLocalDualF=[];
Sw=zeros(M,M);  % M=dimenion of the vectors.

counter=0;
for n=1:N
    c1=counter+1;   % set the beginning column number of the n_th subspace frame elements.
    c2=counter+K(n);  % set the ending column number of the n_th subspace frame elements.
    counter =c2;  % set the counter at the end of n_th subspace elements, and preprae for the n+1_th subspace.
    Fn=F(:, c1:c2);  % extracting the frame elements spanning Wn.  
    weightedF(:, c1:c2)=v(n)*Fn;     % generating the weighted frame elements of Wn.  
    %%%   calculating local dual frame   %%%%
            localDualFn=(pinv(Fn))'; 
            localDualF=[localDualF, localDualFn];  %putting local duals in a matrix.  
            weightedLocalDualF=[weightedLocalDualF, v(n)*localDualFn];  %putting weighted local duals in a matrix.  
    %%%   calculating the n_th components of the fusion frame operator  %%%%%
            Swn=(v(n)^2)*localDualFn*Fn';    
            Sw=Sw+Swn; 
end 
invSw=pinv(Sw);   %  pinv is used in case that the span{Wn} is not the whole space.  

fusionDual = invSw*weightedLocalDualF; 
globalDual = pinv(weightedF)';
