function [ Mo ] = matrix_2_expon (M)
% MATRIX_2_EXPON Returns all the multipliers for the rows of matrix M with dimension nXp in 
% a 3-dimensional matrix Mo
%
% Each row in M  produces a pXp matrix with all the possible multipliers of that vector
% The third dimension in Mo corresponds to the row index  
%
%Copyright 2013-2014, Julia Delos, Philips Research 
%	julia.delos@philips.com	
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.

[n m] =size(M);
if ~isobject(M)
    Mo=zeros(m,m,n);
else
    Mo=sym(zeros(m,m,n));
end

for i=1:n
    Mo(:,:,i)=(M(i,:).'*M(i,:));
end



end

