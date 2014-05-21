function [ Y ] = D3_matrix_addition(A,weight)
%D3_matrix_addition Third dimesnion matrix addition
%
%D3_matrix_addition(A) Adds the 3-D matrix n x m x p along with the
%   dimension p in oedre to return a n x m Y matrix
%
%D3_matrix_addition(A, weight) Adds the 3-D matrix n x m x p along with the
%   dimension p in oedre to return a n x m Y, matrix saling the 2-D matrix
%   with the weight vector 
%
%  Philips Research, Eindhoven,  Netherlands
%  julia.delos@philps.com
%
    
	
	
[n, m, p] =size(A); %Get matrix size

Y=zeros([n m]); % Initialize the ouput matrix for a 2-D size

if nargin<2
    weight(1:p)=1;
else
    if length(weight)~=p
        error('Matlab:D3_matrix_addition','wrong size of weight paramter')
    end
end

for i=1:p
    Y=weight(i).*A(:,:,i) + Y;
end




end

