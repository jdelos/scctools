function [ Q T ] = fun_cutset( A,T )
%fun_cutset Compute the fundamental cutset of the incidence matrix A
%according using the tree specified with the vector T. It contains the
%corresponding column index to the edges forming the tree 
%
% The fundamental cut-set matrix Q can be computed as follows
%
%                   Q= [ I | At\Al]
%where At is the part of the incindence matrix A correponding to the twing
%branches, and Al to the link branches
%
%  Philips Research, Eindhoven,  Netherlands
%  julia.delos@philps.com
%
    

if nargin<2
    T = build_tree(A,1,[]);
end

if iscolumn(T)
    T=T.';
end
%number of twings in the tree
n_twings=length(T);

%get matrix size
[n m]= size(A);

col_order=1:m;
col_order=[T, col_order(~ismember(col_order,T)) ];

%Initialize transformation matrix
Tx=zeros(m);
for i=1:m
    Tx(col_order(i),i)=1;
end
%Transform the matrix A to A =[At Al]
A=A*Tx;

%Compute links cut-set matrix
Ql=A(:,1:n_twings)\A(:,(n_twings+1):end);

Q=[eye(size(Ql,1)) Ql];

%reorder the matrix
Q=Q/Tx;


end

