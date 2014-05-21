function [ Bf T ] = fun_loop(A,T)
%%  [Bf T] = fun_loop(A,T)
%Returns the fundamental loops from the matrix A. Twinges can be optionally
% specified in using the vector T.
%
%  Philips Research, Eindhoven,  Netherlands
%  julia.delos@philps.com
%
    
%Remove  zero rows
A(all(A==0,2),:)=[];

if nargin<2
%find a tree of the matrix
[ T ] = build_tree(A);
end

if iscolumn(T)
    T=T.';
end

%number of twinges in the tree
n_twings=length(T);

%get matrix size
[~, m]= size(A);

if all((T == 1:n_twings)==1)
   Tx=1; %Columns keep order
else
    %Get links col index 
    col_order = 1:m;
    col_order = [T col_order(~ismember(col_order,T))];

    %split matrix as A = [ At | Al ] using the transformation matrix Tx
    Tx=zeros(m);
    for i=1:m
        Tx(col_order(i),i)=1;
    end
    %Transform the matrix A to A =[At Al]
    A=A*Tx;
end

%Twings loop matrix Bt=-Al*At^-1;
Bt=-A(:,1:n_twings)\A(:,n_twings+1:end);
Bt=Bt.';
Bf=[Bt eye(size(Bt,1))]/Tx;
        

end

