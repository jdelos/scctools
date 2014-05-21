function [ T ] = build_tree( A,initial_edge)
%T = build_tree  Finds a tree from the incidence matrix A  
% build_tree(A) from the incidence matrix A including the 
% first edge of the matrix as the frist element of the tree
%                 
% build_tree(A,initial_edge) from the incidence matrix A 
% including the initial_edge as the first element of the tree
%
% The function returns T vector containing the indexes of the 
% edges corresponding to tree
%               
                   

if nargin < 2 %Start edge is the first element
initial_edge=1;
end

%Initialise Tree element vector
T=initial_edge;

%init edge counter
i=1; %Edges are explored from the first column 

[n m] = size(A);
while ~full_tree(A(:,T))  
     
     if i~=initial_edge && ~loop_tree(A(:,T),A(:,i))  
        T=[T i];
     end
    i=i+1;  
end


end



function y = loop_tree(A_tree,edge)
    y = true;
    nodes = find(edge~=0);
    i=1;
    n = length(nodes);
    while (y && i <= n)
      y = all(A_tree(nodes(i),:));
      i=i+1;
    end
    
end

