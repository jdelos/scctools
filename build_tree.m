function [ T ] = build_tree( A,initial_edge,excluded_elem)
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
%  Philips Research, Eindhoven,  Netherlands
%  julia.delos@philps.com
%


if ~isnumeric(A)
    SymList = symvar(A);
    A = subs(A,SymList,ones(1,length(SymList)));
end


n_branches = size(A,2);

switch nargin
    case 1 
        for i=1:n_branches
            if any(A(:,i))
                initial_edge=i;
                break
            end
        end
        %Start edge is the first element
    case 3    
        %Eliminate excluded elements
        A(:,excluded_elem)=[];
end

if ~any(A(:,initial_edge))
    error('build_tree:ArgError','Initial Edge is not a valid branch!')
end

%Function needs the complete incidence matrix
A = [A; -sum(A,1)];
%Number of nodes
n_nodes = size(A,1);


%Initialize nodes lst
int_branc = A(:,initial_edge);
%Look for the nodes of the current branch nodes are where (1,-1)
edge_nodes = find( int_branc ~= 0);


%Start recursive function of tree builder
[~,tree,out_node]= tree_builder(A,initial_edge,0,n_nodes,edge_nodes);

if out_node == -1
    T=sort(tree);
else
    T=-1;
end


end


function [A_branch,tree,in_node,lst_nodes]= tree_builder(A_branch,tree,in_node,n_nodes,lst_nodes)


if in_node == -1 || length(tree) == (n_nodes - 1)
    %tree is completed
    in_node = -1;
    return
else
    %last value in tree vector is actual branch
    current_branc = A_branch(:,tree(end));
    %Find parallel elements
   idx_dup = parallel_elem(tree(end),A_branch);
   if ~isempty(idx_dup)
       %Clean parallel elemens
       A_branch(:,idx_dup)=A_branch(:,idx_dup)-A_branch(:,idx_dup);
   end 
   %Look for the nodes of the current branch nodes are where (1,-1)
   edge_nodes = find( current_branc ~= 0);
   %Pick the oposite node that brought to this branch
   edge_nodes = edge_nodes(edge_nodes~=in_node)';
   
   %Since the current branch is part of the tree is canceled making it the
   %branch 0
   A_branch(:,tree(end))=zeros(n_nodes,1);
   
   %The initial node (0)  will bring 2 options, both must be explored 
   for current_node = edge_nodes
   
   %Look for other branches connected to this node
   %Check for edges (1,-1) connected to the current node
   
   
   %Get the list of adjacent branches 
   adj_branches = find(A_branch(current_node,:)~=0);
   while ~isempty(adj_branches)
    %There are adjacent branches
    %Iterate through them
    branch = adj_branches(1);
    
    adj_branches(1)=[];
    
     %Check if the branch closes the loop 
     branch_nodes = find( A_branch(:,branch)~=0);
     if  ~loop_branch(branch_nodes,lst_nodes) 
      %Add it to the tree and jump to it
      tree = [tree branch]; %#ok<*AGROW>
      
      %Add nodes to the node list
      lst_nodes = unique([lst_nodes; branch_nodes]);
      
      %Check if it completes the tree
      if length(tree) == (n_nodes - 1)
        %tree is completed 
        in_node = -1;
        return
       else
        %Still branches to grow
        %explore it
       [A_branch,tree,in_node,lst_nodes] = ...
           tree_builder(A_branch,tree,current_node,n_nodes,lst_nodes);
         if in_node == -1
           %Tree is completed return
           return
         end
      %Update list of adjacent branches, parallel branches will not be
      %anymore in the matrix
      adj_branches = find(A_branch(current_node,:)~=0);    %#ok<NASGU>
      end
     end
    end
    %All the branches of this node are explored this node is empty\
   end
   %Empty node nothing to explore return
   in_node = current_node;
   end
 
end
   


function edges = parallel_elem(edge,A)

n_columns =size(A,2);
%Alocate vector
edges = zeros(n_columns,1);
Ax=abs(A);

%Make target colum different
for i =1:n_columns
    if i~= edge && all(Ax(:,i)==Ax(:,edge))
     edges(i)=i;
    end
end
%Clean zeor entries
edges(edges==0)=[];     
end

% function y = loop_tree(A_tree,edge)
%     y = true;
%     nodes = find(edge~=0);
%     i=1;
%     n = length(nodes);
%     while (y && i <= n)
%       y = all(A_tree(nodes(i),:));
%       i=i+1;
%     end 
% end
% 
% function y = floating_edge(A_tree,edge)
%     y = true;
%     nodes = find(edge~=0);
%     i=1;
%     n = length(nodes);
%     while (y && i <= n)
%       y = all(A_tree(nodes(i),:));
%       i=i+1;
%     end 
% end
% 
% function y = next_edge(A,T,i)
%     current_edge =  A(:,i);
%     Afree = A;
%     Afree(:,[T i]);
%     
% end

function  y = loop_branch(branch_nodes,lst_nodes)
% In case that both nodes that connects the branch are in the lst of nodes
% this is a loop branch 
y = sum(ismember(lst_nodes,branch_nodes)) == 2;
end
