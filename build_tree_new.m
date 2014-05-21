function [ T ] = build_tree_new( A,initial_edge,excluded_elem)
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

%Eliminate excluded elements
A(:,excluded_elem)=[];

if nargin < 2 %Start edge is the first element
initial_edge=1;
end


%Initialise Tree element vector
%T=initial_edge;

%init edge counter
edge_idx=1; %Edges are explored from the first column 

%Number of nodes
n_node = size(A,1);

%add last row
node_list=[];

A_branch= [A; -sum(A)];
  while n_node < length(node_list)
   %Create auxilary matrix from the remaining brances   
   Aaux=A_branch;
   
   %Get the current element
   current_edge = A(:,edge_idx);
   
   %Find parallel elements
   idx_dup = parallel_elem(current_edge,Aaux);
   if ~isempty(idx_dup)
       %Clean parallel elemens
       Aaux(:,idx_dup)=[];
   end     
   
   %Check ones for the element
   edge_nodes = find( current_edge ~= 0);
   if  ~isempty(edge_nodes)
    %The element is conected to nodes   
    %Check that nodes aren't in the TREE the list
    free_nodes = edge_nodes(~ismember(edge_nodes,node_list));
    if ~isempty(free_node) 
    %The nodes of the current edge are not in the TREE
    %Temporary links matrix current edge 
    Aux=Alinks;                              
    %Make zeros current edge column
    Aux(:,edge_idx) = zeros(n_node,1);
    
       for current_node = free_nodes  %Iterate through the edge nodes
        %Check for edges (1,-1) connected to the current node
        %The first element of the list is used
        edge_cont = find(Aux(1,current_node),1);
        if isempty(edge_cont) 
            %Valid twink edge on node is floating 
            %Node has no contigous elements 
            %Jump back to the previous edge
            next_edge=T(end); 
            %Add current edge to the tree  
            T = [T edge_idx];
            
            
            %The last node hast to be reused to fins a contigous edge
            %This edge won't be found it has been removed from the A_branch
            %So last node us deleted indeed is the other end of this the
            %current edge
            node_list(end)=[];
            
            
            %Add current node to the node list
            node_list = [node_list current_node];
            
            %Update reamining branches matrix
            A_branch=Aaux;
            edge_idx=next_edge;
        else
            %Edge has a contigous node 
            %If the edge is not in the tree add it
            if any(edge_idx==T)
                T=[T edge_idx];
            end
            
            %add node to the list
            node_list = [node_list current_node];
            
            %
            edge_idx=edge_cont;
            
            
            
        end 
       end
    else
        error('build_tree_new:no_solution',...
            'Initial Edge is floating! No contiguous edge!');
    end
   else 
        error('build_tree_new:no_solution',...
            'Initial Edge isn ot connected! Empty row!');
   end
    
  end
end

    

function edges = parallel_elem(edge,A)

n_columns =size(A,2);
%Alocate vector
edges = zeros(n_columns,1);
A=abs(A);

%Make target colum different
for i =1:n_columns
    if i~= edge && all(Ax(:,i)==Ax(:,edge))
     edges(i)=i;
    end
end
%Clean zeor entries
edges(edges==0)=[];     
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

function y = floating_edge(A_tree,edge)
    y = true;
    nodes = find(edge~=0);
    i=1;
    n = length(nodes);
    while (y && i <= n)
      y = all(A_tree(nodes(i),:));
      i=i+1;
    end 
end

function y = next_edge(A,T,i)
    current_edge =  A(:,i);
    Afree = A;
    Afree(:,[T i]);
    
end