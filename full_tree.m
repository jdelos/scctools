function y = full_tree(A)
%A is the incidence matrix of the tree elements
%of a dirceted Graph.
%
%Full_tree checks if the A 
%is a valid tree matrix 
%
%  Philips Research, Eindhoven,  Netherlands
%  julia.delos@philps.com
%
    
	
    y = false;
    if ~iscolumn(A)
        node_vec = sum(abs(A).');
    else
        node_vec = A';
    end
     y = all(node_vec~=0);
end