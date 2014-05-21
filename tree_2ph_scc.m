function [ Tph1 Tph2 ] = tree_2ph_scc(A_ph1 , A_ph2, n_caps,mode )
%tree_2ph_scc Provides the trees of a well proposed SCC
%   
%The tree branches are composed using only the capacitors and the input 
%supply. The load sources are always links of the graph.
%The capacitors contained in the tree of the first phase are discarded for 
%the tree of the second phase.  
%
%Copyright 2013-2014, Julia Delos, Philips Research 
%	julia.delos@philips.com	
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.



if mode
    %Boost mode operation input source must be avoided to compute the trees
    Tph1 =  build_tree( A_ph1(:,2:end))+1;
    Tph2 =  build_tree( A_ph2(:,2:end))+1;
else
    %Buck mode operation ouput loads must be out of the tree generator 
Tph1 = build_tree( A_ph1(:,1:n_caps+1));
Tph2 =  build_tree( A_ph2(:,1:n_caps+1));
end
%chec if it is a valid tree
if ~full_tree(A_ph2(:,Tph2))
    error('MATLAB:badInputMatrix','Non well proposed converter');
end

end

