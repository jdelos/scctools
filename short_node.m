function [ Aconv ] = short_node( Aconv,edge )
%short_node Transforms the incidence matrix Aconv that corresponds to 
% the description of a converter introducing a short-circuit in the node corresponding 
% to edge
%
%Copyright 2013-2014, Julia Delos, Philips Research 
%	julia.delos@philips.com	
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.

 Node_pos = find(A_aux(:,i)==1);
 Node_neg = find(A_aux(:,i)==-1);
 
 if ~isempty(Sw_neg) || ~isempty(Sw_pos)
     Aconv([Sw_neg Sw_pos],:)=[];
 else
     Aconv(Node_pos,:)=Aconv(Node_pos,:)+Aconv(Node_neg,:);
     Aconv(Node_neg,:)=[];
 end
 

end

