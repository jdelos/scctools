function [ Aconv ] = short_edge( Aconv,edge )
%short_edge Transforms the incidence matrix Aconv that corresponds to 
% the description of a converter introducing a short-circuit in the node corresponding 
% to edge
%
%Copyright 2013-2014, Julia Delos, Philips Research 
%	julia.delos@philips.com	
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.

 Node_pos = find(Aconv(:,edge)==1);
 Node_neg = find(Aconv(:,edge)==-1);
 
 if ~isempty(Node_pos) || ~isempty(Node_neg)
     Aconv([Node_pos Node_neg],:)=[];
 else
     Aconv(Node_pos,:)=Aconv(Node_pos,:)+Aconv(Node_neg,:);
     Aconv(Node_neg,:)=[];
 end
 Aconv(:,edge)=[];

end

