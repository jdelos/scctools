function [ A_switched ] = phase_conv(Aelem,Asw)
%PHASE_CONV - Generates the corresponding total converter incidence matrix for a given
% phase defined by the Aelem incidence matrix and the corresponding phase switches Asw
%
%Copyright 2013-2014, Julia Delos, Philips Research 
%	julia.delos@philips.com	
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.

A_aux=[ Asw Aelem ];

[~, m] =size(Asw);

for i=1:m
    Sw_pos = find(A_aux(:,i)==1);
    Sw_neg = find(A_aux(:,i)==-1);
    if isempty(Sw_neg)     
        A_aux(Sw_pos,:)=[];
    elseif isempty(Sw_pos)
        A_aux(Sw_neg,:)=[];
    else
        A_aux(Sw_pos,:)=A_aux(Sw_pos,:)+A_aux(Sw_neg,:);
        A_aux(Sw_neg,:)=[];       
    end
end

if ~all(A_aux(:,1:m))
    A_switched=A_aux(:,m+1:end);
end
    

    
end


