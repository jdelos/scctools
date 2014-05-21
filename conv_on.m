function [ A_conv_sw_on ] = phase_conv(Acaps,Aloads,Asw)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

A_aux=[1 Asw Acaps Aloads ];
[n m] =size(Asw);

for i=1:m
    Sw_pos = find(A_aux(:,i)==1);
    Sw_neg = find(A_aux(:,i)==-1);
    if ~isempty(Sw_neg)
        A_aux(Sw_pos,:)=A_aux(Sw_pos,:)+A_aux(Sw_neg,:);
        A_aux(Sw_neg,:)=[];
    else
        A_aux(Sw_pos,:)=[];
    end
end

if ~all(A_aux(:,1:m))
    A_conv_sw_on=A_aux(:,m+1:end);
end
    

    
end

