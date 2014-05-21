function [ A_caps A_sw1 A_sw2 ] = ladder_matrix(n_stages, in_cap)
%DICKSON_MATRIX Creates the incindence matrixes for an N_STAGES (#capacitors)
% Dikson Ladder converter
% 
% Retunr values 
%   A_caps -> Capacitor Incidnece matrix Inculuding the voltage supply
%             connected to the first node
%   A_sw1  -> Phase 1 ON Switches 
%   A_sw2  -> Phase 2 ON switches
%
%The firs row of the incidence matrix corresponds to the input supply node
%The Phase 1 is the one that drawns power from the load
%
%Copyright 2013-2014, Julia Delos, Philips Research 
%	julia.delos@philips.com	
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.

%Initialize matrix
A_caps = zeros(n_stages+2,n_stages);
A_sw = zeros(n_stages+2,n_stages+2);

%Loop Caps matrix
for i=1:n_stages-1
%    if i<n_stages
     A_caps([i+1 i+3],i)=[1 -1];
%    else
%     A_caps(i+1,i) = 1;  
%    end  
%   
end
%loop Switches Matrix

for i=1:n_stages+1
%    if i<n_stages+2
     A_sw(i:i+1,i) = [1; -1]; 
%    else
%     A_sw(i,i) = 1;    
%    end  
end


A_caps(end-1)=1;
A_sw(end)=1;


if in_cap==1  
A_caps = [zeros(size(A_caps,1),1) A_caps];
A_caps(1)=1;
end

A_sw1 = A_sw(:,1:2:n_stages+2);
A_sw2 = A_sw(:,2:2:n_stages+2);


      
end

