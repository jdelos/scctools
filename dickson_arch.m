function [ ArchDef ] = dickson_arch(n_caps)
%% DICKSON_ARCH Creates the incindence Arch for an N_STAGES (#capacitors)
% Dikson Ladder converter
% 
% Retunr values Arch:
%    ArchDef --> Architecture description file, with the fields: 
%       Acaps   --> Capacitor incidence matrix
%       Asw     --> Switch incidence matrix
%       Asw_atc --> Switch activation matrix
%
%The firs row of the incidence matrix corresponds to the input supply node
%The Phase 1 is the one that drawns power from the load
%
%  Philips Research, Eindhoven,  Netherlands
%  julia.delos@philps.com
%
    
%% Generate the incidence matrixs
[Acaps, A_sw1, A_sw2] = dickson_matrix(n_caps,0);

n_sw_ph1 = size(A_sw1,2);
n_sw_ph2 = size(A_sw2,2);
Asw       = zeros(size(A_sw1,1),n_sw_ph1+n_sw_ph2);
Asw_act = zeros(2,n_sw_ph1+n_sw_ph2);

for i=1:n_sw_ph1
   Asw(:,i*2-1) = A_sw1(:,i); 
   Asw_act(1,i*2-1) = 1;
end

for i=1:n_sw_ph2
   Asw(:,i*2) = A_sw2(:,i); 
   Asw_act(2,i*2) = 1; 
end

ArchDef = struct('Acaps',Acaps,'Asw',Asw,'Asw_act',Asw_act);
    
end





