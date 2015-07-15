function [ ArchDef ] = single_scc_arch
%% SINGLE_SCC_ARCH Creates the incindence Arch for an SINGLE STAGE ARCH
% 
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
Acaps =  [  0  0;
            1  0;
            0  1;];

Asw = [ 1  0;
       -1  1;
        0 -1;];
     
Asw_act = [ 1  0;
            0  1;];
     
ArchDef = struct('Acaps',Acaps,'Asw',Asw,'Asw_act',Asw_act);
    
end





