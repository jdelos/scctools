function [ topology ] = ladder_float_topology(max_ratio,duty,varargin)
%% ladder_float_topology: Create topology of a floating ladder cell
%   Julia Delos, Philips Research
%
%   ladder_float_topology( n_caps,duty_dc)
%       n_stages: conversio ratio
%       duty: duty cylce
%
%   Returns: topology structure with all the relevant functions to model
%   the behaviour of a Hyrbrid Switched Capacitro and compute the losses
%
%   Created 21/03/13, Last modified: ?
%   Copyright 2013-2014, Julia Delos, Philips Research 
%	julia.delos@philips.com	
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.

if nargin == 1 
    duty = 0.5;
end

%%Generate the incidence matrixs
%% Compute number of capacitors
n_caps = max_ratio*2 - 1;

%% Create the incidence matrix 
[A_caps, A_sw1, A_sw2] =  ladder_float_matrix(n_caps);

%Create class 
topo =  generic_switched_capacitor_class(A_caps,A_sw1,A_sw2,'Duty',duty,varargin{:});

%The functions only return the ouput nodes used in hybrid cell with the
%excursion of the all dc-outputs
OutNodes = 1:topo.n_outs;

% if dc_out
%     OutNodes=2:2:topo.n_outs;
% end


%Generate output structures
topology.ratio = topo.m_ratios(OutNodes);
topology.ar = topo.ar(:,OutNodes);
topology.ac = topo.ac(:,OutNodes);
topology.vc = topo.v_caps_norm.'; %Capacitor voltages voltages
topology.vr = topo.v_sw_norm; %Switches voltages
topology.f_ssl = ... %Retrun the symbolic impedance function normalized respect 1Hz 
    topo.r_ssl(OutNodes); 

topology.f_fsl = ... %Retrun the symbolic fsl impedance function of the switches 
    subs(topo.r_fsl(OutNodes),topo.esr_caps,zeros(1,topo.n_caps));

topology.f_esr = ... %Retrun the symbolic fsl impedance function of the esr capacitors 
    subs(topo.r_fsl(OutNodes),topo.ron_switches,zeros(1,topo.n_switches));

topology.f_ssl_in = ... %Retrun the symbolic impedance function normalized respect 1Hz 
    topology.f_ssl./topology.ratio.^2;

topology.f_fsl_in = ... %Retrun the symbolic fsl impedance function of the switches 
    topology.f_fsl./topology.ratio.^2;

topology.f_esr_in = ... %Retrun the symbolic fsl impedance function of the esr capacitors 
    topology.f_esr./topology.ratio.^2;

topology.vo_swing = 1/n_caps;
topology.duty = duty;    
topology.g = topo.k_factors;
topology.g_ssl = topo.k_ssl;
topology.g_fsl = topo.k_fsl;


end

