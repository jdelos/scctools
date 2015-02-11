function [ topology ] = hybrid_topology(Arch,duty,opt)
%% hybrid_topology: Create topology of a floating hybrid cell
%
%   dickson_matrix_hybrid( n_caps,duty,opt)
%       n_stages: number of capacitors
%
%   Returns: topology structure with all the relevant functions to model
%   the behaviour of a Hyrbrid Switched Capacitro and compute the losses
%
%   Created   21/03/13 v0 
%   modified: 20/13/13 v1 Added: Half input step 
%                                Options by struct 
%   Copyright 2013-2014, Julia Delos, Philips Research 
% 	julia.delos@philips.com
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.

dc_out = opt.dc_out;

if (nargin == 1) || isempty(duty) 
    duty = 0.5;
end

%% Create class 
gent_top =  generic_switched_capacitor_class(Arch,'Duty',duty);

%The functions only return the ouput nodes used in hybrid cell with the
%excursion of the all dc-outputs
OutNodes = 1:gent_top.n_outs;

if ~dc_out
    if n_caps > 2
        OutNodes([gent_top.dc_out_cap end])=[];
    else
        OutNodes([gent_top.dc_out_cap])=[];
    end
end


%Generate output structures
topology.ratio = gent_top.m_ratios(OutNodes);
topology.vc = gent_top.v_caps_norm.'; %Capacitor voltages voltages
topology.vr = gent_top.v_sw_norm; %Switches voltages
topology.Y_ssl = gent_top.k_ssl;
topology.Y_fsl = gent_top.k_fsl;

topology.f_ssl = ... %Retrun the symbolic impedance function normalized respect 1Hz 
    gent_top.r_ssl(OutNodes); 

topology.f_fsl = ... %Retrun the symbolic fsl impedance function of the switches 
    subs(gent_top.r_fsl(OutNodes),gent_top.esr_caps,zeros(1,gent_top.n_caps));

topology.f_esr = ... %Retrun the symbolic fsl impedance function of the esr capacitors 
    subs(gent_top.r_fsl(OutNodes),gent_top.ron_switches,zeros(1,gent_top.n_switches));

topology.var_ssl = symvar(topology.f_ssl);

topology.eval_ssl = @(x)... %Returns a function that evaluates the Output Impedance as function
     subs(topology.f_ssl,topology.var_ssl,x); %of flying capacitances 

topology.var_fsl = symvar(topology.f_fsl);
topology.eval_fsl = @(x)... %Returns a function that evaluates the Output Impedance as function
     subs(topology.f_fsl,topology.var_fsl,x); %of flying capacitances 
    
topology.var_fesr = symvar(topology.f_esr);
topology.eval_fesr = @(x)... %Returns a function that evaluates the Output Impedance as function
     subs(topology.f_esr,topology.var_fesr,x); %of flying capacitances 
     
topology.vo_swing = 1/n_caps;
topology.duty = duty;    
topology.g = gent_top.k_factors;

topology.dc_outputs = gent_top.dc_out_cap;
topology.r = cat(3,gent_top.phase{1}.r_vector,gent_top.phase{2}.r_vector); 
topology.r_vars = symvar(topology.r);
topology.eval_r =@(x) subs(topology.r,topology.r_vars,x);
topology.q_dc = [gent_top.phase{1}.r_vector(end,gent_top.dc_out_cap) gent_top.phase{2}.r_vector(end,gent_top.dc_out_cap)];
topology.eval_q_dc = @(x)... %Returns a function that evaluates the Output Impedance as function
     subs(topology.q_dc,symvar(topology.q_dc),x); %of flying capacitances 

topology.N_outs = length(topology.ratio);
topology.N_sw     = gent_top.n_switches;
topology.N_caps   = gent_top.n_caps;

topology.ph = gent_top.phase;
end



function Ao = append_mA(A1,A2)
% Append incidence matrixes 
sze2=size(A2);
sze1=size(A1);
Ao = zeros( sze1+sze2 - [1 0]);
% Append to the converter matrices
Ao(1:sze1(1),1:sze1(2))=A1;
Ao(sze1(1):end,(sze1(2)+1):end)=A2;
end
