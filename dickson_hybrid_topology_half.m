function [ topology ] = dickson_hybrid_topology(n_caps,duty,opt)
%% dickson_matrix_hybrid: Create topology of a floating hybrid cell
%   Julia Delos, Philips Research
%
%   dickson_matrix_hybrid( n_caps )
%       n_stages: number of capacitors
%
%   Returns: topology structure with all the relevant functions to model
%   the behaviour of a Hyrbrid Switched Capacitro and compute the losses
%
%   Created   21/03/13 v0 
%   modified: 20/13/13 v1 Added: Half input step 
%                                Options by struct 
%   Copyright 2013-2014, Julia Delos, Philips Research 
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.

dc_out = opt.dc_out;


if (nargin == 1) || isempty(duty) 
    duty = 0.5;
end

%%Generate the incidence matrixs
[A_caps, A_sw1, A_sw2] = dickson_matrix(n_caps,0);

%Create class 
Dickson =  generic_switched_capacitor_class(A_caps,A_sw1,A_sw2,'Duty',duty);

%The functions only return the ouput nodes used in hybrid cell with the
%excursion of the all dc-outputs
OutNodes = 1:Dickson.n_outs;

if ~dc_out

    if n_caps > 2
        OutNodes([Dickson.dc_out_cap end])=[];
    else
        OutNodes([Dickson.dc_out_cap])=[];
    end
end


%Generate output structures
topology.ratio = Dickson.m_ratios(OutNodes);
topology.ar = Dickson.ar(:,OutNodes);
topology.ac = Dickson.ac(:,OutNodes);
topology.vc = Dickson.v_caps_norm.'; %Capacitor voltages voltages
topology.vr = Dickson.v_sw_norm; %Switches voltages
topology.Y_ssl = Dickson.k_ssl;
topology.Y_fsl = Dickson.k_fsl;

topology.f_ssl = ... %Retrun the symbolic impedance function normalized respect 1Hz 
    Dickson.r_ssl(OutNodes); 

topology.f_fsl = ... %Retrun the symbolic fsl impedance function of the switches 
    subs(Dickson.r_fsl(OutNodes),Dickson.esr_caps,zeros(1,Dickson.n_caps));

topology.f_esr = ... %Retrun the symbolic fsl impedance function of the esr capacitors 
    subs(Dickson.r_fsl(OutNodes),Dickson.ron_switches,zeros(1,Dickson.n_switches));

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
topology.g = Dickson.k_factors;

topology.dc_outputs = Dickson.dc_out_cap;
topology.r = cat(3,Dickson.phase{1}.r_vector,Dickson.phase{2}.r_vector); 
topology.r_vars = symvar(topology.r);
topology.eval_r =@(x) subs(topology.r,topology.r_vars,x);
topology.q_dc = [Dickson.phase{1}.r_vector(end,Dickson.dc_out_cap) Dickson.phase{2}.r_vector(end,Dickson.dc_out_cap)];
topology.eval_q_dc = @(x)... %Returns a function that evaluates the Output Impedance as function
     subs(topology.q_dc,symvar(topology.q_dc),x); %of flying capacitances 
end

