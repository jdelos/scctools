function [ topology ] = dickson_hybrid_topology(n_caps,duty,opt)
%%  topology = dickson_matrix_hybrid(n_caps,duty,opt) 
%Create topology of a floating hybrid cell
%   Input arguments:
%       n_caps --> # of capacitors
%       duty   --> dutcy cycle
%       opt    --> options strucutre
%       opt.dc_out  --> 1 keeps dc outputs in the results 
%       opt.half_point --> 1 adds half point structure in the tolology
%
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

if isfield(opt,'dc_out')
    dc_out = opt.dc_out;
else
   dc_out = 1;  
end

if isfield(opt,'half_point')
    half_point = opt.half_point;
else
    half_point = 0;  
end


if (nargin == 1) || isempty(duty) 
    duty = 0.5;
end

%% Generate the incidence matrixs
%[A_caps, A_sw1, A_sw2] = dickson_matrix(n_caps,0);
 ArchDef  = dickson_arch(n_caps);

%% Add half point conversion
if half_point 
   %Generate half point converter
   A_cap_hp = [0 1 -1 0]';
   A_sw1_hp = [1 0; -1 0; 0 1; 0 -1] ;
   A_sw2_hp = [1 0; 0 1; -1 0; 0 -1];
   
   %Rearrange Dickson matrices
   %Remove Vcc node
   A_caps(1,:) = []; 
   A_sw1(1,:)  = [];  
   A_sw2(1,:)  = []; 
   
   %Remove top swithc is the first swhitch of phase1

   A_sw1(:,1)  = [];   
   A_sw1(:,1)  = [];  
   

   A_caps = append_mA(A_cap_hp,A_caps);
   A_sw1  = append_mA(A_sw1_hp,A_sw1);
   A_sw2  = append_mA(A_sw2_hp,A_sw2);
end

%% Create class 
top =  generic_switched_capacitor_class(ArchDef,'Duty',duty);

%The functions only return the ouput nodes used in hybrid cell with the
%excursion of the all dc-outputs
OutNodes = 1:top.n_outs;

if ~dc_out
    if n_caps > 2
        OutNodes([top.dc_out_cap end])=[];
    else
        OutNodes([top.dc_out_cap])=[];1
    end
end


%Generate output structures
topology.ratio = top.m_ratios(OutNodes);
%topology.ar = Dickson.ar(:,OutNodes);
%topology.ac = Dickson.ac(:,OutNodes);
topology.vc = top.v_caps_norm.'; %Capacitor voltages voltages
topology.vr = top.v_sw_norm; %Switches voltages
topology.Y_ssl = top.k_ssl;
topology.Y_fsl = top.k_fsl;

topology.f_ssl = ... %Retrun the symbolic impedance function normalized respect 1Hz 
    top.r_ssl(OutNodes); 

topology.f_fsl = ... %Retrun the symbolic fsl impedance function of the switches 
    subs(top.r_fsl(OutNodes),top.esr_caps,zeros(1,top.n_caps));

topology.f_esr = ... %Retrun the symbolic fsl impedance function of the esr capacitors 
    subs(top.r_fsl(OutNodes),top.ron_switches,zeros(1,top.n_switches));

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
topology.g = top.k_factors;

topology.dc_outputs = top.dc_out_cap;
topology.r = cat(3,top.phase{1}.r_vector,top.phase{2}.r_vector); 
topology.r_vars = symvar(topology.r);
topology.eval_r =@(x) subs(topology.r,topology.r_vars,x);
topology.q_dc = [top.phase{1}.r_vector(end,top.dc_out_cap) top.phase{2}.r_vector(end,top.dc_out_cap)];
topology.eval_q_dc = @(x)... %Returns a function that evaluates the Output Impedance as function
     subs(topology.q_dc,symvar(topology.q_dc),x); %of flying capacitances 

topology.N_outs   = length(topology.ratio);
topology.N_sw     = top.n_switches;
topology.N_caps   = top.n_caps;

topology.ph       = top.phase;
topology.Rssl     = @(Fsw,Cx,n_outputs) top.Rssl(Fsw,Cx,n_outputs);
topology.Rfsl     = @(Ron,n_outputs)    top.Rfsl(Ron,n_outputs);
topology.Resr     = @(Cesr,n_outputs)   top.Resr(Cesr,n_outputs);
topology.Rscc     = @(Fsw,Cx,Ron,Cesr,n_outputs)    top.Rscc(Fsw,Cx,Ron,Cesr,n_outputs);
% #4: Necessary to pass the entire genertic topology class for the
% parasitic functionality 
topology.g_top         = top;

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
