classdef generic_switched_capacitor_class < handle
    %% objC Switched Capacitor Class
    %  Creates a objC_class object that provides all the equations for a given       % objC structure, of a multi-phase converter
    %
    % The creator takes the following arguments:
    %    ArchDef --> Architecture description file, with the fields: 
    %       Acaps   --> Capacitor incidence matrix
	%       Asw     --> Switch incidence matrix
	%       Asw_atc --> Switch activation matrix
    %  Philips Research, Eindhoven,  Netherlands
    %  julia.delos@philps.com
	% 
	%  iss2: update --> 14 Nov 2014 (J.Delos)
    %  iss3: enh    --> 24 Feb 2015 (J.Delos)  
    %                   + Adding Rssl, Rfsl, Rser, Rscc functions 
    %                   + Changed class referencing to obj has subggested in
    %                   the Matlab class guidelines.
    
    
    
    %% Properties
    properties (SetAccess = private)
        n_caps;      %Number of capacitors
        dc_out_cap;
        n_outs;      % # of loads/ouputs
        n_nodes;     % # nodes
        n_inputs;    %Number of inputs
        n_switches = 0;  %# switches
        
        
        n_phases=0;    %Number of phases
        
        incidence_matrix; %Complete Converter incidence matrix
        inc_caps;    %Incidence matrix of capacitors
        inc_switches; %Incidence matrix switches
        inc_loads;   %Incidence matrix of loads
        
        %iss4: J.Delos variable name changed to lower case. The name is now
        %more clear 
        sw_activation_matrix    %Matrix with the active switches
        inc_sw_act   %Phase activation matrix
        
        m_ratios;    %Conversion Ratio vector
        m_boost;
        
        k_ssl;       %K ssl factor matrix
        k_ssl_a;     %K ssl factor matrix
        k_fsl;       %K fsl factors
        k_factors;   %square summ of the k_ssl and k_fsl
        r_buck;      %Buck equivalent resitance
        k_boost;     %Boost mode matrix
        r_boost;     %Boost Equivalent resitance
        
        r_ssl;
        r_fsl;
        beta;
        
        beta_a;
        
        
        phase;       %Phase infromation
        duty;        %Duty cycle vector
        
        v_caps_norm;  %Capacitor voltage normailzet respect input voltage
        v_sw_norm;    %Switches blocking voltage normailzet respect input voltage
        v_sw_abs;     %Absoulte blocking voltage of the switches
        v_sw_frw;
        v_sw_rev;
        
        
        caps;         %symbolic parameters
        esr_caps;
        ron_switches; %Switches resitance
        fsw_op;
        
        in_node = 1;
        supply_branch;
        input_cap; %Capacitor parallel to the input
        mode= 0; %0-Buck mode (Current output 0), 1- Boost mode (current input)
        
        ac;
        ar;
    end
    
    methods
        %iss3, jdelos: Changed class referecing name from obj to obj, as
        %suggested in the Matlab guidelines.
        function obj = generic_switched_capacitor_class(ArchDesc,varargin)
          %% iss2: input arguments method toallyu changed %%  
            %% Get Capacitor incidence matrix
            if isfield(ArchDesc,'Acaps')
                obj.inc_caps = ArchDesc.Acaps;
            else
                error('gen_scc:argChk','Missing Acaps field in the ArchDesc structure!')
            end
            
            %% Get Switch incidence matrix
            if isfield(ArchDesc,'Asw')
                obj.inc_switches = ArchDesc.Asw;
            else
                error('gen_scc:argChk','Missing Asw field in the ArchDesc structure!')
            end
            
            %% Get Switch incidence matrix
            if isfield(ArchDesc,'Asw_act')
                obj.sw_activation_matrix = ArchDesc.Asw_act;
            else
                error('gen_scc:argChk','Missing Asw field in the ArchDesc structure!')
            end
            
            %% Check matrices size consitancy
            if (size(obj.inc_switches,1)~=size(obj.inc_caps,1)) 
                error('gen_scc:argChk','Different number of nodes between incidence matrices!')
            end 
            
            if (size(obj.inc_switches,2)~=size(obj.sw_activation_matrix,2))
                 error('gen_scc:argChk','Not consitency between Switch Activation and Incidence matrices!')
            end
            
            %% Get number of phases
            obj.n_phases=size(obj.sw_activation_matrix,1);
            
            %% Number of loads
            obj.n_nodes    = size(obj.inc_caps,1);
            obj.n_outs     = obj.n_nodes; %One node is the input
            obj.n_caps     = size(obj.inc_caps,2);
            obj.n_switches = size(obj.inc_switches,2);
            
            %% Parse input arguments
            i=1;
            while i<=(nargin-1)
                if ischar(varargin{i})
                    switch varargin{i}
                        case 'Mode'
                            obj.mode=varargin{i+1};
                            i=i+2;
                        case 'InNode'
                            obj.in_node=varargin{i+1};
                            i=i+2;
                        case 'Duty'
                            try
                                obj.duty = varargin{i+1};
                            catch %#ok<*CTCH>
                                error('objC:argChk','Wrong Duty vector');
                            end
                            i=i+2;
                            
                        otherwise
                            error('objC:argChk','Unknown Parameter');
                            
                    end
                end
            end
            
            %% Generate activation matrices
            obj.inc_sw_act = cell(1,obj.n_phases);
            for j = 1:obj.n_phases
                obj.inc_sw_act{j} = obj.inc_switches(:,logical(obj.sw_activation_matrix(j,:)));
            end
            
            
            %Create converter Incidence Matrix
            obj.incidence_matrix=[obj.inc_caps ...
                obj.inc_sw_act{1:obj.n_phases}];
            
            %Switchinf frequency normalitzed respect 1
            obj.fsw_op = 1;
            
            %Initialize symbolic variables
            obj.ron_switches=sym('Ron',[1 obj.n_switches]);
            obj.caps=sym('C',[1 obj.n_caps]);
            obj.esr_caps=sym('Resr',[1 obj.n_caps]);
            
            %Initialize Symbolir ducty vector
            if isempty(obj.duty)
                obj.duty=sym('D',[1 (obj.n_phases-1)]);
            end
            
            %Add last duty cycle
            obj.duty(obj.n_phases)= 1-sum(obj.duty);
            
            %Create Loads incidence matrix
            obj.inc_loads =  eye(obj.n_nodes);
            
            %Initialize the incidence vector of the Voltage supply
            %The positive terminal is 1 and is defined by the InNode option
            obj.supply_branch=zeros(obj.n_nodes,1);
            if max(obj.in_node) <= obj.n_nodes
                if length(obj.in_node) < 2
                    obj.supply_branch(obj.in_node)=1;  %Input supply grounded
                    %In this case one test load is parallel to the voltage supply
                    % it will be removed from the incidence matrix
                    obj.inc_loads(:,obj.in_node) = [];
                else
                    obj.supply_branch(obj.in_node(1:2),1)=[1 -1]';  %Input supply is floating
                end
            end
            
            %Initialize Capacitor incidence matrix
            %obj.input_cap = parallel_elem(obj.supply_branch,obj.inc_caps); %Identify the input capacitor
            
            %if the input source is not connected through an inductor remove
            %capacitor in parallel to the voltage supply
            if ~obj.mode
                %Buck mode
                %Input capacitor is parallel to the supply
                [~,~,ic]=unique([obj.supply_branch obj.inc_caps]','rows','stable');
                %Parallele capacitor branch
                obj.input_cap=find(ic(2:end)==1,1,'last');
                %Eliminate the capacitor branch
                obj.inc_caps(:,obj.input_cap)=[];
                obj.n_caps = size(obj.inc_caps,2);
                
                %Remove the ESR and the Flying caps
                obj.esr_caps(obj.input_cap)=[];
                obj.caps(obj.input_cap)=[];
                obj.n_inputs=1;
                
                obj.n_nodes = size(obj.inc_caps,1);
                obj.n_outs = obj.n_nodes-1; %One node is the input
                
            else
                %Boost mode
                %All ouput load are treated as inputs and the DC input as ouput
                %The capacitors in parallel to the outputs voltage loads must
                %be deleted from the incidence matrix
                %Loads ar considered to be connected to ground
                
                obj.n_caps = size(obj.inc_caps,2);
                
            end
            
            
            
            %initialise phase
            if obj.n_phases >0
                %obj.n_phases= length(varargin); %Number of phases
                ph_lst = 1:obj.n_phases;
                
                
                for i=ph_lst %Switches matrixs
                       
                    %Current mode control function
                    gload =@(x)(obj.duty(i)*(1-x) + x);
                    gsrc  =@(x)(obj.duty(i)*x + (1-x));
                    
                    %Initialize phase
                    obj.phase{i}=SCC_Phase(obj.inc_sw_act{i},...
                        obj.inc_caps,...
                        obj.inc_switches(:,~logical(obj.sw_activation_matrix(i,:))),...
                        obj.inc_loads*gload(obj.mode),...
                        obj.n_caps,...
                        find(obj.sw_activation_matrix(i,:)),...
                        obj.supply_branch*gsrc(obj.mode));
                    
                end
                %% #1 Multiphase [JD 5/21/2014]:This is the original code that was created
                % to solve only two phase converters, this code will be commented once the
                % a and b vector solvers are rewritten for mutiphase converters
                if obj.n_phases > 1
                    obj.a_vec_multiphase();
                    obj.b_vec_multiphase();
                    obj.gen_k();
                end
            else
                error('objC:argChk','No phase matrixs!');
            end
            
            
            %obj.ar = abs( [obj.phase{1}.ar_vector((obj.n_caps+1):end,:);...
            %    obj.phase{2}.ar_vector((obj.n_caps+1):end,:)]);
            
            %Normalized voltage at the caps respect input voltage
            obj.caps_voltage_ratio();
            
            %Normalized voltage at the switches respect input voltage
            obj.switch_voltage_ratio();
            
            %Get dc caps
            dc_caps = logical(sum(obj.inc_caps,1));   %iss2: Changed find by logical boost speed
            obj.dc_out_cap = find(sum([obj.inc_caps(:,dc_caps)...
                obj.inc_loads],2)==2)-1;
        end
        
        %% [#1,multiphase,JD,5/21/2014]: Function to compute the a vector for a
        % multiphase converter
        function a_vec_multiphase(obj)
        %% function a_vec_multiphase that solves the charge flow vectors from an objC.
        % Is the extesion of the function for 2 phase converters
            Ql = cell(1,obj.n_phases);
            
            %Get cut-set matrix
            for j=1:obj.n_phases
                %Create phase trees
                Tph = tree_ph_scc(obj.phase{j}.get_on_no_sw(),...
                    obj.n_caps,...
                    obj.mode);
                %Generate cut-set matrices
                Ql{j} = fun_cutset(obj.phase{j}.get_on_no_sw(),Tph);
            end
            %% Compute the charge flow vectors
            [ al, m] = solve_charge_vectors(Ql,obj.n_caps,obj.duty);
            
            %% Update the converter structure
            for p=1:obj.n_phases
                obj.phase{p}.set_a_vector(al{p});
            end
            
            %% Update converter conversion ratio
            obj.m_ratios=m;
            obj.m_boost=1./m;
        end
        
        function b_vec_multiphase(obj)
        %% Function b_vec_multiphase
        % Replaced the b_vec_2phase providing a solver for the b vectors for
        % multiphase converters

            for i=1:obj.n_phases
                %Short circuit source
                A_sc = short_edge(obj.phase{i}.get_on_no_sw() ,1);
                
                %Generate a tree using only capacitors
                T_sc = build_tree(A_sc(:,1:obj.n_caps));
                n_twings=length(T_sc);
                
                %Fundamental cut-set matrix for the tree
                Qf_sc = fun_cutset(A_sc,T_sc);
                
                %fundamental loop matrix
                Bf_sc = sym(fun_loop(A_sc,T_sc));
                
                %Eliminate loops corresponding to the load sources
                %Number of cap links
                n_cap_links=obj.n_caps-n_twings;
                Bf_sc(n_cap_links+1:end,:)=[];
                
                %Transform capacitor votages to currents by the cap eq.
                % C dV(t)/dt = i(t) => v(S)=i(s)/(sC) s term is droped out
                Bf_sc(:,1:obj.n_caps)=(Bf_sc(:,1:obj.n_caps))/diag(obj.caps);
                
                %Create a system of linear #Caps equations to sovle the current of the
                %capacitors
                %Use all the linear independent equations from Qf and the necessari from the
                %Bf matrix
                S = [ Qf_sc ; Bf_sc(1:obj.n_caps-size(Qf_sc,1),:)];
                
                %Split matrix with variables and excitations [Sx | So] beta solution is
                % b = -Sx\So;
                
                b = -S(:,1:obj.n_caps)\S(:,obj.n_caps+1:end);
                obj.phase{i}.set_b_vector(b);
            end
        end      
        
        function n_switch = get_n_switches(obj)
            %Number of switches of the converter
            if obj.n_phases > 0
                n_switch = obj.phase{1}.n_on_sw + obj.phase{1}.n_off_sw;
            end
        end
        
        function gen_k_ssl(obj)
            
            %Compute the square terms of the rows and add them for each
            %phase
            %The resulting matrix is presented in the model format matrix
            r_tot = matrix_2_expon(obj.phase{1}.get_r_vector()) ;
            for i=2:obj.n_phases;
                r_tot = matrix_2_expon(obj.phase{i}.get_r_vector()) + r_tot;
            end
            
            %Weight the factors with the Caps values and add them
            obj.beta= D3_matrix_addition(r_tot,1./obj.caps);
            obj.k_ssl = 1/(2*obj.fsw_op)*obj.beta;
            obj.beta= diag(D3_matrix_addition(r_tot,1./obj.caps));
            
            
        end
        
        function k_ssl_a = gen_k_ssl_a(obj)
            
            %Compute the square terms of the rows and add them for each
            %phase
            %The resulting matrix is presented in the model format matrix
            r_tot = matrix_2_expon(obj.phase{1}.get_ac_vector()) ;
            for i=2:obj.n_phases;
                r_tot = matrix_2_expon(obj.phase{i}.get_ac_vector()) + r_tot;
            end
            
            %Weight the factors with the Caps values and add them
            obj.beta_a= D3_matrix_addition(r_tot,1./obj.caps);
            obj.k_ssl_a = 1/(2*obj.fsw_op)*obj.beta_a;
            obj.beta_a= diag(D3_matrix_addition(r_tot,1./obj.caps));
            k_ssl_a = obj.k_ssl_a;
            
        end
                
        function gen_k_fsl(obj)
            %Compute the square terms of the rows and add them for each
            %phase
            %The resulting matrix is presented in the model format matrix
            k_fsl=zeros(obj.n_outs); %#ok<*PROP>
            
                       
            for i=1:obj.n_phases
                if isempty(obj.ron_switches)
                    sw_vector = [obj.esr_caps...
                        sym(['R' num2str(i) '_s'], [1 obj.phase{i}.n_on_sw])];
                else
                    sw_vector = [obj.esr_caps...
                        obj.ron_switches(logical(obj.sw_activation_matrix(i,:)))];
                end
                
                k_fsl =1./(obj.duty(i)).*...
                    D3_matrix_addition(matrix_2_expon(obj.phase{i}.get_ar_vector()),...
                    sw_vector) + ...
                    k_fsl;
            end
            obj.k_fsl = k_fsl;
            %Weight the factors with the Caps values and add them
        end
        
        function gen_k(obj,varargin)
            %Check k factors
            if isempty(obj.k_fsl)
                obj.gen_k_fsl();
            end
            
            if isempty(obj.k_ssl)
                obj.gen_k_ssl();
            end
            
            %Best aproximation between the limtis is square summ
            %this aproximation should be checked in future work
            obj.k_factors = (obj.k_ssl.^2 + obj.k_fsl.^2).^(1/2);
            obj.r_buck = diag(obj.k_factors);
            obj.k_boost = diag(obj.m_boost.^2)*obj.k_factors;
            obj.r_boost = diag(obj.k_boost);
            obj.r_ssl = diag(obj.k_ssl);
            obj.r_fsl = diag(obj.k_fsl);
            
        end
        
        function  caps_voltage_ratio(obj)
        %% Compute normalitzed voltage ration at each capacitor
        %Load effects are not taken into account
            B=[];
            for i=1:obj.n_phases
                B =[B ; inM2loopM(obj.phase{i}.get_on_caps())]; %#ok<AGROW>
            end
            obj.v_caps_norm=(-B(:,2:end)\B(:,1));
        end
         
        function switch_voltage_ratio(obj)
        %% Compute nomralitzed blocking voltage at the devices
        %Load effecets are not taken into account
        %For two phase converters
        
            %sw_volt=sym(zeros(1,obj.get_n_switches()));
            sw_volt = zeros(obj.n_phases,obj.n_switches);
            for i=1:obj.n_phases
                
                idxs = 1:obj.n_switches;
                for x = obj.phase{i}.sw_idxs
                    idxs(x ==idxs) = [];
                end
                B = inM2loopM(obj.phase{i}.get_on_no_load());
                %Preserve only the loops where the switches are links
                B = B((end-obj.phase{i}.n_off_sw)+1:end,:);
                
                %% Bloking voltage of the other phase switches
                sw_volt(i,idxs) = ...
                    -(B(:,(obj.n_caps+2):end)\B(:,1:(obj.n_caps+1))*...
                    [1; obj.v_caps_norm]); %Solve the blocking voltages
                %sw_volt_max = max(sw_volt,sw_volt_max);
                %sw_volt_min = min(sw_volt,sw_volt_min);
            end
            obj.v_sw_frw = max(sw_volt,[],1);
            obj.v_sw_rev = min(sw_volt,[],1);
            [obj.v_sw_abs, I] = max(abs(sw_volt),[],1);
            obj.v_sw_norm = zeros(1,obj.n_switches);
            for j=1:length(I)
                obj.v_sw_norm(1,j)= sw_volt(I(j),j);
            end 
        end
        
        
        
        %% iss3, jdelos: Adding solver functions to the class
        % Date: 24-Feb-2015        
        function rssl = Rssl(obj,Fsw,Cx,n_outputs)
        %%  function  [rssl] = Rssl(Cx,Fsw,[n_output])
        % Returns the slow switching limit resitance (rssl) value for the topology
        % This function is implemnted for non symbolic duty-cycles
        %
        % Input arguments:
        %   + Cx : Value/s of the capacitors. If scalar all the capacitors are assumed to be equal.
        %   + Fsw: Switching frequency. Scalar or array
        %   +  n_outputs: Array with the necessary outputs to be solved. If not specified all outputs are solved.
        %
        
            %Check Symbolic duty-cycle
            if isa(obj.duty,'sym')
                error('Parameters parser:','Duty cycle is a symbolic vairable')
            end
            
            %Chec for the number of outputs
            if nargin > 3
                s_eqs = obj.r_ssl(n_outputs);
            else
                s_eqs = obj.r_ssl;
            end
            
            %check if Cx is scalar 
            if isscalar(Cx)
                Cx= ones(1,obj.n_caps)*Cx;
            else
                if length(Cx)~= obj.n_caps
                    error('Parameters parser:','Capacitor array length does not matche the number of capacitors in the topology!')
                end
            end
            rssl = double(subs(s_eqs,obj.caps,Cx)./Fsw);
            
        end
        
        function rfsl = Rfsl(obj,Ron,n_outputs)
        %%  function  [rfsl] = Rfsl(Ron,[n_output])
        % Returns the fast switching limit resitance (rfsl) value for the topology
        % This function is implemnted for non symbolic duty-cycles
        %
        % Input arguments:
        %   + Ron       : Value/s of the resistors. If scalar all the capacitors are assumed to be equal.
        %   + n_outputs : Array with the necessary outputs to be solved. If not specified all outputs are solved.
        %
        
            %Check Symbolic duty-cycle
            if isa(obj.duty,'sym')
                error('Parameters parser:','Duty cycle is a symbolic vairable')
            end
                       
            %Chec for the number of outputs
            if nargin > 2
                s_eqs = subs(obj.r_fsl(n_outputs),...
                             obj.esr_caps,...
                             zeros(1,obj.n_caps)) ;
            else
                s_eqs = subs(obj.r_fsl,...
                               obj.esr_caps,...
                               zeros(1,obj.n_caps)) ;
            end
            
            %Remove the Cesr variables
           
            
            %check if Cx is scalar 
            if isscalar(Ron)
                Ron = ones(1,obj.n_switches)*Ron;
            else
                if length(Ron)~= obj.n_switches
                    error('Parameters parser:','Capacitor array length does not matche the number of capacitors in the topology!')
                end
            end
            rfsl = double(subs(s_eqs,obj.ron_switches,Ron));
        end
        
        function resr = Resr(obj,Cesr,n_outputs)
        %%  function  [resr] = Resr(Cesr,[n_output])
        % Returns the fast switching limit resitance (rfsl) value due to
        % due to the Series Equivalent Resitance of the capacitors 
        % This function is implemnted for non symbolic duty-cycles
        %
        % Input arguments:
        %   + Cesr : Value/s of the ESR of the capacitors. If scalar all values are assumed to be equal.
        %   +  n_outputs: Array with the necessary outputs to be solved. If not specified all outputs are solved.
        %
        
            %Check Symbolic duty-cycle
            if isa(obj.duty,'sym')
                error('Parameters parser:','Duty cycle is a symbolic vairable')
            end
            
            %Chec for the number of outputs
            if nargin > 2
                s_eqs = subs(obj.r_fsl(n_outputs),...
                             obj.ron_switches,...
                             zeros(1,obj.n_switches)) ;
            else
                s_eqs = subs(obj.r_fsl,...
                               obj.ron_switches,...
                               zeros(1,obj.n_switches)) ;
            end
            
            %check if Cx is scalar 
            if isscalar(Cesr)
                Cesr = ones(1,obj.n_caps)*Cesr;
            else
                if length(Cesr)~= obj.n_caps
                    error('Parameters parser:','Capacitor array length does not matche the number of capacitors in the topology!')
                end
            end
            resr = double(subs(s_eqs,obj.esr_caps,Cesr));
        end      
        
        function [rscc, fcut] = Rscc(obj,Fsw,Cx,Ron,Cesr,n_outputs)
        %%  function  [resr] = Resr(Cesr,[n_output])
        % Returns the fast switching limit resitance (rfsl) value due to
        % due to the Series Equivalent Resitance of the capacitors 
        % This function is implemnted for non symbolic duty-cycles
        %
        % Input arguments:
        %   + Cx        : Value/s of the capacitors. If scalar all the capacitors are assumed to be equal.
        %   + Fsw       : Switching frequency. Scalar or array.
        %   + Ron       : Value/s of the resistors. If scalar all the capacitors are assumed to be equal.
        %   + Cesr      : Value/s of the ESR of the capacitors. If scalar all values are assumed to be equal.
        %   + n_outputs : Array with the necessary outputs to be solved. If not specified all outputs are solved.
        %
        
            %Check Symbolic duty-cycle
            if isa(obj.duty,'sym')
                error('Parameters parser:','Duty cycle is a symbolic vairable')
            end
            
            %if Number of arguments is bigger than 4 
            % Cesr or n_caps has been specified 
            switch nargin 
                case {4,5} % Cesr has been specified 
                    Rssl = obj.Rssl(Fsw,Cx);
                    Rfsl = obj.Rfsl(Ron);
                    
                    if (nargin == 5) && ~isempty(Cesr)
                        Resr = obj.Resr(Cesr);
                    else
                        Resr = [];
                    end
                    
                case 6
                    Rssl = obj.Rssl(Fsw,Cx,n_outputs);
                    Rfsl = obj.Rfsl(Ron,n_outputs);
                    
                    if ~isempty(Cesr)
                        Resr = obj.Resr(Cesr,n_outputs);
                    else
                        Resr = [];
                    end
            end
            
            rscc = sqrt(Rssl.^2 + (Resr+Rfsl).^2);
            
            if nargout == 2
                % Get the specific Slow Switching Limit normalized w.r.t. 
                % the switching frequency
                Fssl = Rssl(1)*Fsw(1);
                
                %The cut frequency is where the Rssl and (Rfsl+Resr) are 
                % equal
                
                fcut = Fssl/(Rfsl+Resr);
                
                
                
            end
            
        end      
        
        
    end
end
