classdef generic_switched_capacitor_class < handle
    %% SCC Switched Capacitor Class
    %  Creates a SCC_class object that provides all the equations for a given       % SCC structure, of a multi-phase converter
    %
    % The creator takes the arguments:
    %    Ac   -> Capacitor Incidence Matrix
    %	 Asw1 -> Phase 1 switches incidence matrix
    %    Asw2 -> Phase 2 switches incidence matrix
    %    Asw3 -> Phase 3 switches incidence matrix
    %    ....
    %    AswN -> Phase N switch incidence matrix
    %    options -> Vecgtor with the relative duty cycle of each phase
    %            [D1 D2 D3 ... Dn-1 ]
    %  Philips Research, Eindhoven,  Netherlands
    %  julia.delos@philps.com
    %
    
    %% Propoerties
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
        
        SW_active    %Matrix with the active switches
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
        
        function SC = generic_switched_capacitor_class(ArchDesc,varargin)
            
            %% Get Capacitor incidence matrix
            if isfield(ArchDesc,'Acaps')
                SC.inc_caps = ArchDesc.Acaps;
            else
                error('gen_scc:argChk','Missing Acaps field in the ArchDesc structure!')
            end
            
            %% Get Switch incidence matrix
            if isfield(ArchDesc,'Asw')
                SC.inc_switches = ArchDesc.Asw;
            else
                error('gen_scc:argChk','Missing Asw field in the ArchDesc structure!')
            end
            
            %% Get Switch incidence matrix
            if isfield(ArchDesc,'Asw_act')
                SC.SW_active = ArchDesc.Asw_act;
            else
                error('gen_scc:argChk','Missing Asw field in the ArchDesc structure!')
            end
            
            %% Get number of phases
            SC.n_phases=size(SC.SW_active,1);
            
            %% Number of loads
            SC.n_nodes    = size(SC.inc_caps,1);
            SC.n_outs     = SC.n_nodes; %One node is the input
            SC.n_caps     = size(SC.inc_caps,2);
            SC.n_switches = size(SC.inc_switches,2);
            
            %% Parse input arguments
            i=1;
            while i<=(nargin-1)
                if ischar(varargin{i})
                    switch varargin{i}
                        case 'Mode'
                            SC.mode=varargin{i+1};
                            i=i+2;
                        case 'InNode'
                            SC.in_node=varargin{i+1};
                            i=i+2;
                        case 'Duty'
                            try
                                SC.duty = varargin{i+1};
                            catch %#ok<*CTCH>
                                error('SCC:argChk','Wrong Duty vector');
                            end
                            i=i+2;
                            
                        otherwise
                            error('SCC:argChk','Unknown Parameter');
                            
                    end
                end
            end
            
            %% Generate activation matrices
            inc_sw_act = cell(1,SC.n_phases);
            for j = 1:SC.n_phases
                inc_sw_act{j} = SC.inc_switches(:,logical(SC.SW_active(j,:)));
            end
            
            
            %Create converter Incidence Matrix
            SC.incidence_matrix=[SC.inc_caps ...
                inc_sw_act{1:SC.n_phases}];
            
            
            %Switchinf frequency normalitzed respect 1
            SC.fsw_op = 1;
            
            %Initialize symbolic variables
            SC.ron_switches=sym('Ron',[1 SC.n_switches]);
            SC.caps=sym('C',[1 SC.n_caps]);
            SC.esr_caps=sym('Resr',[1 SC.n_caps]);
            
            %Initialize Symbolir ducty vector
            if isempty(SC.duty)
                SC.duty=sym('D',[1 (SC.n_phases-1)]);
            end
            
            %Add last duty cycle
            SC.duty(SC.n_phases)= 1-sum(SC.duty);
            
            %Create Loads incidence matrix
            SC.inc_loads =  eye(SC.n_nodes);
            
            %Initialize the incidence vector of the Voltage supply
            %The positive terminal is 1 and is defined by the InNode option
            SC.supply_branch=zeros(SC.n_nodes,1);
            if max(SC.in_node) <= SC.n_nodes
                if length(SC.in_node) < 2
                    SC.supply_branch(SC.in_node)=1;  %Input supply grounded
                    %In this case one test load is parallel to the voltage supply
                    % it will be removed from the incidence matrix
                    SC.inc_loads(:,SC.in_node) = [];
                else
                    SC.supply_branch(SC.in_node(1:2),1)=[1 -1]';  %Input supply is floating
                end
            end
            
            %Initialize Capacitor incidence matrix
            %SC.input_cap = parallel_elem(SC.supply_branch,SC.inc_caps); %Identify the input capacitor
            
            %if the input source is not connected through an inductor remove
            %capacitor in parallel to the voltage supply
            if ~SC.mode
                %Buck mode
                %Input capacitor is parallel to the supply
                [~,~,ic]=unique([SC.supply_branch SC.inc_caps]','rows','stable');
                %Parallele capacitor branch
                SC.input_cap=find(ic(2:end)==1,1,'last');
                %Eliminate the capacitor branch
                SC.inc_caps(:,SC.input_cap)=[];
                SC.n_caps = size(SC.inc_caps,2);
                
                %Remove the ESR and the Flying caps
                SC.esr_caps(SC.input_cap)=[];
                SC.caps(SC.input_cap)=[];
                SC.n_inputs=1;
                
                SC.n_nodes = size(SC.inc_caps,1);
                SC.n_outs = SC.n_nodes-1; %One node is the input
                
            else
                %Boost mode
                %All ouput load are treated as inputs and the DC input as ouput
                %The capacitors in parallel to the outputs voltage loads must
                %be deleted from the incidence matrix
                %Loads ar considered to be connected to ground
                
                SC.n_caps = size(SC.inc_caps,2);
                
            end
            
            
            
            %initialise phase
            if SC.n_phases >0
                %SC.n_phases= length(varargin); %Number of phases
                ph_lst = 1:SC.n_phases;
                
                
                for i=ph_lst %Switches matrixs
                       
                    %Current mode control function
                    gload =@(x)(SC.duty(i)*(1-x) + x);
                    gsrc  =@(x)(SC.duty(i)*x + (1-x));
                    
                    %Initialize phase
                    SC.phase{i}=SCC_Phase(inc_sw_act{i},...
                        SC.inc_caps,...
                        SC.inc_switches(:,~logical(SC.SW_active(i,:))),...
                        SC.inc_loads*gload(SC.mode),...
                        SC.n_caps,...
                        find(SC.SW_active(i,:)),...
                        SC.supply_branch*gsrc(SC.mode));
                    
                end
                %% #1 Multiphase [JD 5/21/2014]:This is the original code that was created
                % to solve only two phase converters, this code will be commented once the
                % a and b vector solvers are rewritten for mutiphase converters
                if SC.n_phases > 1
                    SC.a_vec_multiphase();
                    SC.b_vec_multiphase();
                    SC.gen_k();
                end
            else
                error('SCC:argChk','No phase matrixs!');
            end
            
            
            %SC.ar = abs( [SC.phase{1}.ar_vector((SC.n_caps+1):end,:);...
            %    SC.phase{2}.ar_vector((SC.n_caps+1):end,:)]);
            
            %Normalized voltage at the caps respect input voltage
            SC.caps_voltage_ratio();
            
            %Normalized voltage at the switches respect input voltage
            SC.switch_voltage_ratio();
            
            %Get dc caps
            dc_caps = find(sum(SC.inc_caps,1));
            SC.dc_out_cap = find(sum([SC.inc_caps(:,dc_caps)...
                SC.inc_loads],2)==2)-1;
        end
        %% [#1,multiphase,JD,5/21/2014]: Function to compute the a vector for a
        % multiphase converter
        
        %% function a_vec_multiphase that solves the charge flow vectors from an SCC.
        % Is the extesion of the function for 2 phase converters
        function a_vec_multiphase(SC)
            Ql = cell(1,SC.n_phases);
            
            %Get cut-set matrix
            for j=1:SC.n_phases
                %Create phase trees
                Tph = tree_ph_scc(SC.phase{j}.get_on_no_sw(),...
                    SC.n_caps,...
                    SC.mode);
                %Generate cut-set matrices
                Ql{j} = fun_cutset(SC.phase{j}.get_on_no_sw(),Tph);
            end
            %% Compute the charge flow vectors
            [ al, m] = solve_charge_vectors(Ql,SC.n_caps,SC.duty);
            
            %% Update the converter structure
            for p=1:SC.n_phases
                SC.phase{p}.set_a_vector(al{p});
            end
            
            %% Update converter conversion ratio
            SC.m_ratios=m;
            SC.m_boost=1./m;
        end
        
        %% Function b_vec_multiphase
        % Replaced the b_vec_2phase providing a solver for the b vectors for
        % multiphase converters
        
        function b_vec_multiphase(SC)
            for i=1:SC.n_phases
                %Short circuit source
                A_sc = short_edge(SC.phase{i}.get_on_no_sw() ,1);
                
                %Generate a tree using only capacitors
                T_sc = build_tree(A_sc(:,1:SC.n_caps));
                n_twings=length(T_sc);
                
                %Fundamental cut-set matrix for the tree
                Qf_sc = fun_cutset(A_sc,T_sc);
                
                %fundamental loop matrix
                Bf_sc = sym(fun_loop(A_sc,T_sc));
                
                %Eliminate loops corresponding to the load sources
                %Number of cap links
                n_cap_links=SC.n_caps-n_twings;
                Bf_sc(n_cap_links+1:end,:)=[];
                
                %Transform capacitor votages to currents by the cap eq.
                % C dV(t)/dt = i(t) => v(S)=i(s)/(sC) s term is droped out
                Bf_sc(:,1:SC.n_caps)=(Bf_sc(:,1:SC.n_caps))/diag(SC.caps);
                
                %Create a system of linear #Caps equations to sovle the current of the
                %capacitors
                %Use all the linear independent equations from Qf and the necessari from the
                %Bf matrix
                S = [ Qf_sc ; Bf_sc(1:SC.n_caps-size(Qf_sc,1),:)];
                
                %Split matrix with variables and excitations [Sx | So] beta solution is
                % b = -Sx\So;
                
                b = -S(:,1:SC.n_caps)\S(:,SC.n_caps+1:end);
                SC.phase{i}.set_b_vector(b);
            end
        end
        
        %Number of switches of the converter
        function n_switch = get_n_switches(SC)
            if SC.n_phases > 0
                n_switch = SC.phase{1}.n_on_sw + SC.phase{1}.n_off_sw;
            end
        end
        
        function gen_k_ssl(SC)
            
            %Compute the square terms of the rows and add them for each
            %phase
            %The resulting matrix is presented in the model format matrix
            r_tot = matrix_2_expon(SC.phase{1}.get_r_vector()) ;
            for i=2:SC.n_phases;
                r_tot = matrix_2_expon(SC.phase{i}.get_r_vector()) + r_tot;
            end
            
            %Weight the factors with the Caps values and add them
            SC.beta= D3_matrix_addition(r_tot,1./SC.caps);
            SC.k_ssl = 1/(2*SC.fsw_op)*SC.beta;
            SC.beta= diag(D3_matrix_addition(r_tot,1./SC.caps));
            
            
        end
        
        function k_ssl_a = gen_k_ssl_a(SC)
            
            %Compute the square terms of the rows and add them for each
            %phase
            %The resulting matrix is presented in the model format matrix
            r_tot = matrix_2_expon(SC.phase{1}.get_ac_vector()) ;
            for i=2:SC.n_phases;
                r_tot = matrix_2_expon(SC.phase{i}.get_ac_vector()) + r_tot;
            end
            
            %Weight the factors with the Caps values and add them
            SC.beta_a= D3_matrix_addition(r_tot,1./SC.caps);
            SC.k_ssl_a = 1/(2*SC.fsw_op)*SC.beta_a;
            SC.beta_a= diag(D3_matrix_addition(r_tot,1./SC.caps));
            k_ssl_a = SC.k_ssl_a;
            
        end
        
        
        function gen_k_fsl(SC)
            %Compute the square terms of the rows and add them for each
            %phase
            %The resulting matrix is presented in the model format matrix
            k_fsl=zeros(SC.n_outs); %#ok<*PROP>
            
                       
            for i=1:SC.n_phases
                if isempty(SC.ron_switches)
                    sw_vector = [SC.esr_caps...
                        sym(['R' num2str(i) '_s'], [1 SC.phase{i}.n_on_sw])];
                else
                    sw_vector = [SC.esr_caps...
                        SC.ron_switches(logical(SC.SW_active(i,:)))];
                end
                
                k_fsl =1./(SC.duty(i)).*...
                    D3_matrix_addition(matrix_2_expon(SC.phase{i}.get_ar_vector()),...
                    sw_vector) + ...
                    k_fsl;
            end
            SC.k_fsl = k_fsl;
            %Weight the factors with the Caps values and add them
        end
        
        function gen_k(SC,varargin)
            %Check k factors
            if isempty(SC.k_fsl)
                SC.gen_k_fsl();
            end
            
            if isempty(SC.k_ssl)
                SC.gen_k_ssl();
            end
            
            %Best aproximation between the limtis is square summ
            %this aproximation should be checked in future work
            SC.k_factors = (SC.k_ssl.^2 + SC.k_fsl.^2).^(1/2);
            SC.r_buck = diag(SC.k_factors);
            SC.k_boost = diag(SC.m_boost.^2)*SC.k_factors;
            SC.r_boost = diag(SC.k_boost);
            SC.r_ssl = diag(SC.k_ssl);
            SC.r_fsl = diag(SC.k_fsl);
            
        end
        
        
        %% Compute normalitzed voltage ration at each capacitor
        %Load effects are not taken into account
        function  caps_voltage_ratio(SC)
            B=[];
            for i=1:SC.n_phases
                B =[B ; inM2loopM(SC.phase{i}.get_on_caps())]; %#ok<AGROW>
            end
            SC.v_caps_norm=(-B(:,2:end)\B(:,1));
        end
        
        %% Compute nomralitzed blocking voltage at the devices
        %Load effecets are not taken into account
        %For two phase converters
        function switch_voltage_ratio(SC)
            %sw_volt=sym(zeros(1,SC.get_n_switches()));
            sw_volt = zeros(SC.n_phases,SC.n_switches);
            for i=1:SC.n_phases
                
                idxs = 1:SC.n_switches;
                for x = SC.phase{i}.sw_idxs
                    idxs(x ==idxs) = [];
                end
                B = inM2loopM(SC.phase{i}.get_on_no_load());
                %Preserve only the loops where the switches are links
                B = B((end-SC.phase{i}.n_off_sw)+1:end,:);
                
                %% Bloking voltage of the other phase switches
                sw_volt(i,idxs) = ...
                    -(B(:,(SC.n_caps+2):end)\B(:,1:(SC.n_caps+1))*...
                    [1; SC.v_caps_norm]); %Solve the blocking voltages
                %sw_volt_max = max(sw_volt,sw_volt_max);
                %sw_volt_min = min(sw_volt,sw_volt_min);
            end
            SC.v_sw_frw = max(sw_volt,[],1);
            SC.v_sw_rev = min(sw_volt,[],1);
            [SC.v_sw_abs, I] = max(abs(sw_volt),[],1);
            SC.v_sw_norm = zeros(1,SC.n_switches);
            for j=1:length(I)
                SC.v_sw_norm(j,I(j))= sw_volt(I(j),j);
            end 
        end
        
    end
end
