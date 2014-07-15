classdef generic_switched_capacitor_class < handle
    %SCC Switched Capacitor Class
    %  Creates a SCC_class object that provides all the equations for a given 
	% SCC structure, of a two phase converter
	%
	% The creator takes 3 arguments:\
	%    Ac   -> Capacitor Incidence Matrix
	%	 Asw1 -> Phase 1 switches incidence matrix
	%    Asw2 -> Phase 2 switches incidence matrix 
	%
	% 
	%  Philips Research, Eindhoven,  Netherlands
	%  julia.delos@philps.com
	%
    
	
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
       inc_loads;   %Incidence matrix of loads
       
       m_ratios;    %Conversion Ratio vector
       m_boost;
       
       k_ssl;       %K ssl factor matrix
       k_fsl;       %K fsl factors
       k_factors;   %square summ of the k_ssl and k_fsl  
       r_buck;      %Buck equivalent resitance
       k_boost;     %Boost mode matrix
       r_boost;     %Boost Equivalent resitance  
       
       r_ssl;
       r_fsl;
       beta;
       
       phase;       %Phase infromation
       duty;        %Duty cycle vector
       
       v_caps_norm; %Capacitor voltage normailzet respect input voltage
       v_sw_norm;   %Switches blocking voltage normailzet respect input voltage
                    
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
        function SC = generic_switched_capacitor_class(inc_caps,varargin)
        
        
        %Number of loads
        SC.n_nodes = size(inc_caps,1);
        SC.n_outs = SC.n_nodes; %One node is the input  
        SC.n_caps = size(inc_caps,2);
                
           
       %Parse input arguments
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
           
           else
               SC.n_phases=1+SC.n_phases;
               
               SC.n_switches = SC.n_switches + ...
                   size(varargin{i},2);
               i=i+1;
           end
       end
       
       %Create converter Incidence Matrix
       SC.incidence_matrix=[inc_caps ...
           varargin{1:SC.n_phases}];
       
       
       
       SC.fsw_op = 1;
       
       %Initialize symbolic variables
       SC.ron_switches=sym('Ron',[1 SC.n_switches]);     
       SC.caps=sym('C',[1 SC.n_caps]);
       SC.esr_caps=sym('Resr',[1 SC.n_caps]);
       
       %Initialize Symbolir ducty vector
       if isempty(SC.duty)         
        SC.duty=sym('D',1:(SC.n_phases-1));              
       end
       
       %Add last duty cycle
       SC.duty(SC.n_phases)= 1-sum(SC.duty);
       
        %Create Loads incidence matrix
        SC.inc_loads =  eye(SC.n_nodes);  
       
           %Init Supply branch is a column vector
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
        %SC.input_cap = parallel_elem(SC.supply_branch,inc_caps); %Identify the input capacitor 
        
        %if the input source is not connected through an inductor remove
        %capacitor in parallel to the voltage supply
        if ~SC.mode
            %Buck mode
            %Input capacitor is parallel to the supply    
            [~,~,ic]=unique([SC.supply_branch inc_caps]','rows','stable');
            %Parallele capacitor branch
            SC.input_cap=find(ic(2:end)==1,1,'last');
            %Eliminate the capacitor branch
            inc_caps(:,SC.input_cap)=[]; 
            SC.inc_caps = inc_caps; 
            SC.n_caps = size(inc_caps,2);
            
            %Remove the ESR and the Flying caps
            SC.esr_caps(SC.input_cap)=[];
            SC.caps(SC.input_cap)=[];
            SC.n_inputs=1;
            
            SC.n_nodes = size(inc_caps,1);
            SC.n_outs = SC.n_nodes-1; %One node is the input  
            
        else
            %Boost mode
            %All ouput load are treated as inputs and the DC input as ouput
            %The capacitors in parallel to the outputs voltage loads must
            %be deleted from the incidence matrix
            %Loads ar considered to be connected to ground
            
            SC.inc_caps = inc_caps; 
            SC.n_caps = size(inc_caps,2);
             
        end
          
           
        
        %initialise phase 
            if SC.n_phases >0  
              %SC.n_phases= length(varargin); %Number of phases
              ph_lst = 1:SC.n_phases;
              
              %Init Switch index
              sw_idx=1;
                for i=ph_lst %Switches matrixs
                  idx=ph_lst;
                  idx(i)=[];
                  
                  
                    %Current mode control function
                    gload =@(x)(SC.duty(i)*(1-x) + x); 
                    gsrc  =@(x)(SC.duty(i)*x + (1-x));
                    
                    %Initialize phase 
                    SC.phase{i}=SCC_Phase(varargin{i},...
                    inc_caps,varargin{idx},...
                    SC.inc_loads*gload(SC.mode),...
                    SC.n_caps,...
                    sw_idx,...
                    SC.supply_branch*gsrc(SC.mode));
                    
                    sw_idx= sw_idx + size(varargin{i},2);
                end
                     if SC.n_phases ==2
                        SC.a_vec_2phase(); %Compute a vectors 
                                           %It computes the conversion ratio 
                        SC.b_vec_2phase(); %Compute b vectors
                        SC.gen_k(); %Generate k_matrix

                     end
                     %Capacitor voltage ratios
                     %SC.caps_voltage_ratio();

                     %Switches voltage ratios
                     %SC.switch_voltage_ratio();
            else
                 error('SCC:argChk','No phase matrixs!'); 
            end
            
            
            SC.ar = abs( [SC.phase{1}.ar_vector((SC.n_caps+1):end,:);...
                    SC.phase{2}.ar_vector((SC.n_caps+1):end,:)]);
            
            %Normalized voltage at the caps respect input voltage
            SC.caps_voltage_ratio();
            
            %Normalized voltage at the switches respect input voltage
            SC.switch_voltage_ratio();
           
            %Get dc caps
            dc_caps = find(sum(SC.inc_caps,1));
            SC.dc_out_cap = find(sum([SC.inc_caps(:,dc_caps)...
                            SC.inc_loads],2)==2)-1;
        end
       
        function a_vec_2phase(SC)
          %Create phase trees
          [ Tph1, Tph2 ] = tree_2ph_scc(SC.phase{1}.get_on_no_sw(),...
                 SC.phase{2}.get_on_no_sw(),...
                 SC.n_caps,...
                 SC.mode);   
          %Cut-set matrixs
          
          
          Ql_p1 = fun_cutset(SC.phase{1}.get_on_no_sw(),Tph1);
          Ql_p2 = fun_cutset(SC.phase{2}.get_on_no_sw(),Tph2);
             
          %Solve a vector Ph1 
          [a1, a2, m] = solve_2ph_charge_vectors(Ql_p1,Ql_p2,...
                                                 SC.n_caps,...
                                                 SC.mode,...
                                                 SC.duty);
          %Save a vectors to the phases
          SC.phase{1}.set_a_vector(a1);
          SC.phase{2}.set_a_vector(a2);
          
          %set the ac vectors
          SC.ac=abs(a1(2:end,:));
          
          
          
          %Get conversion ratio from ain1+ain2
          SC.m_ratios=m;  
          SC.m_boost=1./m; 
          
          
         
        end
        
        function b_vec_2phase(SC)
         for i=1:2
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

            b=-S(:,1:SC.n_caps)\S(:,SC.n_caps+1:end);
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
       
        
        function gen_k_fsl(SC)
            %Compute the square terms of the rows and add them for each
            %phase
            %The resulting matrix is presented in the model format matrix
            k_fsl=zeros(SC.n_outs); %#ok<*PROP>
            
            %Offset index switches 
            sw_offset=1;
            
            for i=1:SC.n_phases
               if isempty(SC.ron_switches)
                sw_vector = [SC.esr_caps...
                    sym(['R' num2str(i) '_s'], [1 SC.phase{i}.n_on_sw])];
               else
                sw_vector = [SC.esr_caps...
                    SC.ron_switches(sw_offset:SC.phase{i}.n_on_sw+(sw_offset-1))];
                sw_offset = SC.phase{i}.n_on_sw+1;
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
        
         
        %Compute normalitzed voltage ration at each capacitor
        %Load effects are not taken into account
        function  caps_voltage_ratio(SC)
            B=[];
            for i=1:SC.n_phases
             B =[B ; inM2loopM(SC.phase{i}.get_on_caps())]; %#ok<AGROW>
            end
            SC.v_caps_norm=(-B(:,2:end)\B(:,1));        
        end       
        
        %Compute nomralitzed blocking voltage at the devices
        %Load effecets are not taken into account
        %For two phase converters
        function switch_voltage_ratio(SC)
            %sw_volt=sym(zeros(1,SC.get_n_switches()));
            idxs = {SC.phase{2}.sw_idxs, SC.phase{1}.sw_idxs };
            sw_volt = zeros(1,SC.n_switches);
            for i=1:SC.n_phases
             B = inM2loopM(SC.phase{i}.get_on_no_load());
             %Preserve only the loops where the switches are links
             B = B((end-SC.phase{i}.n_off_sw)+1:end,:);
             
             
             %% Bloking voltage of the other phase switches
             sw_volt( idxs{i} ) = ...
              -(B(:,(SC.n_caps+2):end)\B(:,1:(SC.n_caps+1))*[1; SC.v_caps_norm]);
             
            end
            SC.v_sw_norm = sw_volt;       
        end
        
    end
end