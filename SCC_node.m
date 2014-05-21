classdef SCC_node < handle
    %SCC Switched Capacitor Class
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
       n_caps;      %Number of capacitors
       dc_out_cap;     
       n_outs;      % # of loads/ouputs
       n_nodes;     % # nodes
       n_inputs;    %Number of inputs
       n_switches = 0;  %# switches
       ron_switches; %Switches resitance 
       ds_caps;     %Parasitic Capacitor Drain source 
       caps;        %Capacitor Vector
       esr_caps;    %Capacitor ESR
       n_phases=0;    %Number of phases
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
       fsw_op;      %Switching frequency
       v_caps_norm; %Capacitor voltage normailzet respect input voltage
       v_sw_norm;   %Switches blocking voltage normailzet respect input voltage
       voltage_switches %Votlage of the switches
       current_switches %Current through the switch
       current_caps_pk %Current through the switch
                    %Converter inputs and ouputs
       vin=sym('Vin');
       current_vector;
       load_resistances;
       voltage_switches_on; %Voltage drop on swuitches
       
       incidence_matrix; %Complete Converter incidence matrix 
       in_node = 1;
       supply_branch;
       input_cap; %Capacitor parallel to the input
       mode= 0; %0-Buck mode (Current output 0), 1- Boost mode (current input)
       
       ac;
       ar;
    end
    
    methods
        function SC = SCC_node(inc_caps,varargin)
        
        
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
               case 'Fsw'
                   SC.fsw_op=varargin{i+1};
                   i=i+2;
               case 'Cfly'
                   try
                     if length(varargin{i+1})==1
                        SC.caps = ones(1,SC.n_caps)*varargin{i+1};
                     else
                       SC.caps = varargin{i+1};
                     end
                   
                   catch %#ok<*CTCH>
                   error('SCC:argChk','Wrong Capacitor vector');
                   end
                   i=i+2;
                   
                    
               case 'Duty'
                   try
                       SC.duty = varargin{i+1};                  
                   catch %#ok<*CTCH>
                   error('SCC:argChk','Wrong Duty vector');
                   end
                   i=i+2;
                   
               case 'Ron'
                   try
                     if length(varargin{i+1})==1
                        SC.ron_switches = ones(1,SC.n_switches)*varargin{i+1};
                     else
                       SC.ron_switches = varargin{i+1};  
                     end
                   catch
                   error('SCC:argChk','Wrong Switch ON resistance vector');
                   end
                   i=i+2;
                   
                   
               case 'Cds'
                   try
                     if length(varargin{i+1})==1
                        SC.ds_caps = ones(1,SC.n_switches)*varargin{i+1};
                     else
                       SC.ds_caps = varargin{i+1};
                     end
                   catch
                   error('SCC:argChk','Wrong Drain Source Capacitor vector');
                   end
                   i=i+2;
                   
               case 'ESR'  
                   try
                     if length(varargin{i+1})==1
                        SC.esr_caps = ones(1,SC.n_caps)*varargin{i+1};
                     else
                       SC.esr_caps = varargin{i+1};
                     end
                   
                   catch
                   error('SCC:argChk','Wrong ESR Capacitor vector');
                   end
                   i=i+2;
                   
               case 'Vin'
                   try
                     if length(varargin{i+1})==1
                         
                        SC.vin = varargin{i+1};
                     else
                        error('SCC:argChk','Wrong Vin is a vector. Itmust be scalar value!');
                    end
                   
                   catch
                   error('SCC:argChk','Wrong Vin Parameter');
                   end
                   i=i+2;
               
               case 'Io'
                   if isempty(SC.load_resistances)
                   try
                     if length(varargin{i+1})==1  
                        SC.current_vector = varargin{i+1};
                     else
                        SC.current_vector = varargin{i+1};
                     end
                   
                   catch
                   error('SCC:argChk','Wrong Output Current Vector');
                   end
                   i=i+2;
                   else
                    disp('Load Resistances has been already defined. They will be used for the results');
                   end
                   
               case 'Ro'
                   if ~isempty(SC.current_vector)
                       disp('Load Currents has been already defined. They will ingored for the results');
                       SC.current_vector = [];
                   end
                   
                   try
                     if length(varargin{i+1})==1 
                        SC.load_resistances(1,:) = varargin{i+1};
                     else
                        SC.load_resistances= varargin{i+1};
                     end
                   
                   catch
                   error('SCC:argChk','Wrong Load Resistances Vector');
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
       
       
       if isempty(SC.fsw_op)
        %Switching frequency
        SC.fsw_op = sym('Fsw');
       end  
       
       if isempty(SC.ron_switches)    
        %Symbolic caps
        SC.ron_switches=sym('Ron',[1 SC.n_switches]);
       end
       
       if isempty(SC.caps)
        %Symbolic caps
        SC.caps=sym('C',[1 SC.n_caps]);
       end
       
       if isempty(SC.esr_caps)
          %Symbolic ESR
          SC.esr_caps=sym('Resr',[1 SC.n_caps]);
       end
       
       if isempty(SC.ds_caps)
          %Symbolic ESR
          SC.ds_caps=sym('Cds',[1 SC.n_switches]);
       end

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
          
        %Create symbolic load resistance vector
       if isempty(SC.load_resistances)  
        SC.load_resistances=sym('Rout', [1 SC.n_outs]);
       elseif length(SC.load_resistances)~=SC.n_outs
        SC.load_resistances(1:SC.n_outs)=SC.load_resistances(1);
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
                
            SC.caps_voltage_ratio();
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
        
        %Compute normalitzed voltage ration at each capacitor
        %Load effects are not taken into account
        function caps_voltage_ratio(SC)
            B=[];
            for i=1:SC.n_phases
             B =[B ; inM2loopM(SC.phase{i}.get_on_caps())]; %#ok<AGROW>
            end
            SC.v_caps_norm=B(:,2:end)\B(:,1);        
        end       
        
        %Compute nomralitzed blocking voltage at the devices
        %Load effecets are not taken into account
        function switch_voltage_ratio(SC)
            %sw_volt=sym(zeros(1,SC.get_n_switches()));
            indx=1;
            for i=1:SC.n_phases
             B = inM2loopM(SC.phase{i}.get_on_no_load());
             %Preserve only the loops where the switches are links
             B = B((end-SC.phase{i}.n_off_sw)+1:end,:);
             
             sw_volt(indx:(SC.phase{i}.n_off_sw+(indx-1)) ) = ...
              -(B(:,(SC.n_caps+2):end)\B(:,1:(SC.n_caps+1))*[1; SC.vcaps()]).';
             indx=indx+SC.phase{i}.n_off_sw;
            end
            SC.v_sw_norm = sw_volt;       
        end
        
        %Number of switches of the converter
        function n_switch = get_n_switches(SC)
            if SC.n_phases > 0
                n_switch = SC.phase{1}.n_on_sw + SC.phase{1}.n_off_sw;
            end
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
        
        function x = output_current(SC,varargin)
            if isempty(SC.current_vector)
             
             %Initialize current vector
             x = sym(zeros(SC.n_outs,1));  
                
             
             if isempty(symvar(SC.load_resistances))
             %Get non infinite values
             idx_val =  ~isinf(SC.load_resistances);
             
                 if any(idx_val)
                     %Create reduced matrixs
                     m_red  = SC.m_ratios(idx_val,1);
                     Rx_red = diag(SC.load_resistances(idx_val));
                     k_red  = SC.k_factors(idx_val,idx_val);

                     %Compute current values
                     x(idx_val,1) = (Rx_red+k_red)...
                                     \(m_red*SC.vin);
        
                 end
             
             else
              Rx = diag(SC.load_resistances);  
              k = SC.k_factors;
              m  = SC.m_ratios;
              x = (Rx+k)...
                   \(m*SC.vin);
             end
                         
             %Write current to interal variable            
             SC.current_vector = x;        
            else
                x = SC.current_vector;
            end
            
            %Symbolic functions replacement
            if nargin > 1
                x = subs_func(x,varargin{:});
           end
        end
                
        function x = drop_loss(SC,varargin)
            if isempty(SC.current_vector)
                [~] = SC.output_current;
             end
            x = diag(SC.output_current)*SC.k_factors*SC.output_current;
            
            %Symbolic functions replacement
            if nargin > 1
             x =  subs_func(x,varargin{:});
           end
        end
        
        function x = sw_loss(SC,varargin)
           x = 1/2*SC.fsw_op*((SC.vswitches_blk()).^2*SC.ds_caps.');  
           
         %Symbolic functions replacement
            if nargin > 1
               x = subs_func(x,varargin{:});
           end
            
        end
        
        function x = vout(SC,varargin)
             if isempty(SC.current_vector)
                [~] = SC.output_current;
             end
            x = SC.m_ratios*SC.vin - SC.k_factors*SC.current_vector;    
            
            %Symbolic functions replacement
            if nargin > 1
                x = subs_func(x,varargin{:});
           end
            
        end
        
        function x = vcaps(SC,varargin)
            
            
            x = ([SC.vin SC.vout().']*SC.inc_caps).';    
            
            %Symbolic functions replacement
            if nargin > 1
                x = subs_func(x,varargin{:});
           end
            
        end
        
        function x = vswitches(SC,varargin)
            if isempty(SC.voltage_switches)
                x=sym(zeros(1,SC.get_n_switches()));

                %Switches array indexs
                sw_idx=1:SC.get_n_switches();

                for i=1:SC.n_phases

                 B = inM2loopM(SC.phase{i}.get_on_no_load());
                 %Preserve only the loops where the switches are links
                 %First loop row never includes the voltages of the
                 %switcehs 
                 %The last rows allways includes them, select the last row
                 %switches
                 B = B((end-SC.phase{i}.n_off_sw)+1:end,:);

                 %The voltage convetion for the supply is reversed from the
                 %deffinition of the Capacitor Incidence Matrix there the
                 %conventio was made following the charge transfer ratio

                 %Change sign first column
                 B(:,1)=-B(:,1);

                 %Current B matrix defines [B]*[Vx] = 0 
                 %Where Vx is a column vector as:
                 % Vx = [Vin Vc1 Vc2 ...Vcn Vsw1 ..Vswn]'
                 %
                 %The matrix can be splited as 
                 % [Bc | Bsw ] [ Vc | Vsw]'= 0;
                 %which yelds to
                 %
                 % Bc*Vc +Bsw*Vsw=0
                 %
                 %so the solution of the switches voltages is
                 % 
                 % Vsw = -Bsw\Bc*Vc
                 %
                 %where
                 %
                 % Vc = [Vin Vc1 Vc2 ...Vcn]'
                 % Vsw = [Vsw1 Vsw2 ... Vswn]'


                 %Define the Capacitor Voltage Loop Matrix
                 Bc=B(:,1:(SC.n_caps+1));


                 %Define the Switch Voltage loop Matrix
                 Bsw=B(:,(SC.n_caps+2):end);

                 %Get the indcies of the off switches
                 idx_off=~ismember(sw_idx,SC.phase{i}.sw_idxs);

                 %Solve the voltages for the off switches
                 x(idx_off) = -Bsw\Bc*[SC.vin; SC.vcaps()];

                end
                SC.voltage_switches=x;
            else
                x=SC.voltage_switches;
            end

            %Symbolic functions replacement
            if nargin > 1
                x = subs_func(x,varargin{:});
           end
            
        end
        
        function x = vswitches_blk(SC,varargin)
          %Returns the blocking voltage of the switches
          %For a two-Phase converter
          
            if isempty(SC.voltage_switches)
                %x=sym(zeros(1,SC.get_n_switches()));
                
                
                A_supply = [ 1; zeros(size(SC.incidence_matrix,1)-1,1)];
                %Create fundemental loop matrix of the switchs
                [B ,~] = fun_loop([A_supply SC.incidence_matrix]);
                
                %Create the sub matrix Bx Bsw_I Bsw_I
                % [Bx | Bsw_I | Bsw_II] * [Vx | Vsw_I | Vsw_II]' = 0
                %
                
                Bx = B(:,1:(SC.n_caps+1));
                Bsw_I = B(:,(SC.n_caps+1) +( 1:SC.phase{1}.n_on_sw));
                Bsw_II =B(:,(SC.phase{1}.n_on_sw+SC.n_caps+1)...
                    +( 1:SC.phase{2}.n_on_sw));
                
                %Get switches on switches voltages
                Von_sw = SC.ron_switches .* SC.iswitches();
                
                %Devide it by phases
                Von_sw_I = Von_sw(SC.phase{1}.sw_idxs).';
                Von_sw_II = Von_sw(SC.phase{2}.sw_idxs).';
                
                %Create Vx voltages array
                Vx = [SC.vin; SC.vcaps()];
                                
                %Solve for phase I
                indep_idx = indep_rows(Bsw_I);
                Vbk_I = -Bsw_I(indep_idx,:)\...
                    [Bx(indep_idx,:) Bsw_II(indep_idx,:)]*[Vx; Von_sw_II];
                
                %Solve for phase II
                indep_idx = indep_rows(Bsw_II);
                Vbk_II = -Bsw_II(indep_idx,:)\...
                    [Bx(indep_idx,:) Bsw_I(indep_idx,:)]*[Vx; Von_sw_I];
                
                SC.voltage_switches=[Vbk_I.' Vbk_II.'];
                
            end
                x=SC.voltage_switches;
            
            %Symbolic functions replacement
            if nargin > 1
                x = subs_func(x,varargin{:});
           end
            
        end
        
        function x = vswitches_on(SC,varargin)
          %Returns the blocking voltage of the switches
          %For a two-Phase converter
          
            if isempty(SC.voltage_switches_on)
               %Get switches on switches voltages
                SC.voltage_switches_on = SC.ron_switches .* SC.iswitches();
            end
                x=SC.voltage_switches_on;
            
            %Symbolic functions replacement
            if nargin > 1
                x = subs_func(x,varargin{:});
           end
            
        end
             
        function x = iswitches(SC,varargin)
            if isempty(SC.current_switches)
            %Init. return vector
            x = sym(zeros(1,SC.n_switches));
            
            %Compute the phases
            for i=1:SC.n_phases 
                x(SC.phase{i}.sw_idxs)=...
                    (SC.phase{i}.ar_vector((SC.n_caps+1):end,:)*...
                    SC.output_current())./(SC.duty(i));
            end
            %Store to the class
            SC.current_switches = x;
           
            else
                x= SC.current_switches;
            end
            
            %Symbolic functions replacement
            if nargin > 1
                x = subs_func(x,varargin{:});
            end
        end
        
        function x = icaps_pk(SC,varargin)
            if isempty(SC.current_caps_pk)
            
            %Init. return vector
            x = sym(zeros(1,SC.n_caps*SC.n_phases));
            
            %Compute the phases
            for i=1:SC.n_phases 
                x((i-1)*SC.n_caps+(1:SC.n_caps))=...
                    (SC.phase{i}.ar_vector(1:SC.n_caps,:)*...
                    SC.output_current())./(SC.duty(i));
            end
            SC.current_caps_pk = x;
            
            else
                x= SC.current_caps_pk;
            end
            
             %Symbolic functions replacement
            if nargin > 1
                x = subs_func(x,varargin{:});
            end
        end
        
        
    end
    
end

function x = subs_func(x,varargin)

if ~isempty(symvar(x))
    if strcmp(varargin{1},'Func')                
        x =symfun(x,symvar(x));
    else  
        %Symbolic values must be replaced
        x = subs(x,varargin{:});             
    end
end
   
end

function x = indep_rows(B)
    
    Baux=B;
    [~, n_cols] = size(B);
    x = zeros(1,n_cols);
    
    %Check unique soultions
%     for i=1:n_cols
%         single_cols = sum(abs(Baux),1)==1;
%         if any(single_cols)
%             [x(single_cols), ~] = find(Baux(:,single_cols)~=0);
%             %Eliminate the selected row makeing them 0
%             Baux(x(x~=0),:)=zeros(1,n_cols);
%         else
%             break
%         end
%     end
    
    
    while any(x==0)
        single_cols = sum(abs(Baux),1)==1;
        if any(single_cols)
            [x(single_cols), ~] = find(Baux(:,single_cols)~=0);
            %Eliminate the selected row makeing them 0
            Baux(x(x~=0),:)=0;
        else
            i=find(x==0,1,'first');
            idx_rows = find((Baux(:,i)~=0)');
            %Check if not the first iteration
            if any(x) 
                %Check that the row is linear dependent of the previous
                %sleceted
                for idx = idx_rows
                %Check if the current row is independent    
                    if ~dependent_row(B(x(x~=0),:),...
                            Baux(idx,:))
                        x(i) = idx;
                        Baux(idx,:)=0;
                        Baux(:,i)=0;
                        break
                    end
                end
            end
        end
    end
    
    if length(unique(x)) ~=  n_cols
        disp('indep_rows:MtrixCorrupt','Input matrix has a dependent rows!')
    end
end

function y = dependent_row(A,r_vec)

y=false;

if any(r_vec) && (size(A,2) == size(r_vec,2)) %Cehck zero row
    %r_vec is not a zero row
    for i=1:size(A,1)
     t_row = A(i,:);
     if any(t_row) %Test row is no 0
         test_multiples = t_row./r_vec;
         %Remove NaN zero by zero division
         test_multiples(isnan(test_multiples)) = [];
         y = (length(unique(test_multiples)) == 1);
         if y 
             break
         end
     end
        
    end
else
    y = true;
end

end

function edges = parallel_elem(edge,A) %#ok<DEFNU>

n_columns =size(A,2);
%Alocate vector
edges = zeros(n_columns,1);
Ax=abs(A);

%Make target colum different
for i =1:n_columns
    if i~= edge && all(Ax(:,i)==Ax(:,edge))
     edges(i)=i;
    end
end
%Clean zeor entries
edges(edges==0)=[];     
end