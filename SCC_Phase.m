classdef SCC_Phase < handle
    
    %%SCC Switched Capacitor Phase Class
    % Creates the necessary matrices to solve the converter composed by
	% different phases 
	%
	%
	% 
	%
	%Copyright 2013-2014, Julia Delos, Philips Research 
	%	julia.delos@philips.com	
	%   May be freely used and modified but never sold.  The original author
	%   must be cited in all derivative work
	
    properties (SetAccess = private)
      n_caps;
      n_loads;
      n_on_sw;
      n_off_sw;
      inc_on_sw;
      inc_on_conv;
      inc_on_conv_sw;
      a_vector;
      a_vector_up;
      ar_vector;
      b_vector;
      r_vector;
      sw_idxs;
    end
    
    methods
        %Class creator
        function PH = SCC_Phase(inc_sw_ph,inc_caps,inc_sw,inc_loads,n_caps,idx_sw,supply_branch)
            PH.inc_on_sw = inc_sw_ph;
            PH.n_on_sw=size(inc_sw_ph,2);
            PH.inc_on_conv = phase_conv([supply_branch inc_caps inc_loads inc_sw],PH.inc_on_sw);
            PH.n_caps = n_caps;
            PH.n_off_sw = size(inc_sw,2);
            PH.n_loads = size(inc_loads,2);
            PH.inc_on_conv_sw = [supply_branch inc_caps PH.inc_on_sw inc_loads];
            PH.sw_idxs = idx_sw:(idx_sw-1)+PH.n_on_sw;
        end
        
        %Return incidence matrix only containing capacitor and load elements
        function conv = get_on_no_sw(PH)
            conv = PH.inc_on_conv(:,1:end-PH.n_off_sw);
        end
        
        %Return incidence matrix only containing capacitor and OFF-switch elements
        function conv = get_on_no_load(PH)
            conv = PH.inc_on_conv( :,[1:(PH.n_caps+1)...
                                    (PH.n_caps+PH.n_loads+2):end]);
        end
        
        %Return incidence matrix only containing capacitors and supply
        function conv = get_on_caps(PH)
            conv = PH.inc_on_conv( :,1:(PH.n_caps+1));
        end
        
        %Set a vector
        function set_a_vector(PH,a_vector)
            PH.a_vector=a_vector;
        end
        
        %Set a vector
        function set_a_vector_up(PH,a_vector_up)
            PH.a_vector = a_vector_up;
        end
        
        %Return a vector
        function a = get_a_vector(PH)
            a = PH.a_vector;
        end
        
        %Return a vector
        function a = get_a_vector_up(PH)
            a = PH.a_vector_up;
        end
        
        %Return ac vector --> a vector only for the caps
        function ac = get_ac_vector(PH)
            ac = PH.a_vector(2:end,:);
        end
        
        %Return b vector
        function b = get_b_vector(PH)
            b = PH.b_vector;
        end
        
        %Write the b vector
        function set_b_vector(PH,b_vector)
            PH.b_vector=b_vector;
        end
        
        %Return the b vector
        function r = get_r_vector(PH)
            %Check if exsit a or b vectors
            if ~isempty(PH.a_vector) && ~isempty(PH.a_vector) 
               if  isempty(PH.r_vector) 
                r =PH.get_ac_vector - PH.b_vector ;
                PH.r_vector = r;
               else
                  r =  PH.r_vector;
               end
            else
                error('SCC_Phase:vectors','a or b vector not initialized!');
            end
        end
        
        %create the switch element vector
        function ar_builder(PH)
          if  ~isempty(PH.a_vector)
                           
                %get the cut-set matrix including the on switches 
                %[Qx T] = fun_cutset(PH.inc_on_conv_sw,sort([1 11 12 10 6 9 8 4 2]));
                [Qx ] = fun_cutset(PH.inc_on_conv_sw);
                %Solve the current thorugh the switches
                %generate a vector with the row index that solves the the
                %currents in the switches 
                %If a switch is a twing this row must be in the index
                          
                %Rearrange the matrix  to the form
                %[Qsw | Qa | Qo] * [asw | a | I]'=0
                
                Qx=[Qx(:,(PH.n_caps+1) + (1:PH.n_on_sw)) ...
                    Qx(:,1:(PH.n_caps+1)) Qx(:,(end-(PH.n_loads-1)):end) ];
                
                %Arrange the matrix to echelon form
                Qx = rref(Qx);
                
                %2 first rows have the switches equations
                Qa = Qx(1:PH.n_on_sw,PH.n_on_sw + (1:PH.n_caps+1));
                Qo = Qx(1:PH.n_on_sw,(end-(PH.n_loads-1):end));
                
                %Compute the vector
                %Qx = Qsw(1:PH.n_caps+1);
                %Qo = Qsw(end-PH.n_loads:end);
                asw = -Qa*PH.a_vector - Qo;
                
                PH.ar_vector=[PH.get_ac_vector(); asw];
                
          else
              error('SCC_Phase:vectors','a or b vector not initialized!')
          end
          
        end
        
        function ar = get_ar_vector(PH)
            if isempty(PH.ar_vector)
                PH.ar_builder();
            end
            ar = PH.ar_vector;
        end
        
 
        
    end
    
end

function y = addSource(A,varargin)
Src_node = 1;
if nargin > 1
    Src_node=varargin{1};
end
Source_edge=zeros(size(A,1),1);
Source_edge(Src_node)=-1;
y = [Source_edge A];
end

function y = dependent_row(A,r_vec) %#ok<*DEFNU>

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
