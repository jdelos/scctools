function [Fo,Qo] = get_poles(top,Cx,Lx,varargin)
%% function Fo = get_poles(top,Cx,Lx,[ Ron,Cesr] )
% Returns the list with the natural freqnecies produced by the parastic
% inductaonces in the capacitors. 
%
% Julià Delos julia.delos@philips.com
% 
% Philips Research 2015
% Copyrights 
%

Fo = zeros(top.n_phases,length(Lx));
Zc = 1./Cx;
Qo = -1;
q_flg = 0;
%%Check variable consitency
% If arrays are scalar all the values are considered the same size

if isscalar(Cx)
    Cx = ones(1,top.n_caps)*Cx;
elseif (top.n_caps ~= length(Cx) ) 
   error('Arugments:Parsing',...
'Cx is not cosistent with number of capacitors in the  the topolgy'); 
end

if isscalar(Lx)
    Lx = ones(1,top.n_caps)*Lx;
elseif   (top.n_caps ~= length(Lx) )
   error('Arugments:Parsing',...
'Lx is not cosistent with number of capacitors in the  the topolgy'); 
end


%%
switch nargin
    case 4  %On resistance has been specified
    q_flg = 1;
    Qo    = Fo;
    Ron = varargin{1};
    if isscalar(Ron)
        Ron = ones(1,top.n_switches)*Ron;
    elseif   (top.n_switches ~= length(Ron) )
        error('Arugments:Parsing',...
    'Ron is not cosistent with number of capacitors in the  the topolgy'); 
    end

    
    case 5  %On resistance and Cesr is provided
    q_flg = 1;
    Qo    = Fo;
    Ron = varargin{1};
    Cesr = varargin{2};
    
    if isscalar(Ron)
        Ron = ones(1,top.n_switches)*Ron;
    elseif   (top.n_switches ~= length(Ron) )
        error('Arugments:Parsing',...
    'Ron is not cosistent with number of capacitors in the  the topolgy'); 
    end
    
    if isscalar(Cesr)
        Cesr = ones(1,top.n_caps)*Cesr;
    elseif   (top.n_caps ~= length(Cesr) )
        error('Arugments:Parsing',...
    'Cesr is not cosistent with number of capacitors in the  the topolgy'); 
    end   
end        




for i_ph = 1:top.n_phases 
    
    %% Get capacitor incidence matrix and shortcircuit the source
    %This matrix is used to compute the equivalent capacitor 
    Ac_eq = top.phase{i_ph}.get_on_caps;
    Ac_eq =  phase_conv( Ac_eq(:,2:end) , Ac_eq(:,1));
    
     %% Compute the equivalent R
    % That is the joint Capacitor and switches incidence matrix
    %Build the incidence matrix with switches
    Ar_eq = phase_conv([top.inc_caps  top.inc_sw_act{i_ph}],top.supply_branch);  
    Rc = [Cesr Ron(logical(top.sw_activation_matrix(i_ph,:))) ];
    %The frequency must be computed for each capacitor 
    for i_elem = 1:size(Ac_eq,2) 
        %Compute equivalent capacitor
        Ac_test = Ac_eq(:,i_elem);
        Ac_elem = Ac_eq;
        Ac_elem(:,i_elem)=[];
        Z_elem = Zc;
        Z_elem(i_elem) = [];
        Zeq            = solve_zeq([Ac_test Ac_elem],Z_elem);
        Ceq = Cx(i_elem)/(1+Cx(i_elem)*Zeq);
        Fo(i_ph,i_elem) = 1/(2*pi*sqrt(Ceq*Lx(i_elem))); 
        
        if q_flg
            %Compute equivalent capacitor
            Ar_test = Ar_eq(:,i_elem);
            Ar_elem = Ar_eq;
            Ar_elem(:,i_elem)=[];
            R_elem = Rc;
            R_elem(i_elem) = [];
            Req             = solve_zeq([Ar_test Ar_elem],R_elem);
            Req = Cesr(i_elem)+Req;
            Qo(i_ph,i_elem) = 1/Req*sqrt(Lx(i_elem)/Ceq);
        end 
    end
end

end