function [Fo] = get_poles(top,Cx,Lx)
%% function Fo = get_poles(top,Cx,Lx)
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
for i_ph = 1:top.n_phases 
    %Get capacitor incidence matrix and shortcircuit the source 
    Ax = top.phase{i_ph}.get_on_caps;
    Ax =  phase_conv( Ax(:,2:end) , Ax(:,1));
    
    %The frequency must be computed for each capacitor 
    for i_elem = 1:size(Ax,2) 
        %First compute equivalent capacitor
        A_test = Ax(:,i_elem);
        A_elem = Ax;
        A_elem(:,i_elem)=[];
        Z_elem = Zc;
        Z_elem(i_elem) = [];
        Zeq             = solve_zeq([A_test A_elem],Z_elem);
        Ceq = Cx(i_elem)/(1+Cx(i_elem)*Zeq);
        Fo(i_ph,i_elem) = 1/(2*pi*sqrt(Ceq*Lx(i_elem)));
    end    
end

end