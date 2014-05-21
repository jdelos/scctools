function implementation = implement_topology(topology, Vin, switchTechs,...
    capTechs, compMetric)
% implement_topology: matches components with switches and components in
% the topology.
% implementation = implement_topology(topology, Vin, switchTechs,
%       capTechs, compMetric)
%   topology: structure created by generate_topology
%   Vin: input voltage of converter
%   switchTechs: an array of technology structures available for switch use
%   capTechs: an array of technology structures available for cap use
%   compMetric: a metric (1=area, 2=loss) used for determining the best
%               component (1=default)
% Created 4/15/08, Last Modified: 4/15/09
%   Copyright 2008-2009, Mike Seeman, UC Berkeley
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.

% Break out components of topology structure
ratio = topology.ratio;
ac = topology.ac;
ar = topology.ar;
vc = topology.vc*Vin;
vr = topology.vr*Vin;
vcb = topology.vcb*Vin;
vrb = topology.vrb*Vin;

if (nargin < 5)
    compMetric=1;
end

switch_assign = [];
cap_assign = [];
switch_rel_size = [];
cap_rel_size = [];

% Assign Capacitors
for i=1:size(ac,2),
    Mc = 0;
    Cc = 0;        % cap cost;
    for j=1:size(capTechs,2),
        if (vc(i) <= capTechs(j).rating),
            % Cap could work ... let's see if it's good
            % Use area-limited metric, which is usually applicable
            C = capTechs(j).area;
            M = capTechs(j).capacitance*vc(i)^2/C;
            if (M > Mc)
                if (Mc == 0),
                    cap_assign = [cap_assign capTechs(j)];
                else
                    cap_assign(i) = capTechs(j);
                end
                Mc = M;
                Cc = C;
            end
        end
    end
    % check to make sure a suitable device exists
    if (Mc == 0),
        error(strcat('No capacitors meet the voltage requirement of: ',...
            num2str(vc(i))));
    end
    % determine relative device size
    if (ac(i) == 0),
        cap_rel_size = [cap_rel_size 0]; % avoid divide by 0
    else
        cap_rel_size = [cap_rel_size (ac(i)*vc(i))/...
                (sqrt(Mc)*cap_assign(i).area)];
    end
    
end


% Assign Switches
for i=1:size(ar,2),
    Msw = 0;
    Csw = 0;        % switch cost;
    for j=1:size(switchTechs,2),
        if (vr(i) <= switchTechs(j).drain_rating),
            % Switch could work ... let's see if it's good
            if (compMetric == 2),   % loss metric
                % assume full gate drive
                C = switchTechs(j).gate_cap*switchTechs(j).gate_rating^2 + ...
                    switchTechs(j).drain_cap*vr(i)^2 + ...
                    switchTechs(j).body_cap*vrb(i)^2;
                M = switchTechs(j).conductance*vr(i)^2/C;
            else % area metric
                C = switchTechs(j).area;
                M = switchTechs(j).conductance*vr(i)^2/C;
            end
            if (M > Msw)
                if (Msw == 0),
                    switch_assign = [switch_assign switchTechs(j)];
                else
                    switch_assign(i) = switchTechs(j);
                end
                Msw = M;
                Csw = C;
            end
        end
    end
    % check to make sure a suitable device exists
    if (Msw == 0),
        error(strcat('No switches meet the voltage requirement of: ',...
            num2str(vr(i))));
    end
    % determine relative device size
    if (ar(i) == 0),
        switch_rel_size = [switch_rel_size 0];
    else
        if (compMetric == 2),
            switch_rel_size = [switch_rel_size (ar(i)*vr(i))/...
                    (sqrt(Msw)*switch_assign(i).conductance)];
        else
            switch_rel_size = [switch_rel_size (ar(i)*vr(i))/...
                    (sqrt(Msw)*switch_assign(i).area)];
        end
    end
end

% Scale Caps for unit area:

cap_area = 0;
for i=1:size(ac,2),
    cap_area = cap_area + cap_rel_size(i)*cap_assign(i).area;
end
cap_size = cap_rel_size./(cap_area+1e-30);

% Scale Switches for unit area:

sw_area = 0;
for i=1:size(ar,2),
    sw_area = sw_area + switch_rel_size(i)*switch_assign(i).area;
end
switch_size = (switch_rel_size.*(sw_area > 0))./(sw_area+(sw_area == 0));

%% Add-On J. Delos Philips Research
% Compute the alpha paremters for the SSL and FSL
% Calculate SSL output resistance:
caps = cap_assign;

Rssl_alpha = 0;
for i=1:size(caps,2),
    if (ac(i) > 0),
        Rssl_alpha = Rssl_alpha + (ac(i)^2)/(caps(i).capacitance*cap_size(i));
    end
end

% Calculate FSL output resistance:
switches = switch_assign;
sw_size = switch_size;
Rfsl_alpha = 0;
for i=1:size(switches,2),
    if (ar(i) > 0),
        Rfsl_alpha = Rfsl_alpha + 2*(ar(i)^2)/...
            (switches(i).conductance*sw_size(i));
    end
end


% Create implementation structure
implementation.topology = topology;
implementation.capacitors = cap_assign;
implementation.switches = switch_assign;
implementation.cap_size = cap_size;
implementation.switch_size = switch_size;
implementation.ssl_alpha = Rssl_alpha;
implementation.fsl_alpha = Rfsl_alpha;
