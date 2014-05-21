function performance = evaluate_loss (imp, Vin, Vout, Iout, fsw, Asw, Ac)
% evaluate_loss: evaluates the loss and other performance metrics for a 
% specific size and operating condition of a implemented SC converter
%
%   imp: implementation generated from implement_topology
%   Vin: converter input voltage for this calc [V]
%   Vout: converter output voltage for this calc [V]
%   Iout: converter output current for this calc [A]
%   fsw: switching frequency [Hz]
%   Asw: switch area [m^2]
%   Ac: capacitor area [m^2]
%
%   Either fsw (in case of regulation) or Vout (in case of open-loop
%   control) should be set to [] and will be found by this function.
%
%   Created: 4/15/08, Last Modified: 4/15/09
%   Copyright 2008-2009, Mike Seeman, UC Berkeley
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.
%
%  Philips Research, Eindhoven,  Netherlands
%  julia.delos@philps.com
%
    


% Break implementation into components for brevity
ac = imp.topology.ac;
vc = imp.topology.vc;
vcb = imp.topology.vcb;
ar = imp.topology.ar;
vr = imp.topology.vr;
vrb = imp.topology.vrb;
ratio = imp.topology.ratio;

caps = imp.capacitors;
cap_size = imp.cap_size;
switches = imp.switches;
sw_size = imp.switch_size;

% Process (and expand) input parameters
paramdim = max(size(Vin), max(size(Vout), max(size(Iout), max(size(fsw),...
    max(size(Asw), size(Ac))))));
type = 0;   % undefined: 0-not yet set; 1-vout; 2-fsw

Vin = expand_input(Vin, paramdim);
if (size(Vout) == [0 0])
    type = 1;
    Vout = zeros(paramdim);
else
    Vout = expand_input(Vout, paramdim);
end
Iout = expand_input(Iout, paramdim);
if (size(fsw) == [0 0])
    if (type == 1)
        error('Both fsw and Vout cannot be undefined');
    else
        type = 2;
        fsw = zeros(paramdim);
    end
else
    fsw = expand_input(fsw, paramdim);
end
Asw = expand_input(Asw, paramdim);
Ac = expand_input(Ac, paramdim);

%%%%%%%%%%%%%%%%%%%%%%%%%% Start Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate SSL output resistance:

Rssl_alpha = 0;
for i=1:size(caps,2),
    if (ac(i) > 0),
        Rssl_alpha = Rssl_alpha + (ac(i)^2)/(caps(i).capacitance*cap_size(i));
    end
end

% Calculate FSL output resistance:

Rfsl_alpha = 0;
for i=1:size(switches,2),
    if (ar(i) > 0),
        Rfsl_alpha = Rfsl_alpha + 2*(ar(i)^2)/...
            (switches(i).conductance*sw_size(i));
    end
end
Rfsl = Rfsl_alpha./Asw;

% Calculate additional ESR loss:
Resr_alpha = 0;
for i=1:size(caps,2),
    if (isfield(caps(i),'esr')),
        if (ac(i) > 0),
            Resr_alpha = Resr_alpha + 4*(ac(i)^2)*(caps(i).esr/cap_size(i));
        end
    end
end
Resr = Resr_alpha./Ac;
if (isfield(imp, 'esr')),
    Resr = Resr + imp.esr;
end

% Find the unknown (endogenous) variable
if (type == 1)
    % Vout is unknown
    Rssl = Rssl_alpha ./ (fsw.*Ac);
    
    % Calculate total output resistance:
    Rout = sqrt(Rssl.^2 + (Rfsl + Resr).^2);
    Vout = Vin*ratio - Rout.*Iout;
    Pout = Vout.*Iout;
    is_prac = ones(paramdim);
elseif (type == 2)
    % fsw is unknown
    % Calculate needed output impedance and switching frequency to match
    % output voltage
    % is_prac is 1 if a finite fsw which satisfies Iout, Vin, Vout exists
    
    Rreq = (Vin*ratio - Vout)./Iout;
    is_prac = ((Rreq > 0) & (Rfsl+Resr < Rreq));
    Rssl = real(sqrt(Rreq.^2 - (Rfsl+Resr).^2));
    fsw = Rssl_alpha ./ (Rssl .* Ac);
    
    % Calculate total output resistance:
    Rout = sqrt(Rssl.^2 + (Rfsl + Resr).^2);
    Pout = Vout.*Iout;
else
    error('Either Vout or fsw must be []');
end

% Calculate resistance losses:
Pssl = Rssl.*Iout.^2;
Pfsl = Rfsl.*Iout.^2;
Pesr = Resr.*Iout.^2;
Pres = Rout.*Iout.^2;

% Calculate Cap-related Parasitic loss:

Pc_alpha = 0;
for i=1:size(caps,2),
    Pc_alpha = Pc_alpha + (caps(i).bottom_cap*cap_size(i)) * (vcb(i))^2;
end
Pc = Pc_alpha.*fsw.*Ac.*Vin.^2;

% Calculate Switch-related Parasitic loss:

Psw_alpha = 0;
Pg_alpha = 0;
for i=1:size(switches,2),
    % Assume switch is driven at full gate_rating voltage
    Vgssw = switches(i).gate_rating;
    Pg_alpha = Pg_alpha + (switches(i).gate_cap*sw_size(i)) * (Vgssw)^2;
    Psw_alpha = Psw_alpha + ((switches(i).drain_cap*sw_size(i)) * (vr(i))^2 + ...
                             (switches(i).body_cap*sw_size(i)) * (vrb(i))^2);
end
Psw = (Psw_alpha.*Vin.^2 + Pg_alpha).*fsw.*Asw;

% Calculate total loss, efficiency, etc...
Ploss = Pres + Pc + Psw;
eff = Pout./(Pout+Ploss);
IND = [];
for i=1:size(Pssl,1),
    [temp, ind] = max([Pssl(i,:); Pfsl(i,:); Pesr(i,:); Pc(i,:); Psw(i,:)]);
    IND = [IND; ind];
end
%[temp, ind] = max([Pssl Pfsl Pc Psw]);
texts = [{'SSL Loss'} {'FSL Loss'} {'ESR Loss'} {'Bottom-Plate'}...
        {'Switch Parasitic'}];

% create structure -- list unknown system parameters
performance.Vout = Vout;
performance.fsw = fsw;
performance.is_possible = is_prac;

% list performance
performance.efficiency = eff;
performance.total_loss = Ploss;
performance.impedance = Rout;
performance.dominant_loss = IND;
performance.dominant_text = texts(IND);

% Sub-function that turns input into a (maxsize) matrix
function result = expand_input (input, maxsize)

if (size(input) == [1 1])
    % scalar input
    result = input * ones(maxsize);
elseif (size(input) == [1 maxsize(2)])
    % row vector input
    result = ones(maxsize(1),1) * input;
elseif (size(input) == [maxsize(1) 1])
    % column vector input
    result = input * ones(1, maxsize(2));
elseif (size(input) == maxsize)
    % input already a properly-sized matrix
    result = input;
elseif (size(input) == [0 0])
    error('Only fsw or Vout can be empty');
    result = 0;
else
    % input is not properly sized
    error('All inputs must have the same number of rows and columns (if not 1)');
    result = 0;
end