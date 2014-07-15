function [cx, FoM ] = dickson_optimizer_ssl (topology,mode,cond)
% dickson_optimizer_ssl(topology,mode,Io): Finds the optimal capacitor relative sizing fora a dkison
% topology. The otpimtzation is done for a converter operatin in the ssl
%   
%   Based on Michael Seeman work
%Input paramters: 
%   + topology: dikson topology structure
%   + mode: -1: All outputs Constant Output Current FoM
%           -2: All outputs Constant Input Power FoM
%           -3: All ouputs weighted with the current vector Io
%     scalar N : Optmized ofr the Nth output 
%     vector N : Optimitzed ofr the Nth outputs 
%
%   Returns: performance structure from evaluate_loss and optimal switching
%   frequency and switch area
%
%   Created 4/14/08, Last modified: 4/15/09
%   Copyright 2008-2009, Mike Seeman, UC Berkeley
%   Copyright 2013-2014, J. Delos, Philips Research Netherlands
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.
%
%  Philips Research, Eindhoven,  Netherlands
%  julia.delos@philps.com
%
    
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-11, 'MaxIter', 400000, ...
    'Display', 'on', 'LargeScale', 'on','Algorithm','active-set');

%% Get number of capacitors
N = length(symvar(topology.f_ssl));

%% Get number of outputs
N_outs = length(topology.f_ssl);

%% Standard Optmitzation paramteters
flg = -10;
lb_t=0.01;
itr=0;
fcond = [];

x0(1:N)=1/N;
Aeq = ones(1,N);
Beq = 1;
A = -eye(N);
B = zeros([N 1]);

ub = ones(1,N);
lb(1,1:N)=lb_t;

if (min(mode) > 0) && (max(mode) <= N_outs )
    if isscalar(mode)
        FoM = @(x)subs(topology.f_ssl(mode),...
            symvar(topology.f_ssl),x);
    else
        mode = unique(mode);
        FoM = @(x)subs(sum(topology.f_ssl(mode)),...
            symvar(topology.f_ssl),x);
    end
else
    switch mode
        case -1
            FoM  = @(x)subs(sum(topology.f_ssl),symvar(topology.f_ssl),x);
            
        case -2
            FoM  = @(x)subs(sum(topology.f_ssl./(topology.ratio).^2),...
                symvar(topology.f_ssl),x);
        case -3 
            FoM  = @(x)subs(sum(topology.f_ssl.*cond(1).^2),...
                symvar(topology.f_ssl),x);
        case -4 %DC node optimitzation with restricted output ripple
            fsw = cond(4);
            
            f  = @(x)subs(topology.f_ssl(topology.dc_outputs),...
                symvar(topology.f_ssl),x);
            FoM =@(x) f(x(2:end))/(x(1)*fsw); %Rssl
            
            fcond = @(x)rippCond(x,topology,cond);
            
    end
end

%% Omptimitze iterations
while(flg < 0 && itr <= 5 )

if isempty(fcond)
    [cx, FoM, flg] = fmincon(@(x)FoM(x), x0, A, B,Aeq,Beq,lb,ub,...
    [],options); 

else
    %% DC Node optimitzation
    Ct_o = 10e-9;
    Ct_up = 5e-3;
    Ct_lb = 1e-9;
    
    x0 = [Ct_o ones(1,N)*1/N];
    %lb=[ lb*ones(1:N,1)];
    Aeq =[ 0 ones(1,N)];
    Beq = 1;
    A = -eye(N+1);
    B = zeros([N+1 1]);
    ub = [Ct_up ones(1,N)];
    lb = [Ct_lb ones(1,N)*lb_t];
    [cx, FoM, flg] = fmincon(@(x)FoM(x), x0, A, B,Aeq,Beq,lb,ub,...
    @(x)fcond(x),options); 

end
itr = itr +1 ;
disp(flg)
if any(lb == cx) 
 flg = -1;
 idxs = (lb == cx);
 lb(idxs)=lb(idxs)*.1; %% Update lower boundareis

 disp(['Low bounderi reached changed to:' num2str(lb) ])
end
break
end

function [c , ceq] = rippCond(cx,topology,cond)
    
    dc_node = topology.dc_outputs;
    
    Vr_pp = cond(1);
    eff = cond(2);
    Ro = cond(3);
    fsw= cond(4);
    
    Ct = cx(1);
    
    %% Get Rssl from the toplogies
    Rssl = topology.eval_ssl(cx(2:end));
    Rssl = Rssl(dc_node)/(fsw*Ct);
    %% Get redistribuited output charge
    qc_dc = max(topology.eval_q_dc(cx(2:end)));
    
    c = 0;
    eq1 = (qc_dc/Ct*cx(dc_node+1) - Vr_pp) ; %% Ripple restriction
    eq2 = eff - Ro/(Ro+Rssl); %% Efficency restriction
    ceq = eq1 +eq2;
    
    