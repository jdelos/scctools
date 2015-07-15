function [cx, minFoM ] = dickson_optimizer_ssl (topology,opt)
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


%% Get mode value
if isfield(opt,'mode')
    mode = opt.mode;
else
    mode = -1;
end

%% Get ripple output condition [NOT WORKING]
if isfield(opt,'cond')
    cond = opt.cond;
end

%% Average duty cycle center of the optimitzation
if isfield(opt,'duty')
    duty = opt.duty;
else
    duty = -1;
end

%% Number of points used in the averaged optimitzation 
if isfield(opt,'dst_points')
    dst_points = opt.dst_points;
else
    dst_points = 30;
end

%% Distribution function
if isfield(opt,'distr')
    distr = opt.distr;
else
    distr = 'flat';
end

%% Averaged limits
if isfield(opt,'bounds')
    bounds = opt.bounds;
else
    bounds = [0.1 0.9];
end

%% Averaged optimitzation
if isfield(opt,'averaged')
    avgFoM = opt.averaged;
else
    avgFoM = false; 
end


%% Optimitzation display 
if isfield(opt,'log')
    log = opt.log;
else
    log = 'off';
end


%% Configure optimizer
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-11, 'MaxIter', 400000, ...
    'Display', log, 'LargeScale', 'on','Algorithm','active-set');


%% Get number of capacitors
N = topology.N_caps;

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

%% Generate Weighting Function

if  avgFoM  
    smplX = linspace(bounds(1),bounds(2),dst_points);
    x     =  -1:0.001:1;
    sig   = 1/sqrt(2*pi);
    switch (distr)
        case 'normal'
            norm  = normpdf(x,0,sig);
            weigD = interp1(x+duty,norm,smplX);
        case 'flat'
            weigD = smplX*0+1;           
        case 'cos'
            norm  = cos(pi*x);
            weigD = interp1(x+duty,norm,smplX);
    end
end

if (min(mode) > 0) && (max(mode) <= N_outs )
    if isscalar(mode)
        if ~avgFoM %Single point optimitzation
            FoM = @(x)subs(topology.f_ssl(mode),...
                symvar(topology.f_ssl),x);
        else %Averaged optimitzation
           FoM = sym(0);
           for j = 1:dst_points 
               FoM = FoM + subs(topology.f_ssl(mode),...
                symvar(topology.f_ssl),[sym('x',[1 topology.N_caps]) smplX(j)])*weigD(j);
           end
           FoM = @(x)eval(subs(FoM,symvar(FoM),x)/dst_points);
        end     
    else
        mode = unique(mode);
        if ~avgFoM %Single point optimitzation
            FoM = @(x)subs(sum(topology.f_ssl(mode)),...
                symvar(topology.f_ssl),x);
        else
           FoM = sym(0);
           for j = 1:dst_points 
               FoM = FoM + subs(sum(topology.f_ssl(mode)),...
                symvar(topology.f_ssl),[sym('x',[1 topology.N_caps]) smplX(j)])*weigD(j);
           end
           FoM = @(x)eval(subs(FoM,symvar(FoM),x)/dst_points); 
        end
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

%% Launch optimitzation
if isempty(fcond)
    [cx, minFoM, flg] = fmincon(@(x)FoM(x), x0, A, B,Aeq,Beq,lb,ub,...
    [],options); 
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
    
    