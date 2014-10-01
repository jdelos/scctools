function [rx, minFoM ] = dickson_optimizer_fsl (topology,opt)
% dickson_optimizer_ssl(topology,mode,Io): Finds the optimal capacitor relative sizing fora a dkison
% topology. The otpimtzation is done for a converter operatin in the FSL
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
%	Philips Research, Eindhoven,  Netherlands
%  julia.delos@philps.com
%
    
%% Get mode value
if isfield(opt,'mode')
    mode = opt.mode;
else
    mode = -1;
end

%% Average duty cycle center of the optimitzation
if isfield(opt,'duty')
    duty = opt.duty;
else
    duty = -1;
end

%% Distribution function
if isfield(opt,'distr')
    distr = opt.distr;
else
    distr = 'normal';
end

%% Number of points used in the averaged optimitzation 
if isfield(opt,'dst_points')
    dst_points = opt.dst_points;
else
    dst_points = 30;
end

%% Weighted switch areas 
if isfield(opt,'rxua')
    rxua = opt.rxua/ opt.rxua(1); %Normalize rxua respect the first element
else
    rxua = ones(1,topology.N_sw);
end

%% Log display
if isfield(opt,'log')
    log = opt.log;
else
    log = 'off';
end

%% Averaged limits
if isfield(opt,'bounds')
    bounds = opt.bounds;
else
    bounds = [0.1 0.9];
end


%% Set optimitzation functions
options = optimset('fmincon');
options = optimset(options,'Algorithm','active-set', 'TolFun', 1e-11, 'MaxIter', 400000, ...
    'Display', log, 'LargeScale', 'off','MaxFunEvals',1400);

%% Get number of switches in the topology it corresponds to the total number of
% variables to be optimitzed
N =  topology.N_sw;


%% Generate Weighting Function
avgFoM = 0; 
if duty > 0 
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
    avgFoM = 1;
end

%% Get number of outputs
N_outs = length(topology.f_fsl);

if (min(mode) > 0) && (max(mode) <= N_outs )
    if isscalar(mode) 
        if ~avgFoM %Single point optimitzation
            FoM = @(x)double(subs(topology.f_fsl(mode),...
            symvar(topology.f_fsl),rxua./x));
        else %Average optimitzation
            FoM = sym(0);
            for j = 1:dst_points 
               FoM = FoM + subs(topology.f_fsl(mode),...
                symvar(topology.f_fsl),[smplX(j) sym('x',[1 topology.N_sw])])*weigD(j);
           end
           FoM = @(x)subs(FoM,symvar(FoM),rxua./x)/dst_points;
        end
    else
        mode = unique(mode);
        FoM = @(x)subs(sum(topology.f_fsl(mode)),...
            symvar(topology.f_fsl),rxua./x);
    end
else
    switch mode
        case -1
            FoM  = @(x)subs(sum(topology.f_fsl),symvar(topology.f_fsl),rxua./x);
            
        case -2
            FoM  = @(x)subs(sum(topology.f_fsl./(topology.ratio).^2),...
                symvar(topology.f_fsl),rxua./x);
        case -3 
            FoM  = @(x)subs(sum(topology.f_fsl.*Io.^2),...
                symvar(topology.f_fsl),rxua./x);
    end
end
x0(1:N)=1/N;
Aeq = ones(1,N);
Beq = 1;
A = -eye(N);
B = zeros([N 1]);

lb(1:N,1)=0.005;
ub = ones(1,N);

[rx, minFoM ] = fmincon(@(x)FoM(x), x0, A, B,Aeq,Beq,lb,ub,[],options);
if isfield(opt,'rxua')
    FoM = FoM*opt.rxua(1); % Bring normalitzation back
end
