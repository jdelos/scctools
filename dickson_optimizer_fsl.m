function [rx, FoM ] = dickson_optimizer_fsl (topology,mode, opt)
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
    
if isfield(opt,'rxua')
    rxua = opt.rxua;
else
    rxua = ones(1,length(symvar(topology.f_fsl)));
end

options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-11, 'MaxIter', 400000, ...
    'Display', 'on', 'LargeScale', 'off','MaxFunEvals',1400);

N = length(symvar(topology.f_fsl));


if (min(mode) > 0) && (max(mode) <= N )
    if isscalar(mode)
        FoM = @(x)subs(topology.f_fsl(mode),...
            symvar(topology.f_fsl),rxua./x);
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

[rx, FoM, ] = fmincon(@(x)FoM(x), x0, A, B,Aeq,Beq,lb,ub,[],options); 
