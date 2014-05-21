function [cx, FoM ] = dickson_optimizer (topology,mode,Io,op_mode)
% optimize_loss: Finds the optimal capacitor relative sizing fora a dkison
% topology
%   
%   Based on Michael Seeman work
%
%   dickson_optimizer(topo)
%       topology: dikson topology
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
    'Display', 'off', 'LargeScale', 'off','MaxFunEvals',2500);

N = length(symvar(topology.f_ssl));


if (min(mode) > 0) && (max(mode) <= N )
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
            FoM  = @(x)subs(sum(topology.f_ssl.*Io.^2),...
                symvar(topology.f_ssl),x);
    end
end
x0(1:N)=1/N;
Aeq = ones(1,N);
Beq = 1;
A = -eye(N);
B = zeros([N 1]);

lb(1:N,1)=0.01;
ub = ones(1,N);

[cx, FoM] = fmincon(@(x)FoM(x), x0, A, B,Aeq,Beq,lb,ub); 
