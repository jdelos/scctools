function result = generate_topology (topology_name, num, den)
% generate_topology: Topology charge and voltage multiplier generation
%
%   result = generate_topology (topology_name, num [, den])
%       topology_name: Name of SC topology {Ladder, Dickson, Cockcroft-Walton,
%           Doubler, Series-Parallel, Fibonacci}
%       num, den: the ratio of output voltage to input voltage.  If denom is
%           omitted, num is a rational number.  Otherwise, both are integers.
%           num/den is the step-up conversion ratio of the converter.
%       
%       result: structure with charge multiplier and voltage vectors for
%           capacitors and switches.  
%   
% Created 4/14/08, Last Modified 4/15/09
%   Copyright 2008-2009, Mike Seeman, UC Berkeley
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.

% if a straight ratio passed in, break it up into numerator and denominator
if nargin < 3,
    [num, den] = rat(num, .001);
end

% generate ac's for step-up only, at first, v's in terms of input voltage
flip = 0;
if (num/den < 1),
    t = den;   den = num;  num = t;
    flip = 1;   % indicate that the converter has been 'flipped'
end

%%%%%%%%%%%%%%%%%% Series-Parallel %%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(topology_name,'Series-Parallel')),
    n = num;
    m = den;
    N = n;
    % SSL values
    ac = ones(1,m*(n-m))/m;
    vc = ones(1,m*(n-m))/m;
    vcb = [];
    for i=1:m,
        for j=1:(n-m),
            vcb = [vcb (i+j-1)/m];
        end
    end
    % FSL values
    vr = [];
    vrb = [];
    for i = 1:m,
        for j = 1:(n-m+1),
            if (j == 1),
                vr = [vr i/m];
                vrb = [vrb (i+j-1)/m];
            elseif (j == (n-m+1)),
                vr = [vr (n-m-1+i)/m];
                vrb = [vrb (i+j-2)/m];
            else
                vr = [vr 1/m];
                vrb = [vrb (i+j-1)/m];
            end
        end
    end
    for i=1:m+1,
        for j=1:n-m,
            if (i==1),
                vr = [vr j/m];
            elseif (i==(m+1)),
                vr = [vr (m-1+j)/m];
            else
                vr = [vr 1/m];
            end
            if ((i==1) || (i==(m+1))),
                vrb = [vrb 0];
            else
                vrb = [vrb (i+j-2)/m];
            end
        end
    end
    ar = ones(size(vr))/m;
    
%%%%%%%%%%%%%%%%%% Ladder %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmpi(topology_name,'Ladder')),
    n = num;
    m = den;
    N = n;
    % SSL values
    ac = [];
    vc = [];
    vcb = [];
    
    for j=1:(n-m-1),
        ac = [ac j j];
        vc = [vc 1/m 1/m];
        vcb = [vcb 1/m 0];
    end
    ac = [ac (n-m)];
    vc = [vc 1/m];
    vcb = [vcb 1/m];
    for j=(m-1):-1:1,
        ac = [ac ones(1,2)*(j*(n/m-1))];
        vc = [vc 1/m 1/m];
        vcb = [vcb 1/m 0];
    end
    
    % FSL values
    ar = [ones(1,2*(n-m)) (n/m-1)*ones(1,2*m)];
    vr = ones(1,2*n)/m;
    vrb = mod(0:(2*n-1),2)/m;
    
%%%%%%%%%%%%%%%%%% Dickson %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmpi(topology_name,'Dickson')),
    if (den ~= 1)
        error('SWITCHCAP:nonIntegerRatio',...
            'the Dickson topology supports integer ratios only');
    end
    N = num;
    % SSL values
    ac = ones(1,N-1);
    vc = [];
    vcb = ones(1,N-1);
    for j=1:(N-1),
        vc = [vc j];
    end
    % FSL values
    if (N == 2),
        vr = ones(1,4);
        ar = ones(1,4);
        vrb = [0 1 0 1];
    else,
        vr = [ones(1,5) 2*ones(1,N-2), 1];
        ar = [floor((j+1)/2)*[1 1] floor((j)/2)*[1 1] ones(1,N)];
        vrb = [0 1 0 1 ones(1,N)];
    end
    
%%%%%%%%%%%%%%%%%% Cockcroft-Walton %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmpi(topology_name, 'Cockcroft-Walton')),
    if (den ~= 1)
        error('SWITCHCAP:nonIntegerRatio',...
            'the Cockcroft-Walton topology supports integer ratios only');
    end
    N = num;
    % SSL values
    ac = [];
    vc = [1 2*ones(1,N-2)];
    vcb = ones(1,N-1);
    for j=1:(N-1),
        ac = [floor((j+1)/2) ac];
    end
    % FSL values
    if (N == 2),
        vr = ones(1,4);
        ar = ones(1,4);
        vrb = [0 1 0 1];
    else,
        vr = [ones(1,5) 2*ones(1,N-2) 1];
        ar = [floor((j+1)/2)*[1 1] floor((j)/2)*[1 1] ones(1,N)];
        vrb = [0 1 0 1 ones(1,N)];
    end
    
%%%%%%%%%%%%%%%%%% Doubler %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmpi(topology_name,'Doubler')),
    if (den ~= 1)
        error('SWITCHCAP:nonIntegerRatio',...
            'the Doubler topology supports integer ratios only');
    end
    n = ceil(log2(num));
    N = 2^n;
    if (N ~= num),
        error('SWITCHCAP:badRatio',...
            'the doubler topology supports conversion ratios ~ 2^n');
    end
    
    % SSL values
    ac = [];
    vc = [];
    vcb = [];
    for j=1:(2*n-1),
        ac = [2^floor((j-1)/2) ac];
        vc = [vc 2^floor(j/2)];
        vcb = [vcb 2^floor(j/2)*mod(j,2)];
    end
    % FSL values
    ar = [];
    vr = [];
    vrb = [];
    for j=1:n,
        ar = [ar ones(1,4)*2^(j-1)];
        vr = [vr ones(1,4)*2^(n-j)];
        vrb = [vrb [0 1 0 1]*2^(n-j)];
    end
    
%%%%%%%%%%%%%%%%%% Fibonacci %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmpi(topology_name,'Fibonacci')),
    if (den ~= 1)
        error('SWITCHCAP:nonIntegerRatio',...
            'the Fibonacci topology supports integer ratios only');
    end
    i = 2;
    while (fibfun(i) < num),
        i = i+1;
    end
    if (fibfun(i) > num)
       error('SWITCHCAP:badRatio',...
           'the fibonacci topology supports ratios of F_n or 1/F_n only');
    end
    N = fibfun(i);
    
    % SSL Calculation
    ac = [];
    vc = [];
    vcb = [];
    for j = 2:i-1,
        ac = [fibfun(j-1) ac];
        vc = [vc fibfun(j)];
        vcb = [vcb fibfun(j-1)];
    end
    % FSL Calculation
    ar = [1];
    vr = [];
    vrb = [0];
    for j = 2:i-1,
        ar = [fibfun(j) fibfun(j-1) fibfun(j-1) ar];
        vr = [vr fibfun(j) fibfun(j) fibfun(j-1)];
        vrb = [vrb fibfun(j-1) 0 fibfun(j-1)];
    end
    vr = [vr fibfun(i-2)];
end

ratio = num/den;

% Compute ideal metrics
Mssl = 2*ratio^2/(ac*vc')^2;
Mfsl = ratio^2/(2*(ar*vr')^2);

if (flip == 1),
    ac = ac/ratio;
    vc = vc/ratio;
    vcb = vcb/ratio;
    ar = ar/ratio;
    vr = vr/ratio;
    vrb = vrb/ratio;
    ratio = 1/ratio;
end

result.topName = topology_name;
result.ac = ac;
result.vc = vc;
result.vcb = vcb;       % bottom plate voltage vector
result.ar = ar;
result.vr = vr;
result.vrb = vrb;       % body swing voltage, for NMOS
result.Mssl = Mssl;
result.Mfsl = Mfsl;
result.ratio = ratio;