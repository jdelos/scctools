%% Scrip to generate the incidence matrix from a multiphase converter with
% 3 capacitors
%
%   Copyright 2014, Julia Delos, Philips Research 
%	julia.delos@philips.com	
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work

%% Architecture converter  matrices
% Capacitors Architecture
Acaps = zeros(6,4);
Acaps([2 5],1) = [1 -1]; % C1 
Acaps([3 6],2) = [1 -1]; % C2
Acaps([4 7],3) = [1 -1]; % C3
Acaps(  5  ,4) = 1;      % Co

% Switch architecture
Aswitches = zeros(6,12);
Aswitches([1 2],1) = [1 -1]; % S1
Aswitches([2 3],2) = [1 -1]; % S2
Aswitches([2 6],3) = [1 -1]; % S3
Aswitches([5 3],4) = [1 -1]; % S4
Aswitches([5 6],5) = [1 -1]; % S5
Aswitches(5,6)     = 1;      % S6
Aswitches([3 4],7) = [1 -1]; % S7
Aswitches(3,8)     =  1;     % S8
Aswitches([6 4],9) = [1 -1]; % S9
Aswitches([6 7],10)= [1 -1]; % S10
Aswitches([4 5],11)= [1 -1]; % S11
Aswitches([5 7],12)= [1 -1]; % S12

%% Codes for 1/8

SW18.exb_code  = [ [0  0  0  1]
                   [0  0  1 -1]
                   [1 -1 -1 -1]
                   [0  1 -1 -1]
                   ];
                     
SW18.switches  = [ 
[0 0 0 0 1 1 0 0 0 1 1 0]
[0 0 0 0 1 1 1 0 0 0 0 1]
[1 0 0 1 0 0 0 0 1 0 0 1]
[0 1 0 0 0 1 0 0 1 0 0 1]
];

SW18.ratio = 1/8;
SW18.phases = size(SW18.exb_code,1);
for j = 1:SW18.phases
   SW18.Asw{j} = Aswitches(:,logical(SW18.switches(j,:))) ;  
end





