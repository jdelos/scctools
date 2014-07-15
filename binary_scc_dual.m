%% Scrip to generate the incidence matrix from a multiphase converter with
% 3 capacitors
%
%   Copyright 2014, Julia Delos, Philips Research 
%	julia.delos@philips.com	
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work

%% Architecture converter  matrices
% Capacitors Architecture
Acaps = zeros(6,5);
Acaps([2 6],1) = [1 -1]; % C1 
Acaps([3 7],2) = [1 -1]; % C2
Acaps([4 8],3) = [1 -1]; % C3
Acaps(5 ,4)    = 1;      % Co_I
Acaps(9 ,5)    = 1;      % Co_II

% Switch architecture
Aswitches = zeros(6,16);
Aswitches([1 2], 1 )  = [1 -1]; % S1
Aswitches([2 3], 2 )  = [1 -1]; % S2
Aswitches([2 7], 3 )  = [1 -1]; % S3
Aswitches([3 6], 4 )  = [1 -1]; % S4
Aswitches([6 7], 5 )  = [1 -1]; % S5
Aswitches(  6  , 6 )  = 1 ;    % S6
Aswitches([3 4], 7 )  = [1 -1]; % S7
Aswitches([3 8], 8 )  = [1 -1]; % S8
Aswitches([4 7], 9 )  = [1 -1]; % S9
Aswitches([7 8], 10 ) = [1 -1]; % S10
Aswitches([4 5], 11 ) = [1 -1]; % S11
Aswitches([8 9], 12 ) = [1 -1]; % S12
Aswitches(   8 , 13 ) = 1     ; % S13
Aswitches([1 6], 14 ) = [1 -1]; % S14
Aswitches([4 9], 15 ) = [1 -1]; % S15
Aswitches([5 8], 16 ) = [1 -1]; % S16
%% Codes for 1/8


SW785.exb_code  = [[1  0  1  0]
                   [1  1 -1  1]
                   [1  0  1  1]
                   [1  1  0  0]
                   [1 -1 -1  0]
                   ];
                     
SW785.switches  = [ 
[1 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 ]
[0 1 0 0 0 0 0 0 0 1 1 0 0 1 0 0 ]
[1 0 1 0 0 0 0 1 0 0 0 0 0 0 1 0 ]
[0 1 0 0 0 0 1 0 0 0 0 0 0 1 1 0 ]
[1 0 0 1 0 0 0 0 0 1 0 0 1 0 0 0 ]
];

SW785.ratio = [7/5 8/5];
SW785.phases = size(SW785.exb_code,1);
for j = 1:SW785.phases
   SW785.Asw{j} = Aswitches(:,logical(SW785.switches(j,:)));
   sw_name = {};
   idx = 1;
   for x = find(SW785.switches(j,:))
     sw_name{idx} =  ['S' num2str(x)];
     idx = 1+ idx;
   end
   SW785.Ref_sw{j} = sw_name ;
end

%% Generate to the topology
Ron = 1.2;
Cfly = 4.7e-6;
Cdc = 470e-6;

Fsw = 1/(5*5e-6);

%Generate the topology
t = generic_switched_capacitor_class(Acaps,...
    SW785.Asw{1},SW785.Asw{2},SW785.Asw{3},SW785.Asw{4},SW785.Asw{5},...
    'Duty',[1/5 1/5 1/5 1/5]); 

k_fsl = subs(t.k_fsl,symvar(t.k_fsl), [0 0 0 0 0 ones(1,8)*Ron]);
k_ssl = subs(t.k_ssl,symvar(t.k_ssl), [Cfly Cfly Cfly Cdc Cdc]);

%%  Compute the entire matrix 
k_scc = sqrt((k_ssl./Fsw).^2 + k_fsl.^2);

%% Sumbstract the matrix for the DC nodes
dc_nodes = [4 8];
k = k_scc(dc_nodes,dc_nodes);





