function [ al, m] = solve_charge_vectors(Q,n_caps,duty)
% [al, m ] = solve_charge_vectors(Q,n_caps)
%
% Function that returns the charge flow vectors from a multiphase SCC. Ql is a cell vector
% providing a the fundamental cut-set matrix for each phase period.
%
% A cut-set matrix Q satisfies KCL expressed for each circuit element charge qx:
%
%    Q . qx = 0 [1]
%
% The colummns of the cut-set matrix is organised as follows  [Q_in | Q_c | Q_o ]
% where:
%	 Q_in --> Input supply (colummne vector)
% 	 Q_c --> Capacitor cut-set sub-matrix
% 	 Q_o --> Output load cut-set sub-matrix
%
% All the capacitors satisfies the charge-flow balance when reached the Steady State,
% thus expressed as
%	 q_c_1 + q_c_2 + ... + q_c_(n-1) +q_c_n = 0 [2]
% 	where 1,2,..., n-1, n, are phases.
%
% Using the eq. 1 and 2 the charge flow vectors are obtained by solving
%  	al = -[ Q_in | Q_c ] \ Q_o
%
% The function returns
%  	al --> Charge flow cell array, on cell for each phase
%	m -->  Converter conversion ratio vector, one element for each output
%
%
%Copyright 2013-2014, Julia Delos, Philips Research
%	julia.delos@philips.com
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.


%% Get number of phases
n_phase = length(Q);

%% Create system matrix Qx = [Q_in | Q_c]
m_size = n_phase*(1+n_caps);
Qx = zeros(m_size);
n_outs = size(Q{1},2) - (1+n_caps); %Number of output nodes

if isempty(symvar(duty(1)))
    Qo = zeros(m_size,n_outs);
else
    Qo = sym(zeros(m_size,n_outs));
end

%Cut-set matrix first element is the voltage supply, then the capacitors in
%conscutively order thus first row is skiped

cap_cols=2:n_caps+1;

%% Populate Qx(j_idx,i_idx) matrix
%Init index
j_idx = 1;
i_idx = 1;
for p =1:n_phase
    Qp = Q{p}; %Get phase matrix
    n_rows = size(Qp,1); %Get number of rows
    Qp(:,1) = -Qp(:,1); % Input supply always delivers charge (sign is changed)
    %Add cut-set matrix
    Qx(j_idx+(0:(n_rows-1)), i_idx + (0:(n_caps))) = Qp(:,[1 cap_cols]);
    %Gather output matrix
    Qo(j_idx+(0:(n_rows-1)), : ) = Qp(:, (n_caps+2):end );
    
    %Add charge flow balance equations
    Qx((m_size-(n_caps-1)):m_size,(i_idx-1)+cap_cols) = eye(n_caps);
    
    j_idx = j_idx +  n_rows;
    i_idx = i_idx + (n_caps+1);
end

%%Solve de linear system
ax = -Qx\Qo;

%%Rearrange the charge vector by phases
al = cell(1,n_phase);
j_idx = 1; 
m = zeros(1,n_outs);
for p=1:n_phase
    %Charge matrix per phase
    al{p} = ax(j_idx +(0:n_caps),:);
    
    %Collec input vector
    m = m + ax(j_idx,:);
    j_idx = j_idx + (n_caps+1);
end
%Return a column vector
m=m.';
end





