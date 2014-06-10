function [ a1, a2, m] = solve_2ph_charge_vectors(Qf1,Qf2,n_caps,mode,duty)
%SOLVE_2PH_CHARGE_VECTORS 
%In steady state the charge flow in the capacitors of ph1 is equal and 
%opposite to the ph2. To solve the system the charge flow of capacitors
%is changed to opposite singe. Incidence matrix are written defining the
%direction of the vectors corresponding to the standard current flow of a
%capacitor. 
%The charge flow analysis considers negative flows the ones that come out
%of the capacitor, so as a negative current

%the cut-set matrix corresponds the charge flows during a phase period as
%   [q_in1 + q_c1 + q_c2 + ... + q_cn + q_out_o + q_out_1 ..+  q_out_m] = 0      
%
%Copyright 2013-2014, Julia Delos, Philips Research 
%	julia.delos@philips.com	
%   May be freely used and modified but never sold.  The original author
%   must be cited in all derivative work.



%Cut-set matrix first element is the voltage supply, then the capacitors in
%conscutively order thus first row is skiped  
%
cap_cols=2:n_caps+1;

%Charge flow in the capacitors in Steady State qcx_ph1 = -qcx_ph2
Qf2(:,cap_cols) = -Qf2(:,cap_cols);

%Current flow out of the voltage supply has been defined like a
%capacitor charges since it is delivering charge singes must be changed
Qf1(:,1) = -Qf1(:,1); %Ph1
Qf2(:,1) = -Qf2(:,1); %Ph2


%Since during each phase the contribution of the input source is different 
% a extra zero colum is added in each cut-set matrix

%First column correpsonds to ph1 and Second column correpsonds to ph2 

if mode
    A = [Qf1; Qf2];
    a1=-A(:,1:n_caps+1)\A(:,(n_caps+2):end);
    m = a1(1,:).';
    
    ac = a1(2:end,:);
    a1 = [a1(1,:).*duty(1); ac];
    a2 = [a1(1,:).*duty(2); -ac]; 
else
    %No charge flowing of ph2 source in Qf1
    Qf1=[Qf1(:,1) zeros(size(Qf1,1),1) Qf1(:,2:end)];

    %No charge flowing of ph1 source in Qf2
    Qf2=[zeros(size(Qf2,1),1) Qf2];

    %A charge vector flow matrix can be created staking Qf1 and Qf2
    A=[Qf1;Qf2];
    %Matrix can be spplited in variables and ouput exictations as 
    % A = [Ax | Axo] so the solution of the charge flow vectors is 
    % a= -Ax\Axo

    a1= -A(:,1:n_caps+2)\A(:,n_caps+3:end);

    %Charge flow vectors belonging to the capacitors
    ac=a1(3:end,:);

    a2= [a1(2,:) ; -ac];
    a1= [a1(1,:) ; ac];

    m=(a1(1,:)+a2(1,:)).';
end
    


end

