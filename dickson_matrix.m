function [ A_caps, A_sw1, A_sw2 ] = dickson_matrix(n_stages,in_cap)
%% DICKSON_MATRIX Creates the incindence matrixes for an N_STAGES (#capacitors)
% Dikson Ladder converter
% 
% Retunr values 
%   A_caps -> Capacitor Incidnece matrix Inculuding the voltage supply
%             connected to the first node
%   A_sw1  -> Phase 1 ON Switches 
%   A_sw2  -> Phase 2 ON switches
%
%The firs row of the incidence matrix corresponds to the input supply node
%The Phase 1 is the one that drawns power from the load
%
%  Philips Research, Eindhoven,  Netherlands
%  julia.delos@philps.com
%
    
%% Initialize matrix
A_caps = zeros(n_stages+2,n_stages);

if n_stages>2  
    n_sw1 = 4+floor((n_stages-3)./2);
else
    n_sw1 = 2;
end

n_sw2 =n_sw1- rem(n_stages,2);  
    
A_sw1  = zeros(n_stages+3,n_sw1);
A_sw2  = zeros(n_stages+3,n_sw2);

j=0;
for i=1:n_stages
   %Generate Capacitor Incidence matrix
   if (n_stages - i) < 3         
    A_caps(i,i)=1; %Caps matrix
    if i < n_stages
    A_caps(end-j,i)=-1;
    j=1+j;
    end
   else
    A_caps(i,i)=1;
    A_caps(i+2,i)=-1;
   end
end

%Add the source line
A_caps = [zeros(1,n_stages); A_caps];

for i=1:n_sw2
   %Phase 1 matrix
   A_sw1((2*i-1),i) =  1;
   if i < n_sw2
     A_sw1((2*i),i) = -1;
   end
end

    %Phase 2
    A_sw2(2:end,:)=A_sw1(1:end-1,1:n_sw2);
    
    if n_stages > 2
        if  rem(n_stages,2) 
            %Odd
            A_sw1(:,end-1:end)=zeros(n_stages+3,2);

            A_sw1((end-2):end,(end-1):end) = [0 1; 1 0; 0 -1];
        else
            A_sw2(:,end-1:end)=zeros(n_stages+3,2);
            A_sw2((end-2):end,(end-1):end) = [0 1; 1 0; 0 -1];
        end    
    else
        %Two stage converter
        A_sw1(end,:)=[];
        A_sw1(end,end)=-1;
        
        A_sw2(end,:)=[];  
        
        A_caps(end-1,:)=[];
    end

if in_cap  
A_caps = [zeros(size(A_caps,1),1) A_caps];
A_caps(1)=1;
end
end

