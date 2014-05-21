function [ Bf ] = inM2loopM( A )
%INM2LOOPM Converts incidence matrix A to loop matrix Bf with only the 
% fundamental loops
%  Philips Research, Eindhoven,  Netherlands
%  julia.delos@philps.com
%
    


Aef=rref(A);
Aef(all(Aef==0,2),:)=[];
[n m] = size(Aef);
ln=[];
U=Aef(:,1:n);

isArranged=0;

if ~(isequal(U,eye(n)))     
    for i=1:n
        for j=1:m
        if Aef(i,j)==1 
           tw(i) = j;
           break
        end
       end
    end
    
    for i=1:m
        if ~max(tw==i)
            ln=[ln i];
        end
    end
    Vtrf=[tw ln];
    
    for i=1:m
        Vbtrf(i)=find(Vtrf==i);
    end
    Aef=Aef(:,Vtrf);
    isArranged=1;
end

Bt = -Aef(:,n+1:end)';
[n m] =size(Bt);
Bf = [Bt eye(n)];

if isArranged
    Bf=Bf(:,Vbtrf);
end

