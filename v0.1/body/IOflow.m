% === One-day lake simulation model === % 
% version 0.1 
% By Ming Cheng 
% Last modification on 30/11/2022
% River inflow module 


% INPUTS:
%   vertical arrays:    z is model grid, Vz is layer volume, Tz layer property (e.g. temperature);
%   scalars:            lvlD is the grid depth level *above* which the inflow settles (m), 
%                       Iflw is the inflow volume (m3/day), and 
%                       T_Iflw is the property of the inflow, and dz the grid stepsize (layer thickness)
% OUTPUTS: Cz is the new property profile after inflow

function [Cz]=IOflow(dz, z, Vz, Tz, lvlD, Iflw, T_Iflw)

dum=find(z>=lvlD); %same as lvlD but in grid level numbers
lvlG=dum(1);
Cz=Tz;

if (lvlD==0) %if inflow is lighter than surface water, mix it with first layer   
Cz(1)=(Vz(1)*Tz(1)+Iflw*T_Iflw)/(Vz(1)+Iflw);    
     
else    %otherwise add the inflow to the appropriate 
        %depth level and "lift" the water column above
        
 Trev_pour=flipud([Tz(1:lvlG-1); T_Iflw; 0]); %up-side down Tz to be poured in (plus an extra zero)       
 Vzrev_pour=flipud([Vz(1:lvlG-1); Iflw; 0]);  %up-side down Vz to be poured in (plus an extra zero)     
 Vzrev_fill=flipud([Vz(1:lvlG-1)]);  %up-side down Vz to be filled in     
  
 for i=1:length(Vzrev_fill)
  inxA=find(cumsum(Vzrev_pour)>Vzrev_fill(i)); 
  inxB=inxA(1); %index of the last layer to be poured in (partly)
  ShakerV=[Vzrev_pour(1:inxB-1); Vzrev_fill(i)-sum(Vzrev_pour(1:inxB-1))];
  ShakerT=[Trev_pour(1:inxB)];
  Cz(lvlG-i)=sum(ShakerV.*ShakerT)/sum(ShakerV); %new property after mixing
  Vzrev_pour(1:inxB)=Vzrev_pour(1:inxB)-ShakerV; %subtract poured volumes from the "reserves"
 end
            
end %if

