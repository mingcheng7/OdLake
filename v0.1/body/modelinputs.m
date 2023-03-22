% === One-day lake simulation model === % 
% version 0.1 
% By Ming Cheng 
% Last modification on 30/11/2022
% Module for input data 

function [tt,In_Z,In_A,In_T,In_DO,In_PAR,In_PO4,In_NO3,In_NH3,In_CyB,In_APSB,In_H2S,In_SO4,In_DFe,In_Fe3,In_MFe,In_DMn,In_PMn,...
    In_DFe_sed,In_DMn_sed,In_NO3_sed,In_SO4_sed,In_H2S_sed,In_PO4_sed,Wt,Inflw,Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names] ...
    = modelinputs(Datetime,init_filename,init_sheet,input_filename,input_sheet,param_filename,param_sheet,dt,n)

% Inputs:
%       Datetime        : Model running date [year, month, day]
%       Input filenames and sheetnames
%		dt		        : Time step 
%       n               : Number of iterated days 
% Outputs:
%		tt		        : Solution time domain (day)
%       In_Z            : Depths read from initial profiles file (m)
%       In_A            : Areas read from initial profiles file (m2)
%       In_T            : Initial temperature profile read from initial profiles file (deg C)
%       In_DO           : Initial dissolved oxygen profile read from initial profiles file (mg m-3) 
%       In_PAR          : Initial PAR profile read from initial profiles file (-) 
%       In_PO4          : Initial PO4 profile read from initial profiles file (mol m-3)
%       In_NO3          : Initial dissolved NO3 profile read from initial profiles file (mol m-3) 
%       In_NH3          : Initial dissolved ammonia (NH3) profile read from initial profiles file (mol m-3)
%       In_CyB          : Initial Cyanobacteria profile read from initial profiles file (cells/m3) 
%       In_APSB         : Initial aerobic phototrophic sulfur bacteria (APSB) profile read from initial profiles file (cells m-3) 
%       In_H2S          : Initial H2S profile read from initial profiles file (mol m-3) 
%       In_SO4          : Initial SO4 profile read from initial profiles file (mol m-3) 
%       In_DFe          : Initial dissolved Fe (Fe2+) profile read from initial profiles file (mol m-3) 
%       In_Fe3          : Initial particulate Fe (Fe2O3/Fe(OH)3/FeOOH) profile read from initial profiles file (mol m-3) 
%       In_MFe          : Initial pyrites (FeS2) profile read from initial profiles file (mol m-3) 
%       In_DMn          : Initial dissolved Mn (Mn2+) profile read from initial profiles file (mol m-3) 
%       In_PMn          : Initial particulate Mn (MnO2) profile read from initial profiles file (mol m-3) 
%       In_DFe_sed      : Initial dissolved Fe2+ concentration in porewater (mol m-3) 
%       In_DMn_sed      : Initial dissolved Mn2+ concentration in porewater (mol m-3) 
%       In_NO3_sed      : Initial NO3 concentration in porewater (mol m-3) 
%       In_SO4_sed      : Initial SO4 concentration in porewater (mol m-3) 
%       In_H2S_sed      : Initial H2S concentration in porewater (mol m-3) 
%       In_PO4_sed      : Initial PO4 concentration in porewater (mol m-3)
%		Wt		        : Weather data
%       Inflow          : Inflow data
%       Phys_par        : Main 16 parameters that are more or less fixed
%       Phys_par_range  : Minimum and maximum values for Phys_par (16 * 2)
%       Phys_par_names  : Names for Phys_par
%       Bio_par         : Main 13 parameters that are more or less site specific
%       Bio_par_range   : Minimum and maximum values for Bio_par (13 * 2)
%       Bio_par_names   : Names for Bio_par


global ies80;

% == Read model parameter file
[ParaMx,StrMx]=xlsread(param_filename,param_sheet);

% Main physical parameters 
Phys_par_names=StrMx(3:18,1);
Phys_par=ParaMx(3:18,2);
Phys_par_range=ParaMx(3:18,3:4);

% Main biological parameters 
Bio_par_names=StrMx(19:31,1);
Bio_par=ParaMx(19:31,2);
Bio_par_range=ParaMx(19:31,3:4);

% == Read morphometric and initial profile file

[InitMx,StrMx]=xlsread(init_filename,init_sheet);
In_Z=InitMx(3:end,1);
In_A=InitMx(3:end,2);
In_T=InitMx(3:end,3);
In_DO=InitMx(3:end,4);
In_PAR=InitMx(3:end,5);
In_PO4=InitMx(3:end,6);
In_NO3=InitMx(3:end,7);
In_NH3=InitMx(3:end,8);
In_CyB=InitMx(3:end,9);
In_APSB=InitMx(3:end,10);
In_H2S=InitMx(3:end,11);
In_SO4=InitMx(3:end,12);
In_DFe=InitMx(3:end,13);
In_Fe3=InitMx(3:end,14);
In_MFe=InitMx(3:end,15);
In_DMn=InitMx(3:end,16);
In_PMn=InitMx(3:end,17);
In_DFe_sed=InitMx(3:end,18);
In_DMn_sed=InitMx(3:end,19);
In_NO3_sed=InitMx(3:end,20);
In_SO4_sed=InitMx(3:end,21);
In_H2S_sed=InitMx(3:end,22);
In_PO4_sed=InitMx(3:end,23);
  
t = (datenum(Datetime):dt:(datenum(Datetime)+1))';		% Solution time domain 
tt = [];
for i=1:n
    tt = [tt;t(1:end-1)];
end 

% == Read input forcing data file

[InputMx,StrMx]=xlsread(input_filename,input_sheet);

In_Date=InputMx(3,1:3);
In_Met=InputMx(3,4:10);
In_Inflow=InputMx(3,11:20);


clear Wt
Wt=In_Met;
% Wt(:,1)  Global radiation (MJ (m^2 day)-1)
% Wt(:,2)  Cloud cover (-)
% Wt(:,3)  Air temperature (deg. C, at 2 m height)
% Wt(:,4)  Relative humidity (%, at 2 m height)
% Wt(:,5)  Air pressure (mbar)
% Wt(:,6)  Wind speed (m s-1 at 10 m height)
% Wt(:,7)  Precipitation (mm day-1)

clear Inflw
Inflw=In_Inflow;
% Inflw(:,1) Inflow volume (m3 day-1)
% Inflw(:,2) Inflow temperature (deg C)
% Inflw(:,3) Inflow dissolved oxygen (mol m-3) 
% Inflw(:,4) Inflow PO4 concentration (mol m-3)
% Inflw(:,5) Inflow dissolved NO3 anion concentration (mol m-3) 
% Inflw(:,6) Inflow SO4 concentration (mol m-3) 
% Inflw(:,7) Inflow dissolved Fe concentration (mol m-3) 
% Inflw(:,8) Inflow particulate Fe(III) concentration (mol m-3) 
% Inflw(:,9) Inflow dissolved Mn concentration (mol m-3) 
% Inflw(:,10) Inflow particulate Mn concentration (mol m-3) 


% International Equation of State 1980
% 5-order polynomial for density as function of temperature
ies80 = [6.536332e-9,-1.120083e-6,1.001685e-4,-9.09529e-3,6.793952e-2,999.842594];


% Default turbulence and wind shelter parameterization (Hondzo and Stefan, 1993; Ellis et al., 1991)
if(isnan(Phys_par(2)))
  Phys_par(2) = 0.00706*(In_A(1)/1e6)^0.56; % default diffusion coeff. parameterisation 0.000817*(In_A(1)/1e6)^0.56
end

if(isnan(Phys_par(3)))
  Phys_par(3) = 8.98e-4;		%default value for diffusion coeff. in ice-covered water   
end

if(isnan(Phys_par(4)))
  Phys_par(4) = 7e-5;			% default minimum allowed stability frequency, N2 > N0 <=> Kz < Kmax (1/s2)    		
end

if(isnan(Phys_par(5)))
  Phys_par(5) =  1-exp(-0.3*In_A(1)/1e6);			% default wind sheltering parameterisation		
end
