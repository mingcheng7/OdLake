% === One-day lake simulation model === % 
% version 0.6 
% By Ming Cheng 
% Main module 

function[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt,APSBzt,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,varargin)

% Inputs (to solvemodel function) 
% Datetime: Model run date [year, month, day]
% + input filenames and sheetnames 

% Inputs (from the input module)
%		tt		        : Solution time domain (day)
%       In_Z            : Depths read from initial profiles file (m)
%       In_A            : Areas read from initial profiles file (m2)
%       In_T            : Initial temperature profile read from initial profiles file (deg C)
%       In_DO           : Initial dissolved oxygen profile read from initial profiles file (mg m-3) 
%       In_pH           : Initial pH profile read from initial profiles file (-) 
%       In_PAR          : Initial PAR profile read from initial profiles file (µmol.s-1.m-2) 
%       In_PO4          : Initial PO4 profile read from initial profiles file (mg m-3)
%       In_NO3          : Initial dissolved NO3- profile read from initial profiles file (mg m-3) 
%       In_NH3          : Initial dissolved ammonia (NH3) profile read from initial profiles file (mg m-3)
%       In_CyB          : Initial Cyanobacteria profile read from initial profiles file (cells/m3) 
%       In_APSB         : Initial aerobic phototrophic sulfur bacteria (APSB) profile read from initial profiles file (cells/m3) 
%       In_H2S          : Initial H2S profile read from initial profiles file (mg/m3) 
%       In_SO4          : Initial SO4 profile read from initial profiles file (mg/m3) 
%       In_DFe          : Initial dissolved Fe (DFe) profile read from initial profiles file (mg/m3) 
%       In_Fe3          : Initial Fe(III) oxyhydroxides profile read from initial profiles file (mg/m3) 
%       In_MFe          : Initial pyrites (FeS2) profile read from initial profiles file (mg/m3) 
%       In_DMn          : Initial dissolved Mn (Mn2+) profile read from initial profiles file (mol/m3) 
%       In_PMn          : Initial particulate Mn (MnO2) profile read from initial profiles file (mol/m3) 
%       In_DFe_sed      : Initial dissolved Fe2+ concentration in porewater (mol/m3) 
%       In_DMn_sed      : Initial dissolved Mn2+ concentration in porewater (mol/m3) 
%       In_NOx_sed      : Initial NOx concentration in porewater (mol/m3) 
%       In_SO4_sed      : Initial SO4 concentration in porewater (mol/m3) 
%       In_H2S_sed      : Initial H2S concentration in porewater (mol/m3) 
%       In_PO4_sed      : Initial PO4 concentration in porewater (mol/m3)
%       Ice0            : Initial conditions, ice and snow thicknesses (m) (Ice, Snow)
%		Wt		        : Weather data
%       Inflow          : Inflow data
%       Phys_par        : Main 23 parameters that are more or less fixed
%       Phys_par_range  : Minimum and maximum values for Phys_par (23 * 2)
%       Phys_par_names  : Names for Phys_par
%       Bio_par         : Main 15 parameters that are more or less site specific
%       Bio_par_range   : Minimum and maximum values for Bio_par (15 * 2)
%       Bio_par_names   : Names for Bio_par

% Outputs
%		Qst         : Estimated surface heat fluxes ([sw, lw, sl] * tt) (W m-2)
%		Kzt	        : Predicted vertical diffusion coefficient (tt * zz) (m2 d-1)
%		Tzt	        : Predicted temperature profile (tt * zz) (deg C)
%       rhozt       : Predicted density profile (tt * zz) (kg m-3)
%       DOzt        : Predicted dissolved oxygen profile (tt * zz) (mol m-3)
%       pHzt        : Predicted pH profile (tt * zz) (-) 
%		PO4zt	    : Predicted PO4 profile (tt * zz) (mol m-3)
%       NO3zt       : Predicted dissolved NO3 profile (tt * zz) (mol m-3) 
%       NH3zt       : Predicted dissolved ammonia (NH3) and ammonium (NH4+) profile (tt * zz) (mol m-3) 
%       CyBzt       : Predicted Cyanobacteria profile (tt * zz) (cells/m3) 
%       APSBzt      : Predicted APSB profile (tt * zz) (cells/m3)
%       H2Szt       : Predicted H2S profile (tt * zz) (mol/m3) 
%       SO4zt       : Predicted SO4 prifile (tt * zz) (mol/m3) 
%       DFezt       : Predicted dissolved Fe profile (tt * zz) (mol/m3) 
%       PFezt       : Predicted particulate Fe profile (tt * zz) (mol/m3) 
%       delta_DFezt : Predicted delta 56 Fe profile in dissolved Fe (tt * zz) (‰) 
%       delta_PFezt : Predicted delta 56 Fe profile in particulate Fe (tt * zz) (‰) 
%       DMnzt       : Predicted dissolved Mn profile (tt * zz) (mol/m3) 
%       PMnzt       : Predicted particulate Mn profile (tt * zz) (mol/m3) 
%		Qz_sed      : Predicted sediment-water heat flux (tt * zz) (W m-2, normalised to lake surface area)
%       PARzt       : Predicted PAR (tt * zz) (mol.s-1.m-2)
%       G_cybzt     : Predicted growth rate profile for cyanobacteria (tt * zz) (day-1)
%       G_apsbzt    : Predicted growth rate profile for APSB (tt * zz) (day-1)

% These variables are still global and not transferred by functions
global ies80;

tic  % Start stopwatch timer 
disp(['Running MyLake on ' datestr(datenum(Datetime))  ' ...']);  % Display model running date 

% ===Switches===
river_inflow_switch=1;          %river inflow: 0=no, 1=yes
selfshading_switch=1;           %light attenuation by chlorophyll a: 0=no, 1=yes
nfixation_switch=0;             %nitrogen fixation: 0=no, 1=yes
% ==============

dt=0.1; % setting the time grid 

% Input data from input module 
if (nargin>8) %if optional command line parameter input is used 
  disp('Bypassing input files...Running with input data & parameters given on command line');
 [tt,In_Z,In_A,In_T,In_DO,In_pH,In_PAR,In_PO4,In_NO3,In_NH3,In_CyB,In_APSB,In_H2S,In_SO4,In_DFe,In_Fe3,In_MFe,In_DMn,In_PMn,...
    In_DFe_sed,In_DMn_sed,In_NO3_sed,In_SO4_sed,In_H2S_sed,In_PO4_sed,Wt,Inflw,Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names] ...
    = deal(varargin{:});
else
% Read input data
  [tt,In_Z,In_A,In_T,In_DO,In_pH,In_PAR,In_PO4,In_NO3,In_NH3,In_CyB,In_APSB,In_H2S,In_SO4,In_DFe,In_Fe3,In_MFe,In_DMn,In_PMn,...
    In_DFe_sed,In_DMn_sed,In_NO3_sed,In_SO4_sed,In_H2S_sed,In_PO4_sed,Wt,Inflw,Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names] ...
    = modelinputs(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,dt);
end

%load H:\MyLake\air_sea\albedot1.mat; %load albedot1 table, in order to save execution time
load albedot1.mat; %load albedot1 table, in order to save execution time


d2s=60*60*24; % day to seconds parameter 

% Load the parameters 
% Unpack the more fixed parameter values from input array "Phys_par"
dz = Phys_par(1); %grid stepsize (m)

   zm = In_Z(end); %max depth
   zz = [0:dz:zm-dz]'; %solution depth domain

Kz_K1 = Phys_par(2);          % open water diffusion parameter (-)
Kz_N0 = Phys_par(3);          % min. stability frequency (s-2)
C_shelter = Phys_par(4);      % wind shelter parameter (-)
lat = Phys_par(5);            % latitude (decimal degrees)
lon = Phys_par(6);            % longitude (decimal degrees)
PAR_C_sat = Phys_par(7);      % PAR saturation level for cyanobacteria growth (mol(quanta) m-2 s-1) 
PAR_S_sat = Phys_par(8);      % PAR saturation level for sulfur bacteria growth (mol(quanta) m-2 s-1) 
f_par = Phys_par(9);          % Fraction of PAR in incoming solar radiation (-)
beta_chl = Phys_par(10);      % Optical cross_section of chlorophyll (m2/mg)
w_ad = Phys_par(11);          % Advection velocity (m day-1)
w_cyb = Phys_par(12);         % settling velocity for cyanobacteria (m day-1)
w_apsb = Phys_par(13);        % settling velocity for APSB (m day-1)
w_Fe3 = Phys_par(14);         % settling velocity for FeOOH (m day-1)
w_PMn = Phys_par(15);         % settling velocity for MnO2 (m day-1)
w_MFe = Phys_par(16);         % settling velocity for FeS2 (m day-1)

% Unpack the more site specific parameter values from input array "Bio_par"

swa_b0 = Bio_par(1);          % non-PAR light atteneuation coeff. (m-1)
swa_b1 = Bio_par(2);          % PAR light atteneuation coeff. (m-1)
m_cyb = Bio_par(3);           % loss rate (1/day) of Cyanobacteria at 20 deg C
m_apsb = Bio_par(4);          % loss rate (1/day) of APSB at 20 deg C
g_cyb = Bio_par(5);           % specific growth rate (1/day) of Cyanobacteria at 20 deg C
g_apsb = Bio_par(6);          % specific growth rate (1/day) of APSB at 20 deg C
r_OM_cyb = Bio_par(7);        % OM parameter of Cyanobacteria (mol OM/cell) 
r_OM_apsb = Bio_par(8);       % OM parameter of APSB (mol OM/cell) 
PO4_half = Bio_par(9);       % Half saturation growth PO4 level (mol/m3)
NO3_half = Bio_par(10);       % Half saturation growth NHx level (mol/m3)
PAR_C_half = Bio_par(11);     % Half saturation growth PAR level for Cyanobacteria (mol m-2 s-1)
PAR_S_half = Bio_par(12);     % Half saturation growth PAR level for sulfur bacteria (mol m-2 s-1)
H2S_half = Bio_par(13);       % Half saturation growth H2S level (mol/m3)
Nfix = Bio_par(14);           % Nitrogen fixation rate of cyanobacteria (mol N per cell per day) at 4 degree C

% Other parameters not appear in the input module 

Nz=length(zz); % total number of layers of the water column 

theta_m=exp(0.1*log(2)); % phytoplankton growth rate base 
e_par=240800; % average energy of PAR photons (J mol-1)

% diffusion parameterisation exponents 
Kz_exp=0.2; %0.43 %0.2

% ice & snow parameter values
Tf=0;               % water freezing point temperature (deg C)
ksw=1e-3;           % Porewater_water mass transfer coefficient (m/d)

% Initialise the output data matrices 
Qst = zeros(3,length(tt)); 
Kzt = zeros(Nz,length(tt));  
Tzt = zeros(Nz,length(tt)); 
rhozt = zeros(Nz,length(tt));
DOzt = zeros(Nz,length(tt));
pHzt = zeros(Nz,length(tt));
PO4zt = zeros(Nz,length(tt)); 
NO3zt = zeros(Nz,length(tt));
NH3zt = zeros(Nz,length(tt));
CyBzt = zeros(Nz,length(tt));
APSBzt = zeros(Nz,length(tt)); 
H2Szt = zeros(Nz,length(tt));
SO4zt = zeros(Nz,length(tt));
DFezt = zeros(Nz,length(tt));
PFezt = zeros(Nz,length(tt));
DMnzt = zeros(Nz,length(tt));
PMnzt = zeros(Nz,length(tt));
Qzt_sed = zeros(Nz,length(tt));
PARzt = zeros(Nz,length(tt));  
G_cybzt = zeros(Nz,length(tt)); 
G_apsbzt = zeros(Nz,length(tt));

% Initial profiles 
Az = interp1(In_Z,In_A,zz);                      % Initial crossing area in water layers (m3) 
Vz = dz * (Az + [Az(2:end); 0])/2;               % Initial volume in water layers (m3)
Tz = interp1(In_Z,In_T,zz+dz/2);                 % Initial temperature distribution (deg C)
DOz = interp1(In_Z,In_DO,zz+dz/2);               % Initial dissolved oxygen distribution (mol m-3) 
pHz = interp1(In_Z,In_pH,zz+dz/2);               % Initial pH distribution (-) 
PARz = 10.^(interp1(In_Z,log10(In_PAR),zz+dz/2));% Initial PAR distribution (µmol.s-1.m-2)
PO4z = interp1(In_Z,In_PO4,zz+dz/2);	         % Initial PO4 distribution (mol m-3)
NO3z = interp1(In_Z,In_NO3,zz+dz/2);             % Initial dissolved NOx distribution (mol m-3) 
NH3z = interp1(In_Z,In_NH3,zz+dz/2);             % Initial dissolved ammonia distribution (mol m-3)
CyBz = interp1(In_Z,In_CyB,zz+dz/2);             % Initial Cyanobacteria distribution (cells m-3) 
APSBz = interp1(In_Z,In_APSB,zz+dz/2);           % Initial aerobic phototrophic sulfur bacteria distribution (cells m-3) 
H2Sz = interp1(In_Z,In_H2S,zz+dz/2);             % Initial H2S distribution (mol m-3) 
SO4z = interp1(In_Z,In_SO4,zz+dz/2);             % Initial SO4 distribution (mol m-3) 
DFez = interp1(In_Z,In_DFe,zz+dz/2);             % Initial dissolved Fe distribution (mol m-3) 
Fe3z = interp1(In_Z,In_Fe3,zz+dz/2);             % Initial particulate Fe(III) distribution (mol m-3) 
MFez = interp1(In_Z,In_MFe,zz+dz/2);             % Initial pyrite (FeS2) distribution (mol m-3) 
DMnz = interp1(In_Z,In_DMn,zz+dz/2);             % Initial dissolved Mn distribution (mol m-3) 
PMnz = interp1(In_Z,In_PMn,zz+dz/2);             % Initial particulate Mn distribution (mol m-3) 
DFe_sed = interp1(In_Z,In_DFe_sed,zz+dz/2);      % Initial dissolved Fe2+ distribution in porewater (mol m-3) 
DMn_sed = interp1(In_Z,In_DMn_sed,zz+dz/2);      % Initial dissolved Mn2+ distribution in porewater (mol m-3) 
NO3_sed = interp1(In_Z,In_NO3_sed,zz+dz/2);      % Initial NOx distribution in porewater (mol m-3) 
SO4_sed = interp1(In_Z,In_SO4_sed,zz+dz/2);      % Initial SO4 distribution in porewater (mol m-3) 
H2S_sed = interp1(In_Z,In_H2S_sed,zz+dz/2);      % Initial H2S distribution in porewater (mol m-3) 
PO4_sed = interp1(In_Z,In_PO4_sed,zz+dz/2);      % Initial PO4 distribution in porewater (mol m-3) 

clear Tzs_sed
for j=1:Nz
    Tzs_sed(:,j) = interp1([0.2 10], [Tz(j) 4], [0.2:0.2:2 2.5:0.5:10])';
end

% >>>>>> Start the time loop >>>>>>
for i=1:length(tt)
    % Calculate surface heat flux (W m-2), wind stress (N m-2), ...
    [Qsw,Qlw,Qsl,tau,DayFrac,DayFracHeating] = heatflux(tt(i),Wt(1),Wt(2),Wt(3),Wt(4),Wt(5),Wt(6),Tz(1),lat,albedot1);     %Qlw and Qsl are functions of Tz(1) 

    if (i==1)
    % Calculate the estimated Chlorophyll concentration 
    Chlz = (1.1166e-6*CyBz+2.85e-8*APSBz)./(316.7+9.6*Tz+133*PARz.*exp(-0.2*Tz));
    
    % Calculate total mean PAR and non-PAR light attenuation coefficient in water 
    epsilonz_ave = zeros(Nz,1);
    epsilonz_NP_ave = zeros(Nz,1);

    if (selfshading_switch==1)
        epsilonz = swa_b1*ones(Nz,1) + beta_chl*Chlz; 
        epsilonz_NP = swa_b0*ones(Nz,1) + beta_chl*Chlz; 
        for j=1:Nz
            epsilonz_ave(j) = mean(swa_b1*ones(j,1) + beta_chl*Chlz(1:j)); 
            epsilonz_NP_ave(j) = mean(swa_b0*ones(j,1) + beta_chl*Chlz(1:j)); 
        end
    else
        epsilonz=swa_b1*ones(Nz,1);         % PAR light attenuation coeff.
        epsilonz_NP=swa_b0*ones(Nz,1);      % non-PAR light attenuation coeff.
        epsilonz_ave=swa_b1*ones(Nz,1);     % PAR light attenuation coeff. in average 
        epsilonz_NP_ave=swa_b0*ones(Nz,1);  % non-PAR light attenuation coeff. in average 
    end

    % solve the light irradiance profile 
    PARz=(3/2)*f_par*Qsw/(DayFrac*e_par)*exp(-epsilonz_ave.*zz); % light irridiance profile at noon
    end

    % Cyanobacteria 
    H_sw_C_z=NaN*zeros(Nz,1); 

    U_sw_z=PARz/PAR_C_sat;    % relative PAR 
    idx_u=find(U_sw_z<=1);  % undersaturated 
    idx_s=find(U_sw_z>1);   % saturated 

    H_sw_C_z(idx_u)=(2/3)*U_sw_z(idx_u);  % undersaturated
 
    dum_a=sqrt(U_sw_z);
    dum_b=sqrt(U_sw_z-1);
    H_sw_C_z(idx_s)=(2/3)*U_sw_z(idx_s) + log((dum_a(idx_s) + dum_b(idx_s))./(dum_a(idx_s) ...  % saturated
        - dum_b(idx_s))) - (2/3)*(U_sw_z(idx_s)+2).*(dum_b(idx_s)./dum_a(idx_s));

    % Sulfur bacteria 
    H_sw_S_z=NaN*zeros(Nz,1); 

    U_sw_z=PARz/PAR_S_sat;    % relative PAR 
    idx_u=find(U_sw_z<=1);  % undersaturated 
    idx_s=find(U_sw_z>1);   % saturated 

    H_sw_S_z(idx_u)=(2/3)*U_sw_z(idx_u);  % undersaturated
 
    dum_a=sqrt(U_sw_z);
    dum_b=sqrt(U_sw_z-1);
    H_sw_S_z(idx_s)=(2/3)*U_sw_z(idx_s) + log((dum_a(idx_s) + dum_b(idx_s))./(dum_a(idx_s) ...  % saturated
        - dum_b(idx_s))) - (2/3)*(U_sw_z(idx_s)+2).*(dum_b(idx_s)./dum_a(idx_s));

    % Vertical turbulent diffusion 
    g=9.81;  % Gravity acceleration (m s-2) 
    rho=polyval(ies80,max(0,Tz(:)))+min(Tz(:),0); % Water density (kg m-3) 
    N2=g*(diff(log(rho))./diff(zz)); % Brunt_Vaisala frequency (s-2) 
    
    Kz=Kz_K1*max(Kz_N0,N2).^(-Kz_exp);  % Vertical diffusion coeff. for ice free (m2 day-1)
   
    % Update biological rates for cyanobacteria 
    limpp_cyb = (min([(PO4z./(PO4_half+PO4z))';(NO3z./(NO3_half+NO3z))';(PARz./(PAR_C_half+PARz))']))'; % limiting term for phytoplankton production
    Growth_bioz=g_cyb*dt*theta_m.^(Tz-20).*limpp_cyb.*(DayFrac./(dz*epsilonz)).*diff([-H_sw_C_z;-H_sw_C_z(end)]).*min(DOz/0.1,ones(Nz,1)); 
    R_bio_cybz=(min([Growth_bioz';(PO4z./(CyBz*r_OM_cyb))';(NO3z./(CyBz*r_OM_cyb*16))']))'; 
        % growth rate is limited by available phosphorus, nitrogen 

    dCyB = CyBz.*R_bio_cybz; 
    dNO3 = dCyB*r_OM_cyb*16; 
    dPO4 = dCyB*r_OM_cyb; 
    dDO = dCyB*r_OM_cyb*150; 

    % Solve the profile after regeneration 
    CyBz = CyBz+dCyB;
    NO3z = NO3z-dNO3;
    PO4z = PO4z-dPO4;
    
    DOsatz=SolDO(Tz);
    DOz = DOz + dDO; 
    DOratioz = DOz./DOsatz; % Dissolved oxygen dissolved ratio (-) 
    idx = find(DOratioz>1.05); % Oversaturated DO ratio cannot be larger than 105% 
    DOratioz(idx) = 1.05;
    DOz(idx) = DOratioz(idx).*DOsatz(idx);

    % Nitrogen fixation of cyanobacteria 
    if (nfixation_switch==1)
        R_Nfix = Nfix*dt*theta_m.^(Tz-4);  % Nitrogen fixation rate in each layer 
        dNH3 = CyBz.*R_Nfix;
        NH3z = NH3z + dNH3;
    end

    % Update biological rates for APSB 
    limpp_apsb = (min([(PO4z./(PO4_half+PO4z))';(NO3z./(NO3_half+NO3z))';(H2Sz./(H2S_half+H2Sz))';(PARz./(PAR_S_half+PARz))']))'; % limiting term for phytoplankton production
    Growth_bioz=g_apsb*dt*theta_m.^(Tz-20).*limpp_apsb.*(DayFrac./(dz*epsilonz)).*diff([-H_sw_S_z;-H_sw_S_z(end)]); 
    R_bio_apsbz=(min([Growth_bioz';(PO4z./(APSBz*r_OM_apsb))';(NO3z./(APSBz*r_OM_apsb*16))';(H2Sz./(APSBz*r_OM_apsb*75))']))'; 
        % growth rate is limited by available phosphorus, nitrogen 
    
    dAPSB = APSBz.*R_bio_apsbz;  
    dNO3 = dAPSB*r_OM_apsb*16; 
    dPO4 = dAPSB*r_OM_apsb;
    dH2S = dAPSB*r_OM_apsb*75;
    dS = dH2S; 

    % Solve the profile after regeneration 
    APSBz = APSBz+dAPSB;
    NO3z = NO3z-dNO3;
    PO4z = PO4z-dPO4;
    H2Sz = H2Sz-dH2S;
    SO4z = SO4z + dS*0.25; 
    H2Sz = H2Sz + dS*0.75; 

    % Remineralisation 
    [CyBz,DOz,NH3z,NO3z,PO4z,SO4z,H2Sz,DFez,Fe3z,DMnz,PMnz] = remineralisation(dt,zz,Tz,CyBz,DOz,NH3z,NO3z,PO4z,SO4z,H2Sz,DFez,Fe3z,DMnz,PMnz,m_cyb,r_OM_cyb,theta_m);
    [APSBz,DOz,NH3z,NO3z,PO4z,SO4z,H2Sz,DFez,Fe3z,DMnz,PMnz] = remineralisation(dt,zz,Tz,APSBz,DOz,NH3z,NO3z,PO4z,SO4z,H2Sz,DFez,Fe3z,DMnz,PMnz,m_apsb,r_OM_apsb,theta_m);

    % Chemical reactions 
    [DOz,H2Sz,SO4z,DFez,Fe3z,MFez,DMnz,PMnz] = Chemfun(dt,zz,DOz,H2Sz,SO4z,DFez,Fe3z,MFez,DMnz,PMnz);

    Fi_ad=tridiag_HAD_v11([NaN;Kz],w_ad,Vz,Az,dz,dt);  % Tridiagonal matrix for advection and diffusion 

    CyBz = Fi_ad \ CyBz;
    APSBz = Fi_ad \ APSBz;
    DOz = Fi_ad \ DOz; % Solving for new DO profile (diffusion) 
    PO4z = Fi_ad \ PO4z; 
    NO3z = Fi_ad \ NO3z; 
    NH3z = Fi_ad \ NH3z; 
    H2Sz = Fi_ad \ H2Sz; 
    SO4z = Fi_ad \ SO4z; 
    DFez = Fi_ad \ DFez;
    Fe3z = Fi_ad \ Fe3z; % Different sinking velocity 'cause the different density of particles 
    MFez = Fi_ad \ MFez;
    DMnz = Fi_ad \ DMnz;
    PMnz = Fi_ad \ PMnz; % Different sinking velocity 'cause the different density of particles 

    % Porewater exchange 
    % Porewater to water 
    PwwFrac=ksw*dt*(-diff([Az;0]))./Vz; % Fraction between resuspended porewater and water layer volumes 
    DFez = (1-PwwFrac).*DFez + PwwFrac.*DFe_sed; 
    DMnz = (1-PwwFrac).*DMnz + PwwFrac.*DMn_sed; 
    SO4z = (1-PwwFrac).*SO4z + PwwFrac.*SO4_sed; 
    NO3z = (1-PwwFrac).*NO3z + PwwFrac.*NO3_sed;
    H2Sz = (1-PwwFrac).*H2Sz + PwwFrac.*H2S_sed; 
    PO4z = (1-PwwFrac).*PO4z + PwwFrac.*PO4_sed;

    % Sediments exchange calculation 
    % Calculate the thickness ratio of newly settled net sedimentation and
    % mix these two to get new sediment 
    delCyB=NaN*ones(Nz,1); % Initialise 
    delAPSB=NaN*ones(Nz,1); % Initialise 
    delFe3=NaN*ones(Nz,1); % Initialise 
    delMFe=NaN*ones(Nz,1); % Initialise 
    delPMn=NaN*ones(Nz,1); % Initialise 

    delA=diff([Az;0]); % Area difference for i layer (OBS: negative) 
    meanA=0.5*(Az+[Az(2:end);0]);

    % Sedimentation is calculated from "Funnelling-NonFunnelling" difference 
    delCyB(1)=(0-CyBz(1)*delA(1)/meanA(1))/(dz/(dt*w_cyb)+1);
    delAPSB(1)=(0-APSBz(1)*delA(1)/meanA(1))/(dz/(dt*w_apsb)+1); 
    delFe3(1)=(0-Fe3z(1)*delA(1)/meanA(1))/(dz/(dt*w_Fe3)+1);
    delMFe(1)=(0-MFez(1)*delA(1)/meanA(1))/(dz/(dt*w_MFe)+1);
    delPMn(1)=(0-PMnz(1)*delA(1)/meanA(1))/(dz/(dt*w_PMn)+1);

    for ii=2:Nz
        delCyB(ii)=(delCyB(ii-1)-CyBz(ii)*delA(ii)/meanA(ii))/(dz/(dt*w_cyb)+1);
        delAPSB(ii)=(delAPSB(ii-1)-APSBz(ii)*delA(ii)/meanA(ii))/(dz/(dt*w_apsb)+1);
        delFe3(ii)=(delFe3(ii-1)-Fe3z(ii)*delA(ii)/meanA(ii))/(dz/(dt*w_Fe3)+1);
        delMFe(ii)=(delMFe(ii-1)-MFez(ii)*delA(ii)/meanA(ii))/(dz/(dt*w_MFe)+1);
        delPMn(ii)=(delPMn(ii-1)-PMnz(ii)*delA(ii)/meanA(ii))/(dz/(dt*w_PMn)+1);
    end

    CyBz=CyBz-delCyB;
    APSBz=APSBz-delAPSB;
    Fe3z=Fe3z-delFe3;
    MFez=MFez-delMFe;
    PMnz=PMnz-delPMn;

    % Inflow Calculation 

    % Inflw(:,1) Inflow volume (m3 day-1)
    % Inflw(:,2) Inflow temperature (deg C)
    % Inflw(:,3) Inflow dissolved oxygen (mg m-3) 
    % Inflw(:,4) Inflow PO4 concentration (mg m-3)
    % Inflw(:,5) Inflow dissolved NOx anion concentration (mg m-3) 
    % Inflw(:,6) Inflow Cyanobacteria concentration (cells m-3) 
    % Inflw(:,7) Inflow SO4 concentration (mg m-3) 
    % Inflw(:,8) Inflow dissolved Fe concentration (mg m-3) 
    % Inflw(:,9) Inflow particulate Fe(III) concentration (mg m-3) 
    % Inflw(:,10) Inflow dissolved Mn concentration (mol m-3) 
    % Inflw(:,11) Inflow particulate Mn concentration (mol m-3) 

    if (river_inflow_switch==1)
        Iflw = Inflw(1)*dt;   % inflow rate 
        Iflw_T = Inflw(2);    % inflow temperature
        if (Iflw_T<Tf)        % negative temperatures changed to Tf
             Iflw_T=Tf;
        end
        Iflw_DO = Inflw(3);   % inflow DO concentration 
        Iflw_PO4 = Inflw(4);  % inflow PO4 concentration
        Iflw_NOx = Inflw(5);  % inflow NOx concentration
        Iflw_SO4 = Inflw(6);  % inflow SO4 concentration 
        Iflw_DFe = Inflw(7);  % inflow dissolved Fe concentration 
        Iflw_Fe3 = Inflw(8);  % inflow particulate Fe concentration 
        Iflw_DMn = Inflw(9); % inflow dissolved Mn concentration 
        Iflw_PMn = Inflw(10); % inflow particulate Mn concentration 

        if (Iflw>0)
            if (isnan(Iflw_T))
                lvlD=0;
                Iflw_T=0;
            else
                rho=polyval(ies80,max(0,Tz(:)))+min(Tz(:),0);
                rho_Iflw=polyval(ies80,max(0,Iflw_T))+min(Iflw_T,0);
                lvlG=find(rho>=rho_Iflw);
                if (isempty(lvlG))
                    lvlG=length(rho);
                end
                lvlD=zz(lvlG(1));  % In which layer of water will be put by the inflow 
            end

            % Update the profiles after the inflow input 
            DOz = IOflow(dz,zz,Vz,DOz,lvlD,Iflw,Iflw_DO); % Dissolved oxygen 
            PO4z = IOflow(dz,zz,Vz,PO4z,lvlD,Iflw,Iflw_PO4);  % Dissolved PO4 
            NO3z = IOflow(dz,zz,Vz,NO3z,lvlD,Iflw,Iflw_NOx);  % Dissolved NOx
            SO4z = IOflow(dz,zz,Vz,SO4z,lvlD,Iflw,Iflw_SO4);  % Dissolved SO4 
            DFez = IOflow(dz,zz,Vz,DFez,lvlD,Iflw,Iflw_DFe);  % Dissolved Fe2+
            Fe3z = IOflow(dz,zz,Vz,Fe3z,lvlD,Iflw,Iflw_Fe3);  % Particulate Fe(III)
            DMnz = IOflow(dz,zz,Vz,DMnz,lvlD,Iflw,Iflw_DMn);  % Dissolved Mn2+ 
            PMnz = IOflow(dz,zz,Vz,PMnz,lvlD,Iflw,Iflw_PMn);  % Particulate Mn  

            % Dissolved oxygen dissolved ratio test 
            DOsatz = SolDO(Tz); % Get new DO solubility 
            DOratioz = DOz./DOsatz; 
            idx = find(DOratioz>1.05); % Oversaturated DO ratio cannot be larger than 104% 
            DOratioz(idx) = 1.05;
            DOz(idx) = DOratioz(idx).*DOsatz(idx); 
        else 
        end % if (Iflw>0)

    else 
    end % if (river_inflow_switch==1) 

    TKE=C_shelter*Az(1)*sqrt(tau^3/rho(1))*(24*60*60*dt);  % Turbulent kinetic energy (J) over the whole lake 

    % wind mixing 
    WmixIndicator = 1; 
    while (WmixIndicator==1)
        d_rho = diff(rho);
        idx = find(d_rho>0); 
        if (isempty(idx)==0); % The water layers are not fully mixed 
            zb=idx(1);
            mld=dz*zb; % Mixed layer depth 
            dD = d_rho(zb); % Dendity difference 
            Zg = sum(Az(1:zb+1).*zz(1:zb+1))/sum(Az(1:zb+1)); % Depth of mass center of mixed layer 
            V_weight=Vz(zb+1)*sum(Vz(1:zb))/(Vz(zb+1)+sum(Vz(1:zb))); 
            POE = dD*g*V_weight*(mld+dz/2-Zg); % Potential energy 
            KP_ratio = TKE/POE; 
            if (KP_ratio>=1)
                % Wind mixing DO (simple method)
                DOsatz = SolDO(Tz);
                DOratioz = DOz./DOsatz; 
                idx = find(DOratioz(1:zb+1)<0.9); % DO concentration in mixed layer cannot less than 90% 
                DOratioz(idx) = 0.9;
                DOz(idx) = DOsatz(idx).*DOratioz(idx); 
                DOmix=sum(Vz(1:zb+1).*DOz(1:zb+1))/sum(Vz(1:zb+1));
                DOz(1:zb+1)=DOmix; 
                % (Test for high DO)
                idx = find(DOratioz(1:zb+1)<0.9); % DO concentration in mixed layer cannot less than 90% 
                DOratioz(idx) = 0.9;
                DOz(idx) = DOsatz(idx).*DOratioz(idx); 
                PO4mix=sum(Vz(1:zb+1).*PO4z(1:zb+1))/sum(Vz(1:zb+1));
                PO4z(1:zb+1)=PO4mix;

                NOxmix=sum(Vz(1:zb+1).*NO3z(1:zb+1))/sum(Vz(1:zb+1));
                NO3z(1:zb+1)=NOxmix;

                NH3mix=sum(Vz(1:zb+1).*NH3z(1:zb+1))/sum(Vz(1:zb+1));
                NH3z(1:zb+1)=NH3mix;

                CyBmix=sum(Vz(1:zb+1).*CyBz(1:zb+1))/sum(Vz(1:zb+1));
                CyBz(1:zb+1)=CyBmix;

                APSBmix=sum(Vz(1:zb+1).*APSBz(1:zb+1))/sum(Vz(1:zb+1));
                APSBz(1:zb+1)=APSBmix;
                   
                H2Smix=sum(Vz(1:zb+1).*H2Sz(1:zb+1))/sum(Vz(1:zb+1));
                H2Sz(1:zb+1)=H2Smix;

                SO4mix=sum(Vz(1:zb+1).*SO4z(1:zb+1))/sum(Vz(1:zb+1));
                SO4z(1:zb+1)=SO4mix;

                DFemix=sum(Vz(1:zb+1).*DFez(1:zb+1))/sum(Vz(1:zb+1));
                DFez(1:zb+1)=DFemix;

                Fe3mix=sum(Vz(1:zb+1).*Fe3z(1:zb+1))/sum(Vz(1:zb+1));
                Fe3z(1:zb+1)=Fe3mix;

                DMnmix=sum(Vz(1:zb+1).*DMnz(1:zb+1))/sum(Vz(1:zb+1));
                DMnz(1:zb+1)=DMnmix;

                PMnmix=sum(Vz(1:zb+1).*PMnz(1:zb+1))/sum(Vz(1:zb+1));
                PMnz(1:zb+1)=PMnmix;

                rho = polyval(ies80,max(Tz(:),0)) + min(0,Tz(:));
                TKE = TKE - POE; 
   
            else % if KP_ratio < 1, then mix the TKE part of the underlying layer 
                DOsatz = SolDO(Tz);
                DOratioz = DOz./DOsatz; 
                idx = find(DOratioz(1:zb+1)<0.9); % DO concentration in mixed layer cannot less than 90% 
                DOratioz(idx) = 0.9;
                DOz(idx) = DOsatz(idx).*DOratioz(idx); 
                DOmix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*DOz(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                DOz(1:zb)=DOmix; 
                DOz(zb+1)=KP_ratio*DOmix + (1-KP_ratio)*DOz(zb+1); 
                % (Test for high DO) 
                idx = find(DOratioz(1:zb+1)<0.9); % DO concentration in mixed layer cannot less than 90% 
                DOratioz(idx) = 0.9;
                DOz(idx) = DOsatz(idx).*DOratioz(idx); 

                PO4mix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*PO4z(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                PO4z(1:zb)=PO4mix; 
                PO4z(zb+1)=KP_ratio*PO4mix + (1-KP_ratio)*PO4z(zb+1); 

                NOxmix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*NO3z(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                NO3z(1:zb)=NOxmix; 
                NO3z(zb+1)=KP_ratio*NOxmix + (1-KP_ratio)*NO3z(zb+1); 

                NH3mix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*NH3z(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                NH3z(1:zb)=NH3mix; 
                NH3z(zb+1)=KP_ratio*NH3mix + (1-KP_ratio)*NH3z(zb+1); 

                CyBmix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*CyBz(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                CyBz(1:zb)=CyBmix; 
                CyBz(zb+1)=KP_ratio*CyBmix + (1-KP_ratio)*CyBz(zb+1); 

                APSBmix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*APSBz(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                APSBz(1:zb)=APSBmix; 
                APSBz(zb+1)=KP_ratio*APSBmix + (1-KP_ratio)*APSBz(zb+1); 

                H2Smix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*H2Sz(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                H2Sz(1:zb)=H2Smix; 
                H2Sz(zb+1)=KP_ratio*H2Smix + (1-KP_ratio)*H2Sz(zb+1); 

                SO4mix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*SO4z(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                SO4z(1:zb)=SO4mix; 
                SO4z(zb+1)=KP_ratio*SO4mix + (1-KP_ratio)*SO4z(zb+1); 

                DFemix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*DFez(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                DFez(1:zb)=DFemix; 
                DFez(zb+1)=KP_ratio*DFemix + (1-KP_ratio)*DFez(zb+1); 

                Fe3mix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*Fe3z(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                Fe3z(1:zb)=Fe3mix; 
                Fe3z(zb+1)=KP_ratio*Fe3mix + (1-KP_ratio)*Fe3z(zb+1); 

                DMnmix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*DMnz(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                DMnz(1:zb)=DMnmix; 
                DMnz(zb+1)=KP_ratio*DMnmix + (1-KP_ratio)*DMnz(zb+1); 

                PMnmix=sum([Vz(1:zb);KP_ratio*Vz(zb+1)].*PMnz(1:zb+1))/sum([Vz(1:zb);KP_ratio*Vz(zb+1)]);
                PMnz(1:zb)=PMnmix; 
                PMnz(zb+1)=KP_ratio*PMnmix + (1-KP_ratio)*PMnz(zb+1); 

                rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0); 
                TKE = 0;
                WmixIndicator=0;
            end % if (KP_ratio>=1) 
            WmixIndicator=0; 
        end % if (isempty(idx)==0); the water column is not fully mixed 
    end % while 

    % Calculate pycnocline depth 
    pycno_thres=0.1; % Threshold density gradient value (kg m-3 m-1)
    rho = polyval(ies80,max(0,Tz(:))) + min(0,Tz(:));
    drho = [NaN;abs(diff(rho))];
    di=find((drho < pycno_thres*dz) | isnan(drho));
    drho(di) = 0;
    pycno_z = sum(zz.*drho)./sum(drho);

    % Output matrices 
    Qst(:,i) = [Qsw Qlw Qsl]';
    Kzt(:,i) = [0;Kz];
    Tzt(:,i) = Tz;
    rhozt(:,i) = rho; 
    DOzt(:,i) = DOz;
    pHzt(:,i) = pHz;
    PO4zt(:,i) = PO4z; 
    NO3zt(:,i) = NO3z;
    NH3zt(:,i) = NH3z; 
    CyBzt(:,i) = CyBz;
    APSBzt(:,i) = APSBz;
    H2Szt(:,i) = H2Sz;
    SO4zt(:,i) = SO4z;
    DFezt(:,i) = DFez;
    PFezt(:,i) = Fe3z+MFez;
    DMnzt(:,i) = DMnz;
    PMnzt(:,i) = PMnz;
    PARzt(:,i) = PARz;
    G_cybzt(:,i) = R_bio_cybz;
    G_apsbzt(:,i) = R_bio_apsbz;

end

runtime = toc; 

disp(['Total model runtime: ' int2str(floor(runtime/60)) ' min ' int2str(round(mod(runtime,60))) ' s']); 

% >>>>>> End of the time loop >>>>>>

% Below are the two functions for calculating tridiagonal matrix Fi for solving the 
% 1) diffusion equation (tridiag_DIF_v11), and 
% 2) advection-diffusion equation (tridiag_HAD_v11) by fully implicit hybrid exponential numerical scheme, 
% based on Dhamotharan et al. 1981, 
%'Unsteady one-dimensional settling of suspended sediments', Water Resources Research 17(4), 1125-1132
% code checked by TSA, 16.03.2004


%Inputs:
% Kz    diffusion coefficient at layer interfaces (plus surface) N (N,1)
% U     vertical settling velocity (scalar)
% Vz    layer volumes (N,1)
% Az    layer interface areas (N,1)
% dz    grid size
% dt    time step

%Output:
% Fi    tridiagonal matrix for solving new profile Cz


%=== DIFFUSIVE EQUATION ===
function[Fi]=tridiag_DIF_v11(Kz,Vz,Az,dz,dt)

Nz=length(Vz); %number of grid points/layers

% Linearized heat conservation equation matrix (diffusion only)
az = (dt/dz) * Kz .* (Az ./ Vz);                                        %coefficient for i-1
cz = (dt/dz) * [Kz(2:end); NaN] .* ([Az(2:end); NaN] ./ Vz);            %coefficient for i+1
bz = 1 + az + cz;                                                       %coefficient for i+1

%Boundary conditions, surface

az(1) = 0;
%cz(1) remains unchanged 
bz(1)= 1 + az(1) + cz(1);
   

%Boundary conditions, bottom

%az(end) remains unchanged 
cz(end) = 0;
bz(end) = 1 + az(end) + cz(end);
  
Gi = [-cz bz -az];
Fi = spdiags(Gi,-1:1,Nz,Nz)';
%end of function


%=== ADVECTIVE-DIFFUSIVE EQUATION ===
function[Fi]=tridiag_HAD_v11(Kz,U,Vz,Az,dz,dt)

if (U<0)
    error('only positive (downward) velocities allowed')
end

if (U==0)
    U=eps; %set Vz next to nothing (=2.2204e-016) in order to avoid division by zero
end

Nz=length(Vz); %number of grid points/layers

theta=U*(dt/dz);

az = theta.*(1 + (1./(exp( (U*Vz)./(Kz.*Az) ) - 1)));                   %coefficient for i-1
cz = theta./(exp( (U*Vz)./([Kz(2:end); NaN].*[Az(2:end); NaN]) ) - 1);  %coefficient for i+1
bz = 1 + az + cz;                                                       %coefficient for i

%Boundary conditions, surface

az(1) = 0;
%cz(1) remains unchanged 
bz(1) = 1 + theta + cz(1);

%Boundary conditions, bottom

%az(end) remains unchanged 
cz(end) = 0;
bz(end) = 1 + az(end);

Gi = [-cz bz -az];
Fi = spdiags(Gi,-1:1,Nz,Nz)';
%end of function


% Below is a function for calculating the partitioning between
% dissolved and inorganic particle bound phosphorus. 
%================================

function [DIP, PIPf]=Ppart(vf,TIP,Psat,Fmax,rho_sed,Fstable)
% Function for calculating the partitioning between
% dissolved and inorganic particle bound phosphorus. 
% Based on Langmuir isotherm approach 
%vf:    volume fraction of suspended inorganic matter (m3 m-3); S/rho_sed OR (1-porosity)
%TIP:   Total inorganic phosphorus (mg m-3)
%Psat, mg m-3 - Langmuir half-saturation parameter
%Fmax, mg kg-1  - Langmuir scaling parameter
%rho_sed, kg m-3 - Density of dry inorganic sediment mass
%Fstable, mg kg-1 - Inactive P conc. in inorg. particles

N=length(TIP);
DIP=NaN*ones(N,1);

for w=1:N
 a = vf(w)-1;
 b = TIP(w) + (vf(w)-1)*Psat - vf(w)*rho_sed*(Fmax+Fstable);
 c = Psat*TIP(w) - vf(w)*rho_sed*Fstable*Psat ;
DIP(w) = max(real(roots([a b c])));
end

PIPf = (TIP - (1-vf).*DIP)./(rho_sed*vf); %inorg. P conc. in sediment particles(mg kg-1 dry w.) 
%end of function 


% Below is a function for calculating dissolved oxygen solubility of a
% given temperature 

function [satDO] = SolDO(T)
% Source from: https://www.colby.edu/chemistry/CH331/O2%20Solubility.html
% Function for calculating dissolved oxygen (DO) solubility of a given
% temperature 
% T: Temperature (deg C) 

TK = T + 273.15; % Transfer the temperature to kelvin (K) 
first = -173.9894 + (255.5907*100./TK); 
second = first + (146.4813 * log(TK/100)) - (22.204*(TK/100)); 
% sal = sal*(-0.037362 + 0.016504*TK/100-0.0020564*(tempr/100)^2); % Salinity contribution component 
sal = 0; % Default salinity contribution 
third = second + sal; 
fourth = exp(third); 
pre = 760; % Default pressure (mmHg) 
sol = fourth*pre/760; % Get the solubility in uM 
satDO = sol/1000; % DO solubility in mol m-3 
% end of function 
