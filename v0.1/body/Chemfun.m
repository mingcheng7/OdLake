% === One-day lake simulation model === % 
% version 0.1 
% By Ming Cheng 
% Last modification on 30/11/2022
% 
% Chemical reactions module 


function [DOz,H2Sz,SO4z,DFez,Fe3z,MFez,DMnz,PMnz] = Chemfun(dt,zz,DOz,H2Sz,SO4z,DFez,Fe3z,MFez,DMnz,PMnz)

% Variables: 
% DOz : Dissolved oxygen distribution (mol m-3) 
% H2Sz: H2S distribution (mol m-3)
% SO4z: SO4 distribution (mol m-3)
% DFez: Dissolved Fe(II) distribution (mol m-3)
% Fe3z: Particulate Fe(III) (FeOOH, Fe(OH)3,Fe2O3...) distribution (mol m-3)
% MFez: Fe(II)-sulfide minerals (FeS,FeS2...) distribution (mol m-3)
% DMnz: Dissolved Mn(II) distribution (mol m-3)
% PMnz: Particulate Manganese dioxides (MnO2) distribution (mol m-3) 

% Set half-saturation parameters  
kDO = 1e-3; 
kH2S = 1e-3;
kDFe = 1e-4;
kFe3 = 1e-5; 
kDMn = 1e-3;
kPMn = 1e-3;

% Set Fe(II) half saturation concentration on forming pyrite 
Fe_half = 1e-4; 

% FeS formation rate 
p0 = 0.0005; 
% Set default elemental S 
S = 0;

% Depth iteration 
for i =1:length(zz)
    % General redox reactions 
    % Oxic zone 
    v_Mn_O = dt*DMnz(i)/(DMnz(i)+kDMn)*DOz(i)/(DOz(i)+kDO);

    dDMn = min(DMnz(i)*v_Mn_O,DOz(i)*v_Mn_O*2);
    dDO = 0.5*dDMn;
    dPMn = dDMn;

    DMnz(i) = DMnz(i) - dDMn;
    DOz(i) = DOz(i) - dDO;
    PMnz(i) = PMnz(i) + dPMn;

    v_Fe_O = dt*DFez(i)/(DFez(i)+kDFe)*DOz(i)/(DOz(i)+kDO);

    dDFe = min(DFez(i)*v_Fe_O,DOz(i)*v_Fe_O*4);
    dDO = 0.25*dDFe;
    dFe3 = dDFe;

    DFez(i) = DFez(i) - dDFe;
    DOz(i) = DOz(i) - dDO;
    Fe3z(i) = Fe3z(i) + dFe3;
        
    % Chemicoline 
    v_Fe_Mn = dt*DFez(i)/(DFez(i)+kDFe)*PMnz(i)/(PMnz(i)+kPMn);

    dDFe = min(DFez(i)*v_Fe_Mn,PMnz(i)*v_Fe_Mn*2);
    dPMn = 0.5*dDFe;
    dFe3 = dDFe;
    dDMn = dPMn; 

    DFez(i) = DFez(i) - dDFe;
    PMnz(i) = PMnz(i) - dPMn;
    Fe3z(i) = Fe3z(i) + dFe3;
    DMnz(i) = DMnz(i) + dDMn; 

    % Anoxic zone 
    v_Mn_S = dt*PMnz(i)/(PMnz(i)+kPMn)*H2Sz(i)/(H2Sz(i)+kH2S);

    dPMn = min(PMnz(i)*v_Mn_S,H2Sz(i)*v_Mn_S);
    dH2S = dPMn;
    dDMn = dPMn;
    dS = dH2S;

    PMnz(i) = PMnz(i) - dPMn;
    H2Sz(i) = H2Sz(i) - dH2S;
    DMnz(i) = DMnz(i) + dDMn;
    S = S + dS;
    
    v_Fe_S = dt*Fe3z(i)/(Fe3z(i)+kFe3)*H2Sz(i)/(H2Sz(i)+kH2S);

    dFe3 = min(Fe3z(i)*v_Fe_S,H2Sz(i)*v_Fe_S*2);
    dH2S = 0.5*dFe3;
    dDFe = dFe3;
    dS = dH2S;

    Fe3z(i) = Fe3z(i) - dFe3;
    H2Sz(i) = H2Sz(i) - dH2S;
    DFez(i) = DFez(i) + dDFe;
    S = S + dS; 

    % The pyrite formation in minimolimnion 
    if (Fe3z(i)<kFe3)
        v_FeS = p0*dt*DFez(i)/(DFez(i)+Fe_half)*H2Sz(i)/(H2Sz(i)+kH2S);
        
        dDFe = min(DFez(i),H2Sz(i))*v_FeS;
        dH2S = dDFe;
        dMFe = dDFe;

        DFez(i) = DFez(i) - dDFe;
        H2Sz(i) = H2Sz(i) - dH2S;
        MFez(i) = MFez(i) + dMFe;
    end 

    % Terminal of elemental S 
    dS = S;
    dH2S = 0.75*dS;
    dSO4 = 0.25*dS;

    S = S-dS;
    H2Sz(i) = H2Sz(i) + dH2S;
    SO4z(i) = SO4z(i) + dSO4;

end

