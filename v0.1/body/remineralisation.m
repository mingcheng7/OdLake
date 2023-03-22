% === One-day lake simulation model === % 
% version 0.1 
% By Ming Cheng 
% Last modification on 30/11/2022
% Remineralisation module 

function [Bacz,DOz,NH3z,NO3z,PO4z,SO4z,H2Sz,DFez,Fe3z,DMnz,PMnz] = remineralisation(dt,zz,Tz,Bacz,DOz,NH3z,NO3z,PO4z,SO4z,H2Sz,DFez,Fe3z,DMnz,PMnz,m_p,r_OM,theta)

% Variables 
% dt     : Time grid (day) 
% zz     : Depth range (m) 
% Bacz   : Bacteria biomass (cells m-3)
% DOz    : Dissolved oxygen distribution (mol m-3) 
% NH3z   : NH3 distribution (mol m-3)
% NO3z   : NO3 distribution (mol m-3)
% PO4z   : PO4 distribution (mol m-3)
% SO4z   : SO4 distribution (mol m-3)
% H2Sz   : H2S distribution (mol m-3)
% DFez   : Dissolved Fe distribution (mol m-3)
% PFez   : Particulate Fe distribution (mol m-3) 
% DMnz   : Dissolved Mn distribution (mol m-3) 
% PMnz   : Particulate Mn distribution (mol m-3)
% m_p    : Metabolism rate parameter (day-1) in 20 degree C. 
% r_OM   : OM parameter (mol OM cell-1) 
% theta  : Phytoplankton growth rate base 

Bac_prev = Bacz;  % An initial bacteria size 

% Set half-saturation parameters 
kNH3=1e-5;
kH2S=1e-3;
kDO=1e-3;

for i = 1:length(zz)
    % Mandatory nitrification and sulfurasation in oxic zone 
    v_N_O = dt*NH3z(i)/(NH3z(i)+kNH3)*DOz(i)/(DOz(i)+kDO);

    dNH3 = min(NH3z(i)*v_N_O,DOz(i)*v_N_O/2);
    dDO = 2*dNH3;
    dNO3 = dNH3;

    NH3z(i) = NH3z(i) - dNH3;
    DOz(i) = DOz(i) - dDO;
    NO3z(i) = NO3z(i) + dNO3;

    v_S_O = dt*H2Sz(i)/(H2Sz(i)+kH2S)*DOz(i)/(DOz(i)+kDO);

    dH2S = min(H2Sz(i)*v_S_O,DOz(i)*v_S_O/2);
    dDO = 2*dH2S;
    dSO4 = dH2S;

    H2Sz(i) = H2Sz(i) - dH2S;
    DOz(i) = DOz(i) - dDO;
    SO4z(i) = SO4z(i) + dSO4;

    % Metabolism  
    m = m_p*theta.^(Tz(i)-20)*dt; % Total metabolism rate in i layer 

    % Oxic zone 
    dBac = min(Bac_prev(i)*m,DOz(i)/(150*r_OM));
    dDO = dBac*r_OM*150;
    dNO3 = dBac*r_OM*16;
    dPO4 = dBac*r_OM;

    Bacz(i) = Bacz(i) - dBac;
    DOz(i) = DOz(i) - dDO;
    NO3z(i) = NO3z(i) + dNO3;
    PO4z(i) = PO4z(i) + dPO4;

    if ((dBac/Bac_prev(i)) < m)  % Pending if oxygen is enough for metabolism 
        m = m - dBac/Bac_prev(i); % Metabolism rate in i layer 

        % Sulfate reduction  
        dBac1 = min(Bac_prev(i)*m*0.69,SO4z(i)/(59*r_OM));
        dSO4 = dBac1*r_OM*59;
        dH2S = dSO4;
        dNH3 = dBac1*r_OM*16;
        dPO4 = dBac1*r_OM;

        Bacz(i) = Bacz(i) - dBac1;
        SO4z(i) = SO4z(i) - dSO4;
        H2Sz(i) = H2Sz(i) + dH2S;
        NH3z(i) = NH3z(i) + dNH3;
        PO4z(i) = PO4z(i) + dPO4;
        
        % Denitrification zone 
        dBac = min(Bac_prev(i)*m*0.31,NO3z(i)/(104*r_OM));
        dNO3 = dBac*r_OM*104;
        dPO4 = dBac*r_OM;

        Bacz(i) = Bacz(i) - dBac;
        NO3z(i) = NO3z(i) - dNO3;
        PO4z(i) = PO4z(i) + dPO4;

        dBac = dBac + dBac1; 

        if((dBac/Bac_prev(i) < m))  % Pending if nitrate is enough for metabolism 
            m = m - dBac/Bac_prev(i); % Metabolism rate in i layer 
            
            % Manganese reduction zone 
            dBac = min(Bac_prev(i)*m,PMnz(i)/(260*r_OM));
            dPMn = dBac*r_OM*260;
            dDMn = dPMn;
            dPO4 = dBac*r_OM;

            Bacz(i) = Bacz(i) - dBac;
            PMnz(i) = PMnz(i) - dPMn;
            DMnz(i) = DMnz(i) + dDMn;
            PO4z(i) = PO4z(i) + dPO4;

            if((dBac/Bac_prev(i) < m))  % Pending if MnO2 is enough for metabolism 
                m = m - dBac/Bac_prev(i); % Metabolism rate in i layer 
            
                % Iron reduction zone 
                dBac = min(Bac_prev(i)*m,Fe3z(i)/(472*r_OM));
                dFe3 = dBac*r_OM*472;
                dDFe = dFe3;
                dNH3 = dBac*r_OM*16;
                dPO4 = dBac*r_OM;

                Bacz(i) = Bacz(i) - dBac;
                Fe3z(i) = Fe3z(i) - dFe3;
                DFez(i) = DFez(i) + dDFe;
                NH3z(i) = NH3z(i) + dNH3;
                PO4z(i) = PO4z(i) + dPO4;

                if((dBac/Bac_prev(i) < m))  % Pending if Fe2O3 is enough for metabolism 
                    m = m - dBac/Bac_prev(i); % Metabolism rate in i layer 
            
                    % Sulfate reduction zone 
                    dBac = min(Bac_prev(i)*m,SO4z(i)/(59*r_OM));
                    dSO4 = dBac*r_OM*59;
                    dH2S = dSO4;
                    dNH3 = dBac*r_OM*16;
                    dPO4 = dBac*r_OM;

                    Bacz(i) = Bacz(i) - dBac;
                    SO4z(i) = SO4z(i) - dSO4;
                    H2Sz(i) = H2Sz(i) + dH2S;
                    NH3z(i) = NH3z(i) + dNH3;
                    PO4z(i) = PO4z(i) + dPO4;

                    if((dBac/Bac_prev(i) < m))  % Pending if sulfate is enough for metabolism 
                        m = m - dBac/Bac_prev(i); % Metabolism rate in i layer 
            
                        % Methane fermentation zone 
                        dBac = Bac_prev(i)*m;
                        dNH3 = dBac*r_OM*16;
                        dPO4 = dBac*r_OM;

                        Bacz(i) = Bacz(i) - dBac;
                        NH3z(i) = NH3z(i) + dNH3;
                        PO4z(i) = PO4z(i) + dPO4;
                    end
                end
            end
        end
    end
end 

