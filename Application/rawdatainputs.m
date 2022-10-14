% OdLake model 
% by Ming Cheng 
% Module for input collecting data 

function[Z,T,PAR,DO,pH,PO4,NOx,NH3,CyB,APSB,H2S,SO4,DFe,PFe,DMn,PMn]=rawdatainputs(raw_filename,raw_sheet)

% Inputs 
% Input filename and sheet 

% Output 
% Z  : Depth of each measuring point (m) 
% T  : Temperature profile in the raw dataset (degree C)
% PAR: PAR (photosynthetically active radiation) profile in the raw dataset (µmol.s-1.m-2)
% DO : Dissolved oxygen profile in the raw dataset (µmol/L)
% pH : pH profile in the raw dataset (-)
% PO4: PO4 profile in the raw dataset (µmol/L)
% NOx: Nitrate + nitrite profile in the raw dataset (µmol/L)
% NH3: Ammonia + ammonium profile in the raw dataset (µmol/L) 
% CyB: Cyanobacteria profile in the raw dataset (cells/mL)
% APSB: Aerobic phototrophic sulfur bacteria profile in the raw dataset (cells/mL) 
% H2S: H2S profile in the raw dataset (µmol/L)
% SO4: SO4 profile in the raw dataset (mmol/L) 
% DFe: Dissolved Fe profile in the raw dataset (µmol/L) 
% PFe: Particulate Fe profile in the raw dataset (µmol/L)
% DMn: Dissolved Mn profile in the raw dataset (µmol/L) 
% PMn: Particulate Mn profile in the raw dataset (µmol/L)

[RawMx,StrMx]=xlsread(raw_filename,raw_sheet);
Z=RawMx(3:end,1);
T=RawMx(3:end,2);
PAR=RawMx(3:end,3);
DO=RawMx(3:end,4);
pH=RawMx(3:end,5);
PO4=RawMx(3:end,6);
NOx=RawMx(3:end,7);
NH3=RawMx(3:end,8);
CyB=RawMx(3:end,9);
APSB=RawMx(3:end,10);
H2S=RawMx(3:end,11);
SO4=RawMx(3:end,12);
DFe=RawMx(3:end,13);
PFe=RawMx(3:end,14);
DMn=RawMx(3:end,15);
PMn=RawMx(3:end,16);
