% This script is used to do inflow volume experiment. 
% 
% By Ming Cheng 
% 
% Last modified on 18 AUG 2022 
% 
% This experiment uses different inflow volume (IV) to show the effect of
% input nutrient variation on cyanobacteria and APSB growth. 
% 
% There are 5 experiment groups (different inflow volume): % control (200
% m3/d), EXP 1 (100), EXP 2 (0), EXP 3 (300), and EXP 4 (500). 
% 
clear all; clc; 

%% Set up each experiment 

Datetime = "27-Jun-2015 00:00:00";

Initfile = "Init_data.xlsx";

Initsheet = 1; 

Inputfile = "Input_data_inflow_volume_experiment.xlsx";

Parafile = "Para_data_growth_rate_experiment.xlsx";

Parasheet = 1; 

Inputsheet1 = 1; % Control 

Inputsheet2 = 2; % IV-

Inputsheet3 = 3; % IV--

Inputsheet4 = 4; % IV+

Inputsheet5 = 5; % IV++

%% Do experiments 

% Experiment results are saved with an extra number 

% Control 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt1,APSBzt1,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet1,Parafile,Parasheet);

% Experiment 1 (growth rate 0.15)  

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt2,APSBzt2,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet2,Parafile,Parasheet);

% Experiment 2 (growth rate 0.1) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt3,APSBzt3,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet3,Parafile,Parasheet);

% Experiment 3 (growth rate 0.25) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt4,APSBzt4,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet4,Parafile,Parasheet);

% Experiment 4 (growth rate 0.3) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt5,APSBzt5,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet5,Parafile,Parasheet);

%% Analyse the experiment results 

% Calculate the max of the cyanobacteria and APSB population in few days
% (0, 500, 1000, 1500, 2000, 2500) 
CyBt1 = [max(CyBzt1(:,1)) max(CyBzt1(:,5000)) max(CyBzt1(:,10000)) max(CyBzt1(:,15000)) max(CyBzt1(:,20000)) max(CyBzt1(:,25000))]/1000000; 

CyBt2 = [max(CyBzt2(:,1)) max(CyBzt2(:,5000)) max(CyBzt2(:,10000)) max(CyBzt2(:,15000)) max(CyBzt2(:,20000)) max(CyBzt2(:,25000))]/1000000; 

CyBt3 = [max(CyBzt3(:,1)) max(CyBzt3(:,5000)) max(CyBzt3(:,10000)) max(CyBzt3(:,15000)) max(CyBzt3(:,20000)) max(CyBzt3(:,25000))]/1000000; 

CyBt4 = [max(CyBzt4(:,1)) max(CyBzt4(:,5000)) max(CyBzt4(:,10000)) max(CyBzt4(:,15000)) max(CyBzt4(:,20000)) max(CyBzt4(:,25000))]/1000000; 

CyBt5 = [max(CyBzt5(:,1)) max(CyBzt5(:,5000)) max(CyBzt5(:,10000)) max(CyBzt5(:,15000)) max(CyBzt5(:,20000)) max(CyBzt5(:,25000))]/1000000; 

APSBt1 = [max(APSBzt1(:,1)) max(APSBzt1(:,5000)) max(APSBzt1(:,10000)) max(APSBzt1(:,15000)) max(APSBzt1(:,20000)) max(APSBzt1(:,25000))]/1000000; 

APSBt2 = [max(APSBzt2(:,1)) max(APSBzt2(:,5000)) max(APSBzt2(:,10000)) max(APSBzt2(:,15000)) max(APSBzt2(:,20000)) max(APSBzt2(:,25000))]/1000000; 

APSBt3 = [max(APSBzt3(:,1)) max(APSBzt3(:,5000)) max(APSBzt3(:,10000)) max(APSBzt3(:,15000)) max(APSBzt3(:,20000)) max(APSBzt3(:,25000))]/1000000; 

APSBt4 = [max(APSBzt4(:,1)) max(APSBzt4(:,5000)) max(APSBzt4(:,10000)) max(APSBzt4(:,15000)) max(APSBzt4(:,20000)) max(APSBzt4(:,25000))]/1000000; 

APSBt5 = [max(APSBzt5(:,1)) max(APSBzt5(:,5000)) max(APSBzt5(:,10000)) max(APSBzt5(:,15000)) max(APSBzt5(:,20000)) max(APSBzt5(:,25000))]/1000000; 

%% Plot the results 

figure(1)

set(gcf,'Units','pixels','position',[100 100 1600 1000]) 

% Cyanobacteria 

subplot(121)

bar([CyBt1;CyBt2;CyBt3;CyBt4;CyBt5])

set(gca,'FontSize',12,'YLim',[0 120000],'YTick',0:20000:120000,'YTickLabel',...
    {'0' '20000' '40000' '60000' '80000' '100000' '120000'},'XTick',1:5,'XTickLabel',...
    {'Control' 'IV-' 'IV--' 'IV+' 'IV++'}) 

xlabel('Treatment','FontSize',18,'FontWeight','bold')

ylabel('Cyanobacteria peak population (cells mL^{-1})','FontSize',18,'FontWeight','bold')

text(0.4,126000,'A','FontSize',18,'FontWeight','bold')

subplot(122)

bar([APSBt1;APSBt2;APSBt3;APSBt4;APSBt5])

set(gca,'FontSize',12,'YLim',[0 1200000],'YTick',0:200000:1200000,'YTickLabel',...
    {'0' '200000' '400000' '600000' '800000' '1000000' '1200000'},...
    'XTick',1:5,'XTickLabel',{'Control' 'IV-' 'IV--' 'IV+' 'IV++'})

xlabel('Treatment','FontSize',18,'FontWeight','bold')

ylabel('APSB peak population (cells mL^{-1})','FontSize',18,'FontWeight','bold')

text(0.4,1260000,'B','FontSize',18,'FontWeight','bold')

lgd = legend({'0d' '500d' '1000d' '1500d' '2000d' '2500d'},'Location','northwest','FontSize',12,'Color','none');

title(lgd,'Experiment days')

saveas(gcf,'Results\IV_experiment.jpg')


