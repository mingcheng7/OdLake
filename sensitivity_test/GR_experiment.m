% This script is used to do growth rate experiment. 
% 
% By Ming Cheng 
% 
% Last modified on 17 AUG 2022 
% 
% This experiment uses different cyanobacteria growth rate parameter to
% show the change in cyanobacteria growth and the correlated effect on APSB
% growth. 
% 
% There are 5 experiment groups (different cyanobacteria growth rate):
% control (0.2 day-1), EXP 1 (0.15), EXP 2 (0.1), EXP 3 (0.25), and EXP 4
% (0.3). 
% 
clear all; clc; 

%% Set up each experiment 

Datetime = "27-Jun-2015 00:00:00";

Initfile = "Init_data.xlsx";

Initsheet = 1; 

Inputfile = "Input_data.xlsx";

Inputsheet = 1;

Parafile = "Para_data_growth_rate_experiment.xlsx";

Parasheet1 = 1; % Control 

Parasheet2 = 2; % Experiment 1 (growth rate 0.15)  

Parasheet3 = 3; % Experiment 2 (growth rate 0.1) 

Parasheet4 = 4; % Experiment 3 (growth rate 0.25) 

Parasheet5 = 5; % Experiment 4 (growth rate 0.3) 

%% Do experiments 

% Experiment results are saved with an extra number 

% Control 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt1,APSBzt1,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet1);

% Experiment 1 (growth rate 0.15)  

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt2,APSBzt2,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet2);

% Experiment 2 (growth rate 0.1) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt3,APSBzt3,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet3);

% Experiment 3 (growth rate 0.25) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt4,APSBzt4,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet4);

% Experiment 4 (growth rate 0.3) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt5,APSBzt5,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet5);

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

set(gca,'FontSize',12,'YLim',[30000 80000],'YTick',30000:10000:80000,'YTickLabel',...
    {'30000' '40000' '50000' '60000' '70000' '80000'},'XTick',1:5,'XTickLabel',...
    {'Control' 'GR-' 'GR--' 'GR+' 'GR++'}) 

xlabel('Treatment','FontSize',18,'FontWeight','bold')

ylabel('Cyanobacteria peak population (cells mL^{-1})','FontSize',18,'FontWeight','bold')

text(0.4,83000,'A','FontSize',18,'FontWeight','bold')

subplot(122)

bar([APSBt1;APSBt2;APSBt3;APSBt4;APSBt5])

set(gca,'FontSize',12,'YLim',[0 700000],'YTick',0:100000:700000,'YTickLabel',...
    {'0' '100000' '200000' '300000' '400000' '500000' '600000' '700000'},...
    'XTick',1:5,'XTickLabel',{'Control' 'GR-' 'GR--' 'GR+' 'GR++'})

xlabel('Treatment','FontSize',18,'FontWeight','bold')

ylabel('APSB peak population (cells mL^{-1})','FontSize',18,'FontWeight','bold')

text(0.4,745000,'B','FontSize',18,'FontWeight','bold')

lgd = legend({'0d' '500d' '1000d' '1500d' '2000d' '2500d'},'Location','northeast','FontSize',12,'Color','none');

title(lgd,'Experiment days')

saveas(gcf,'Results\cyb_GR_experiment.jpg')


