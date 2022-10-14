% This script is used to do inflow temperature experiment. 
% 
% By Ming Cheng 
% 
% Last modified on 18 AUG 2022 
% 
% Inflow temperature can affect the inflow density, impacting the input
% position of the lake. 
% 
% This experiment uses different inflow temperature (IT) to show the effect
% of input nutrient position variation on cyanobacteria and APSB
% distribution. 
% 
% There are 5 experiment groups (different inflow temperature): % control:
% 6 degree C (9 m) , EXP 1: 4.7 degree C (14 m), EXP 2: 4.5 degree C
% (bottom), EXP 3: 10 degree C (4.5 m), and EXP 4: 14 degree C (surface) 
% 
clear all; clc; 

%% Set up each experiment 

Datetime = "27-Jun-2015 00:00:00";

Initfile = "Init_data.xlsx";

Initsheet = 1; 

Inputfile = "Input_data_inflow_temperature_experiment.xlsx";

Parafile = "Para_data.xlsx";

Parasheet = 1; 

Inputsheet1 = 1; % Control 

Inputsheet2 = 2; % IV-

Inputsheet3 = 3; % IV--

Inputsheet4 = 4; % IV+

Inputsheet5 = 5; % IV++

%% Do experiments 

% Experiment results are saved with an extra number 

% Control 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt1,pHzt,PO4zt,NO3zt,NH3zt,CyBzt1,APSBzt1,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet1,Parafile,Parasheet);

% Experiment 1 (growth rate 0.15)  

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt2,pHzt,PO4zt,NO3zt,NH3zt,CyBzt2,APSBzt2,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet2,Parafile,Parasheet);

% Experiment 2 (growth rate 0.1) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt3,pHzt,PO4zt,NO3zt,NH3zt,CyBzt3,APSBzt3,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet3,Parafile,Parasheet);

% Experiment 3 (growth rate 0.25) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt4,pHzt,PO4zt,NO3zt,NH3zt,CyBzt4,APSBzt4,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet4,Parafile,Parasheet);

% Experiment 4 (growth rate 0.3) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt5,pHzt,PO4zt,NO3zt,NH3zt,CyBzt5,APSBzt5,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet5,Parafile,Parasheet);

%% Analyse the experiment results 

% Obtain the final profiles of each experiment and initial condition for
% both cyanobacteria and APSB 

CyBz0 = CyBzt1(:,1)/1000000; 

CyBz1 = CyBzt1(:,end)/1000000; 

CyBz2 = CyBzt2(:,end)/1000000; 

CyBz3 = CyBzt3(:,end)/1000000; 

CyBz4 = CyBzt4(:,end)/1000000; 

CyBz5 = CyBzt5(:,end)/1000000; 

APSBz0 = APSBzt1(:,1)/1000000;

APSBz1 = APSBzt1(:,end)/1000000;

APSBz2 = APSBzt2(:,end)/1000000;

APSBz3 = APSBzt3(:,end)/1000000;

APSBz4 = APSBzt4(:,end)/1000000;

APSBz5 = APSBzt5(:,end)/1000000;

DOz0 = DOzt1(:,1)*1000;

DOz1 = DOzt1(:,end)*1000; 

DOz2 = DOzt2(:,end)*1000; 

DOz3 = DOzt3(:,end)*1000; 

DOz4 = DOzt4(:,end)*1000; 

DOz5 = DOzt5(:,end)*1000; 

%% Plot the results 

figure(1)

set(gcf,'Units','pixels','position',[100 200 1700 900]) 

subplot(131)

plot(CyBz0,zz,'k','LineWidth',2); hold on; 

plot(CyBz1,zz,'r','LineWidth',2)

plot(CyBz2,zz,'b','LineWidth',2)

plot(CyBz3,zz,'g','LineWidth',2)

plot(CyBz4,zz,'m','LineWidth',2)

plot(CyBz5,zz,'c','LineWidth',2)

set(gca,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',12,'XTick',0:20000:80000, ...
    'XTickLabel',{'0' '20000' '40000' '60000' '80000'})

xlabel('Cyanobacteria (cells mL^{-1})','Fontsize',18,'FontWeight','bold')

ylabel('Depth (m)','FontSize',18,'FontWeight','bold')

text(-2000,-1.5,'A','FontSize',18,'FontWeight','bold')

lgd = legend({'Initial' 'Control' 'IT-' 'IT--' 'IT+' 'IT++'},'FontSize',14,'Location','southeast');

title(lgd,'Treatment','FontSize',14,'FontWeight','bold')

subplot(132)

plot(APSBz0,zz,'k','LineWidth',2); hold on; 

plot(APSBz1,zz,'r','LineWidth',2)

plot(APSBz2,zz,'b','LineWidth',2)

plot(APSBz3,zz,'g','LineWidth',2)

plot(APSBz4,zz,'m','LineWidth',2)

plot(APSBz5,zz,'c','LineWidth',2)

set(gca,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',12,'XTick',0:200000:800000, ...
    'XTickLabel',{'0' '200000' '400000' '600000' '800000'})

xlabel('APSB (cells mL^{-1})','Fontsize',18,'FontWeight','bold')

ylabel('Depth (m)','FontSize',18,'FontWeight','bold')

text(-20000,-1.5,'B','FontSize',18,'FontWeight','bold')

subplot(133)

plot(DOz0,zz,'k','LineWidth',2); hold on; 

plot(DOz1,zz,'r','LineWidth',2)

plot(DOz2,zz,'b','LineWidth',2)

plot(DOz3,zz,'g','LineWidth',2)

plot(DOz4,zz,'m','LineWidth',2)

plot(DOz5,zz,'c','LineWidth',2)

set(gca,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',12)

xlabel('Dissolved oxygen (μmol L^{-1})','Fontsize',18,'FontWeight','bold')

ylabel('Depth (m)','FontSize',18,'FontWeight','bold')

text(-12,-1.5,'C','FontSize',18,'FontWeight','bold')

saveas(gcf,'Results\IT_experiment.jpg')

