% This script is used to do light attenuation coefficient experiment. 
% 
% By Ming Cheng 
% 
% Last modified on 18 AUG 2022 
% 
% This experiment uses different light attenuation coefficient (ε) to show
% the effect of ε variation on cyanobacteria and APSB distribution. 
% 
% There are 5 experiment groups (different inflow temperature): control:
% 0.1/m, EXP 1: 0.08/m, EXP 2: 0.06/m, EXP 3: 0.15/m, and EXP 4: 0.2/m 
% 
clear all; clc; 

%% Set up each experiment 

Datetime = "27-Jun-2015 00:00:00";

Initfile = "Init_data.xlsx";

Initsheet = 1; 

Inputfile = "Input_data.xlsx";

Inputsheet = 1; 

Parafile = "Para_data_attenuation_coefficient_experiment.xlsx";

Parasheet1 = 1; % Control 

Parasheet2 = 2; % ε-

Parasheet3 = 3; % ε--

Parasheet4 = 4; % ε+

Parasheet5 = 5; % ε++

%% Do experiments 

% Experiment results are saved with an extra number 

% Control (0.1)

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt1,pHzt,PO4zt,NO3zt,NH3zt,CyBzt1,APSBzt1,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt1,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet1);

% Experiment 1 (0.08)  

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt2,pHzt,PO4zt,NO3zt,NH3zt,CyBzt2,APSBzt2,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt2,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet2);

% Experiment 2 (0.06) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt3,pHzt,PO4zt,NO3zt,NH3zt,CyBzt3,APSBzt3,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt3,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet3);

% Experiment 3 (0.15) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt4,pHzt,PO4zt,NO3zt,NH3zt,CyBzt4,APSBzt4,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt4,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet4);

% Experiment 4 (0.2) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt5,pHzt,PO4zt,NO3zt,NH3zt,CyBzt5,APSBzt5,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt5,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet5);

%% Analyse the experiment results 

% Obtain the final profiles of each experiment and initial condition for
% both cyanobacteria and APSB 

PARz1 = PARzt1(:,1)*1000000;

PARz2 = PARzt2(:,1)*1000000;

PARz3 = PARzt3(:,1)*1000000;

PARz4 = PARzt4(:,1)*1000000;

PARz5 = PARzt5(:,1)*1000000;

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

set(gcf,'Units','pixels','position',[100 200 2400 900]) 

subplot(141)

plot(log10(PARz1),zz,'r','LineWidth',2); hold on; 

plot(log10(PARz2),zz,'b','LineWidth',2)

plot(log10(PARz3),zz,'g','LineWidth',2)

plot(log10(PARz4),zz,'m','LineWidth',2)

plot(log10(PARz5),zz,'c','LineWidth',2)

set(gca,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',12,'XLim',[-4,3],...
    'XTick',-4:3,'XTickLabel',{'0.0001' '0.001' '0.01' '0.1' '1' '10' '100' '1000'})

xlabel('PAR (µmol s^{-1} m^{-2})','Fontsize',18,'FontWeight','bold')

ylabel('Depth (m)','FontSize',18,'FontWeight','bold')

text(-4.2,-2,'A','FontSize',18,'FontWeight','bold')

subplot(142)

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

text(-2000,-2,'B','FontSize',18,'FontWeight','bold')

subplot(143)

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

text(-20000,-2,'C','FontSize',18,'FontWeight','bold')

lgd = legend({'Initial' 'Control' 'ε-' 'ε--' 'ε+' 'ε++'},'FontSize',14);

title(lgd,'Treatment','FontSize',14,'FontWeight','bold')

subplot(144)

plot(DOz0,zz,'k','LineWidth',2); hold on; 

plot(DOz1,zz,'r','LineWidth',2)

plot(DOz2,zz,'b','LineWidth',2)

plot(DOz3,zz,'g','LineWidth',2)

plot(DOz4,zz,'m','LineWidth',2)

plot(DOz5,zz,'c','LineWidth',2)

set(gca,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',12)

xlabel('Dissolved oxygen (μmol L^{-1})','Fontsize',18,'FontWeight','bold')

ylabel('Depth (m)','FontSize',18,'FontWeight','bold')

text(-12,-2,'D','FontSize',18,'FontWeight','bold')

saveas(gcf,'Results\AC_experiment.jpg')

