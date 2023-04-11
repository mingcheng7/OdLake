% This script is used to do shortwave radiation experiment. 
% 
% By Ming Cheng 
% 
% Last modified on 18 AUG 2022 
% 
% This experiment uses different shortwave radiation to show the effect of
% light attenuation variation on cyanobacteria and APSB distribution. 
% 
% There are 5 experiment groups (different shortwave radiation): control 
% (simulated shortwave radiation), EXP 1: 80% control, EXP 2: 60% control, 
% EXP 3: 120% control, and EXP 4: 140% control  
% 
clear all; clc; 

%% Set up each experiment 

Datetime = "27-Jun-2015 00:00:00";

Initfile = "Init_data.xlsx";

Initsheet = 1; 

Inputfile = "Input_data.xlsx";

Inputsheet = 1;

Parafile = "Para_data.xlsx";

Parasheet = 1; 

%% Do experiments 

% Experiment results are saved with an extra number 

% Control 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt1,PO4zt,NO3zt,NH3zt,CyBzt1,APSBzt1,H2Szt,SO4zt,DFezt1,PFezt,DMnzt1,PMnzt,PARzt1,Wt]...
    = odlake_v0_1_SR_exp(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,1);

% Experiment 1 (80%)  

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt2,PO4zt,NO3zt,NH3zt,CyBzt2,APSBzt2,H2Szt,SO4zt,DFezt2,PFezt,DMnzt2,PMnzt,PARzt2,Wt]...
    = odlake_v0_1_SR_exp(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,0.8);

% Experiment 2 (60%) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt3,PO4zt,NO3zt,NH3zt,CyBzt3,APSBzt3,H2Szt,SO4zt,DFezt3,PFezt,DMnzt3,PMnzt,PARzt3,Wt]...
    = odlake_v0_1_SR_exp(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,0.6);

% Experiment 3 (120%) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt4,PO4zt,NO3zt,NH3zt,CyBzt4,APSBzt4,H2Szt,SO4zt,DFezt4,PFezt,DMnzt4,PMnzt,PARzt4,Wt]...
    = odlake_v0_1_SR_exp(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,1.2);

% Experiment 4 (140%) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt5,PO4zt,NO3zt,NH3zt,CyBzt5,APSBzt5,H2Szt,SO4zt,DFezt5,PFezt,DMnzt5,PMnzt,PARzt5,Wt]...
    = odlake_v0_1_SR_exp(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,1.4);

%% Analyse the experiment results 

% Obtain the final profiles of each experiment and initial condition for
% both cyanobacteria and APSB 

PARz1 = PARzt1(:,1)*1000000;

PARz2 = PARzt2(:,1)*1000000;

PARz3 = PARzt3(:,1)*1000000;

PARz4 = PARzt4(:,1)*1000000;

PARz5 = PARzt5(:,1)*1000000;

CyBz1 = CyBzt1(:,end)/1000000; 

CyBz2 = CyBzt2(:,end)/1000000; 

CyBz3 = CyBzt3(:,end)/1000000; 

CyBz4 = CyBzt4(:,end)/1000000; 

CyBz5 = CyBzt5(:,end)/1000000; 

APSBz1 = APSBzt1(:,end)/1000000;

APSBz2 = APSBzt2(:,end)/1000000;

APSBz3 = APSBzt3(:,end)/1000000;

APSBz4 = APSBzt4(:,end)/1000000;

APSBz5 = APSBzt5(:,end)/1000000;

DOz1 = DOzt1(:,end)*1000; 

DOz2 = DOzt2(:,end)*1000; 

DOz3 = DOzt3(:,end)*1000; 

DOz4 = DOzt4(:,end)*1000; 

DOz5 = DOzt5(:,end)*1000; 

DFez1 = DFezt1(:,end)*1000; 

DFez2 = DFezt2(:,end)*1000; 

DFez3 = DFezt3(:,end)*1000; 

DFez4 = DFezt4(:,end)*1000; 

DFez5 = DFezt5(:,end)*1000; 

DMnz1 = DMnzt1(:,end)*1000;

DMnz2 = DMnzt2(:,end)*1000;

DMnz3 = DMnzt3(:,end)*1000;

DMnz4 = DMnzt4(:,end)*1000;

DMnz5 = DMnzt5(:,end)*1000;

%% Plot the results 

figure(1)

set(gcf,'Units','pixels','position',[50 50 775 800]) 

% PAR 

axesPosition = [50 450 200 300]; % Axes position, in pixels 

h1 = axes('Units','pixels','Position',axesPosition,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',7,...
    'XTick',-3:3,'XTickLabel',{'0.001' '0.01' '0.1' '1' '10' '100' '1000'},'NextPlot','add');

plot(h1,log10(PARz1),zz,'r','LineWidth',1); hold on; 

plot(h1,log10(PARz2),zz,'b','LineWidth',1)

plot(h1,log10(PARz3),zz,'g','LineWidth',1)

plot(h1,log10(PARz4),zz,'m','LineWidth',1)

plot(h1,log10(PARz5),zz,'c','LineWidth',1)

xlabel('PAR (µmol s^{-1} m^{-2})','Fontsize',9,'FontWeight','bold')

ylabel('Depth (m)','FontSize',9,'FontWeight','bold')

text(-4.2,-2.55,'a','FontSize',9,'FontWeight','bold')

lgd = legend({'Control' 'Q_{sw}-' 'Q_{sw}--' 'Q_{sw}+' 'Q_{sw}++'},'FontSize',6,'Location','northwest');

title(lgd,'Treatment','FontSize',6,'FontWeight','bold') 

box on;

% Cyanobacteria 

axesPosition = [300 450 200 300]; % Axes position, in pixels 

h2 = axes('Units','pixels','Position',axesPosition,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',7,'XTick',0:20000:80000, ...
    'XTickLabel',{'0' '20000' '40000' '60000' '80000'},'NextPlot','add');

plot(h2,CyBz1,zz,'r','LineWidth',1); hold on; 

plot(h2,CyBz2,zz,'b','LineWidth',1)

plot(h2,CyBz3,zz,'g','LineWidth',1)

plot(h2,CyBz4,zz,'m','LineWidth',1)

plot(h2,CyBz5,zz,'c','LineWidth',1)

xlabel('Cyanobacteria (cells mL^{-1})','Fontsize',9,'FontWeight','bold')

ylabel('Depth (m)','FontSize',9,'FontWeight','bold')

text(-2000,-2.55,'b','FontSize',9,'FontWeight','bold')

box on;

% APSB

axesPosition = [550 450 200 300]; % Axes position, in pixels 

h3 = axes('Units','pixels','Position',axesPosition,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',7,'XTick',0:200000:800000, ...
    'XTickLabel',{'0' '200000' '400000' '600000' '800000'},'NextPlot','add');

plot(h3,APSBz1,zz,'r','LineWidth',1); hold on; 

plot(h3,APSBz2,zz,'b','LineWidth',1)

plot(h3,APSBz3,zz,'g','LineWidth',1)

plot(h3,APSBz4,zz,'m','LineWidth',1)

plot(h3,APSBz5,zz,'c','LineWidth',1)

xlabel('APSB (cells mL^{-1})','Fontsize',9,'FontWeight','bold')

ylabel('Depth (m)','FontSize',9,'FontWeight','bold')

text(-20000,-2.55,'c','FontSize',9,'FontWeight','bold')

box on;

% Dissolved oxygen

axesPosition = [50 50 200 300]; % Axes position, in pixels 

h4 = axes('Units','pixels','Position',axesPosition,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',7,'NextPlot','add');

plot(h4,DOz1,zz,'r','LineWidth',1); hold on; 

plot(h4,DOz2,zz,'b','LineWidth',1)

plot(h4,DOz3,zz,'g','LineWidth',1)

plot(h4,DOz4,zz,'m','LineWidth',1)

plot(h4,DOz5,zz,'c','LineWidth',1)

xlabel('Dissolved oxygen (μmol L^{-1})','Fontsize',9,'FontWeight','bold')

ylabel('Depth (m)','FontSize',9,'FontWeight','bold')

text(-12,-2.55,'d','FontSize',9,'FontWeight','bold')

box on;

% Dissolved Fe

axesPosition = [300 50 200 300]; % Axes position, in pixels 

h5 = axes('Units','pixels','Position',axesPosition,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',7,'NextPlot','add');

plot(h5,DFez1,zz,'r','LineWidth',1); hold on; 

plot(h5,DFez2,zz,'b','LineWidth',1)

plot(h5,DFez3,zz,'g','LineWidth',1)

plot(h5,DFez4,zz,'m','LineWidth',1)

plot(h5,DFez5,zz,'c','LineWidth',1)

xlabel('Dissolved Fe (μmol L^{-1})','Fontsize',9,'FontWeight','bold')

ylabel('Depth (m)','FontSize',9,'FontWeight','bold')

text(0,-2.55,'e','FontSize',9,'FontWeight','bold')

box on;

% Dissolved Mn

axesPosition = [550 50 200 300]; % Axes position, in pixels 

h6 = axes('Units','pixels','Position',axesPosition,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',7,'NextPlot','add');

plot(h6,DMnz1,zz,'r','LineWidth',1); hold on; 

plot(h6,DMnz2,zz,'b','LineWidth',1)

plot(h6,DMnz3,zz,'g','LineWidth',1)

plot(h6,DMnz4,zz,'m','LineWidth',1)

plot(h6,DMnz5,zz,'c','LineWidth',1)

xlabel('Dissolved Mn (μmol L^{-1})','Fontsize',9,'FontWeight','bold')

ylabel('Depth (m)','FontSize',9,'FontWeight','bold')

text(0,-2.55,'f','FontSize',9,'FontWeight','bold')

box on;

exportgraphics(gcf,'Results\Test_4.jpg','Resolution',600)


