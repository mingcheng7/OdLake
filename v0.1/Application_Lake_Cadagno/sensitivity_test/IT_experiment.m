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

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt1,PO4zt,NO3zt1,NH3zt,CyBzt1,APSBzt1,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,PARzt,Wt]...
    = odlake_v0_1(Datetime,Initfile,Initsheet,Inputfile,Inputsheet1,Parafile,Parasheet);

% Experiment 1 (growth rate 0.15)  

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt2,PO4zt,NO3zt2,NH3zt,CyBzt2,APSBzt2,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,PARzt,Wt]...
    = odlake_v0_1(Datetime,Initfile,Initsheet,Inputfile,Inputsheet2,Parafile,Parasheet);

% Experiment 2 (growth rate 0.1) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt3,PO4zt,NO3zt3,NH3zt,CyBzt3,APSBzt3,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,PARzt,Wt]...
    = odlake_v0_1(Datetime,Initfile,Initsheet,Inputfile,Inputsheet3,Parafile,Parasheet);

% Experiment 3 (growth rate 0.25) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt4,PO4zt,NO3zt4,NH3zt,CyBzt4,APSBzt4,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,PARzt,Wt]...
    = odlake_v0_1(Datetime,Initfile,Initsheet,Inputfile,Inputsheet4,Parafile,Parasheet);

% Experiment 4 (growth rate 0.3) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt5,PO4zt,NO3zt5,NH3zt,CyBzt5,APSBzt5,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,PARzt,Wt]...
    = odlake_v0_1(Datetime,Initfile,Initsheet,Inputfile,Inputsheet5,Parafile,Parasheet);

%% Analyse the experiment results 

% Obtain the final profiles of each experiment and initial condition for
% both cyanobacteria and APSB 

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

NO3z1 = NO3zt1(:,end)*1000;

NO3z2 = NO3zt2(:,end)*1000;

NO3z3 = NO3zt3(:,end)*1000;

NO3z4 = NO3zt4(:,end)*1000;

NO3z5 = NO3zt5(:,end)*1000;

%% Plot the results 

figure(1)

set(gcf,'Units','pixels','position',[50 50 1025 400]) 

axesPosition = [50 50 200 300]; % Axes position, in pixels 

h1 = axes('Units','pixels','Position',axesPosition,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',7,'XTick',0:20000:80000, ...
    'XTickLabel',{'0' '20000' '40000' '60000' '80000'},'NextPlot','add');

plot(h1,CyBz1,zz,'r','LineWidth',1); hold on; 

plot(h1,CyBz2,zz,'b','LineWidth',1)

plot(h1,CyBz3,zz,'g','LineWidth',1)

plot(h1,CyBz4,zz,'m','LineWidth',1)

plot(h1,CyBz5,zz,'c','LineWidth',1)

xlabel('Cyanobacteria (cells mL^{-1})','Fontsize',9,'FontWeight','bold')

ylabel('Depth (m)','FontSize',9,'FontWeight','bold')

text(-2000,-2.55,'a','FontSize',9,'FontWeight','bold')

lgd = legend({'Control' 'IT-' 'IT--' 'IT+' 'IT++'},'FontSize',6,'Location','southeast');

title(lgd,'Treatment','FontSize',6,'FontWeight','bold')

box on;

axesPosition = [300 50 200 300]; % Axes position, in pixels 

h2 = axes('Units','pixels','Position',axesPosition,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',7,'XTick',0:200000:800000, ...
    'XTickLabel',{'0' '200000' '400000' '600000' '800000'},'NextPlot','add');

plot(APSBz1,zz,'r','LineWidth',1); hold on; 

plot(APSBz2,zz,'b','LineWidth',1)

plot(APSBz3,zz,'g','LineWidth',1)

plot(APSBz4,zz,'m','LineWidth',1)

plot(APSBz5,zz,'c','LineWidth',1)

xlabel('APSB (cells mL^{-1})','Fontsize',9,'FontWeight','bold')

ylabel('Depth (m)','FontSize',9,'FontWeight','bold')

text(-20000,-2.55,'b','FontSize',9,'FontWeight','bold')

box on;

axesPosition = [550 50 200 300]; % Axes position, in pixels 

h3 = axes('Units','pixels','Position',axesPosition,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',7,'NextPlot','add');

plot(DOz1,zz,'r','LineWidth',1); hold on; 

plot(DOz2,zz,'b','LineWidth',1)

plot(DOz3,zz,'g','LineWidth',1)

plot(DOz4,zz,'m','LineWidth',1)

plot(DOz5,zz,'c','LineWidth',1)

xlabel('Dissolved oxygen (μmol L^{-1})','Fontsize',9,'FontWeight','bold')

ylabel('Depth (m)','FontSize',9,'FontWeight','bold')

text(-12,-2.55,'c','FontSize',9,'FontWeight','bold')

box on;

axesPosition = [800 50 200 300]; % Axes position, in pixels 

h4 = axes('Units','pixels','Position',axesPosition,'YDir','reverse','XAxisLocation','top','YLim',[0 20],'Fontsize',7,'NextPlot','add');

plot(NO3z1,zz,'r','LineWidth',1); hold on; 

plot(NO3z2,zz,'b','LineWidth',1)

plot(NO3z3,zz,'g','LineWidth',1)

plot(NO3z4,zz,'m','LineWidth',1)

plot(NO3z5,zz,'c','LineWidth',1)

xlabel('NO_{3}^{-} (μmol L^{-1})','Fontsize',9,'FontWeight','bold')

ylabel('Depth (m)','FontSize',9,'FontWeight','bold')

text(-12,-2.55,'d','FontSize',9,'FontWeight','bold')

box on;

exportgraphics(gcf,'Results\Test_3.jpg','Resolution',600)

