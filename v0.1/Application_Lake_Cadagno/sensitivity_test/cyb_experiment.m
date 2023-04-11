% This script is used to make a timeline that describes cyanobacteria
% growth process. 
% 
% Last modified on 17 AUG 2022 
% 
% Here, I will do 5 experiments to test the effect of different initial
% cyanobacteria population on the later growth and the final cyanobacteria
% population. In this experiment, I use the peak value of cyanobacteria
% population as the asseement standard. 
% 
% 5 experiments are: control, 80% compared to control, 60%..., 120%...
% 140%... . 

clear all; clc; 

%% Set up each experiment 

Datetime = "27-Jun-2015 00:00:00";

Initfile = "Init_data_cyb_experiment.xlsx";

Inputfile = "Input_data.xlsx";

Inputsheet = 1;

Parafile = "Para_data.xlsx";

Parasheet = 1;

Initsheet1 = 1; % Control 

Initsheet2 = 2; % Experiment 1 (80% initial cyanobacteria)  

Initsheet3 = 3; % Experiment 2 (60% initial cyanobacteria) 

Initsheet4 = 4; % Experiment 3 (120% initial cyanobacteria) 

Initsheet5 = 5; % Experiment 4 (140% initial cyanobacteria) 

%% Do experiments 

% Experiment results are saved with an extra number 

% Control 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,PO4zt,NO3zt,NH3zt,CyBzt1,APSBzt,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,PARzt,Wt]...
    = odlake_v0_1(Datetime,Initfile,Initsheet1,Inputfile,Inputsheet,Parafile,Parasheet);

% Experiment 2 (80% initial cyanobacteria) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,PO4zt,NO3zt,NH3zt,CyBzt2,APSBzt,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,PARzt,Wt]...
    = odlake_v0_1(Datetime,Initfile,Initsheet2,Inputfile,Inputsheet,Parafile,Parasheet);

% Experiment 3 (60% initial cyanobacteria) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,PO4zt,NO3zt,NH3zt,CyBzt3,APSBzt,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,PARzt,Wt]...
    = odlake_v0_1(Datetime,Initfile,Initsheet3,Inputfile,Inputsheet,Parafile,Parasheet);

% Experiment 4 (120% initial cyanobacteria) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,PO4zt,NO3zt,NH3zt,CyBzt4,APSBzt,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,PARzt,Wt]...
    = odlake_v0_1(Datetime,Initfile,Initsheet4,Inputfile,Inputsheet,Parafile,Parasheet);

% Experiment 5 (140% initial cyanobacteria) 

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,PO4zt,NO3zt,NH3zt,CyBzt5,APSBzt,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,PARzt,Wt]...
    = odlake_v0_1(Datetime,Initfile,Initsheet5,Inputfile,Inputsheet,Parafile,Parasheet);

%% Analyse the experiment results 

% Calculate the max of the cyanobacteria population in each day 

CyBt1 = max(CyBzt1); 

CyBt2 = max(CyBzt2); 

CyBt3 = max(CyBzt3); 

CyBt4 = max(CyBzt4); 

CyBt5 = max(CyBzt5); 

%% Plot 

% Fix time 

dt=0.1;

zlim = [0 max(zz)];

tlim = [0 length(tt)*dt];

tt_std = [0:dt:(length(tt)-1)*dt]';

% Make plot 

figure(1)

set(gcf,'Units','pixels','Position',[50 50 600 400])

plot(tt_std,CyBt1'/1000000,'k','LineWidth',1); hold on; 

plot(tt_std,CyBt2'/1000000,'b','LineWidth',1); hold on; 

plot(tt_std,CyBt3'/1000000,'m','LineWidth',1); hold on; 

plot(tt_std,CyBt4'/1000000,'g','LineWidth',1); hold on; 

plot(tt_std,CyBt5'/1000000,'r','LineWidth',1); hold on; 

set(gca,'XLim',[0 2500],'FontSize',7,'YLim',[40000 100000],'YTick',40000:10000:100000,...
    'YTickLabel',{'40000' '50000' '60000' '70000' '80000' '90000' '100000'});

xlabel('Time (day)','FontSize',9,'FontWeight','bold')

ylabel('Peak cyanobacteria population (cells mL^{-1})','FontSize',9,'FontWeight','bold')

legend({'Control' 'Cyanobacteria-' 'Cyanobacteria--' 'Cyanobacteria+' 'Cyanobacteria++'},'Location','northwest','Color','none','FontSize',6) 

exportgraphics(gcf,'Results\Test_1.jpg','Resolution',600)







