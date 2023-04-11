% This is to run the OdLake model using real data and plot the results 
% Plots are used in JGR: Biogeosciences
clear all;clc;

%% Run the model for a specific day 

Datetime = "27-Jun-2015 00:00:00";
Initfile = "Init_data.xlsx";
Initsheet = 1;
Inputfile = "Input_data.xlsx";
Inputsheet = 1;
Parafile = "Para_data.xlsx";
Parasheet = 1;

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,PO4zt,NO3zt,NH3zt,CyBzt,APSBzt,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,PARzt,Wt]...
    = odlake_v0_1(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet);

%% Write output file 
save 'OdLake_result_base.mat' zz Az Vz tt Qst Kzt Tzt rhozt DOzt PO4zt NO3zt NH3zt CyBzt APSBzt H2Szt SO4zt DFezt PFezt DMnzt PMnzt PARzt Wt;

%% Input raw dataset 
[Z,T,PAR,DO,pH,PO4,NO3,NH3,CyB,APSB,H2S,SO4,DFe,PFe,DMn,PMn]=rawdatainputs('LakeCadagnoData.xlsx',1);

%% Result figure 1 

fontsize1 = 9; % General font size 

fontsize2 = 7; % x/y ticks font size 

fontsize3 = 6; % Legend font size 

markersize = 6; 

linewidth = 1; 

yWidth = 20; % y axes spacing 

yLimit = [0 20];

axesPosition = [50 950 200 300]; % Axes position, in pixels 

yOffset = -yWidth*diff(yLimit)/axesPosition(3);

figure(1)

set(gcf,'Units','pixels','position',[50 50 1025 1300]) 

% Physical variables 

% plot observation 

% plot PAR

axesPosition = [50 950 200 300]; % Axes position, in pixels 

h1 = axes('Units','pixels','Position',axesPosition,'Color','w','XColor','g','YColor','k',...
    'YDir','reverse','XLim',[-3 3],'XTick',[-3 -2 -1 0 1 2 3],'XTickLabel',{'0.001' '0.01' '0.1' '1' '10' '100' '1000'},...
    'YLim',yLimit,'NextPlot','add');

plot(h1,log10(PAR(PAR>0)),Z(PAR>0),'-sg','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',linewidth); hold on; 

set(h1,'FontSize',fontsize2,'XColor','k','YColor','k');

xlabel(h1,'PAR (µmol s^{-1} m^{-2})','FontSize',fontsize1,'FontWeight','bold')

ylabel(h1,'Depth (m)','FontSize',fontsize1,'FontWeight','bold')

% plot dissolved oxygen 

h2 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 500],'YLim',yLimit,'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h2,DO,Z,'-^b','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','b','Linewidth',linewidth); hold on; 

set(h2,'FontSize',fontsize2,'XColor','k','YColor','k');

xlabel(h2,'Dissolved Oxygen (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

% plot temperature

h3 = axes('Units','pixels','Position',axesPosition+yWidth*[0 -2.5 0 2.5],'YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[4 14],'YLim',yLimit+[1.1*yOffset 0],'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h3,[-1 0],[-2 -1],'-sg','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',linewidth)

plot(h3,[-1 0],[-2 -1],'-^b','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','b','Linewidth',linewidth)

plot(h3,T,Z/1.03-2.2,'-or','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','r','Linewidth',linewidth); hold on; 

set(h3,'FontSize',fontsize2,'XColor','k','YColor','none');

xlabel(h3,'Temperature (^{o}C)','FontSize',fontsize1,'FontWeight','bold')

lgd = legend(h3,{'PAR','Dissolved Oxygen','Temperature'},'Location','northwest','Color','none','FontSize',fontsize3);

title(lgd,'Observation')

text(3.5,-4.7,'a','FontSize',fontsize1,'FontWeight','bold')

h4 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','none','YColor','none','YLim',yLimit,'XTick',[],'XTickLabel',[],'YTickLabel',[],'NextPlot','add');

set(h4,'XColor','k','YColor','k')

box on;

% plot simulation 

axesPosition = [300 950 200 300]; % Axes position, in pixels 

% plot PAR 

h5 = axes('Units','pixels','Position',axesPosition,'Color','w','XColor','g','YColor','k',...
    'YDir','reverse','XLim',[-3 3],'XTick',[-3 -2 -1 0 1 2 3],'XTickLabel',{'0.001' '0.01' '0.1' '1' '10' '100' '1000'},...
    'YLim',yLimit,'NextPlot','add');

plot(h5,log10(PARzt(:,end)*1e6),zz,'g','Linewidth',linewidth); hold on; 

set(h5,'FontSize',fontsize2,'XColor','k','YColor','k');

xlabel(h5,'PAR (µmol s^{-1} m^{-2})','FontSize',fontsize1,'FontWeight','bold')

ylabel(h5,'Depth (m)','FontSize',fontsize1,'FontWeight','bold')

% plot dissolved oxygen 

h6 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 500],'YLim',yLimit,'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h6,[-1 0],[-2 -1],'g','Linewidth',linewidth)

plot(h6,DOzt(:,end)*1000,zz,'b','Linewidth',linewidth); hold on; 

set(h6,'FontSize',fontsize2,'XColor','k','YColor','k');

xlabel(h6,'Dissolved Oxygen (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

lgd = legend(h6,{'PAR','Dissolved Oxygen'},'Location','northwest','Color','none','FontSize',fontsize3);

title(lgd,'Simulation')

text(-30,-2.55,'b','FontSize',fontsize1,'FontWeight','bold')

h7 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','none','YColor','none','YLim',yLimit,'XTick',[],'XTickLabel',[],'YTickLabel',[],'NextPlot','add');

set(h7,'XColor','k','YColor','k')

box on;

% Biological variables 

% Observation 

axesPosition = [550 950 200 300]; % Axes position, in pixels 

h1 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse',...
    'YLim',yLimit,'XLim',[0 600000],'XTick',[200000 400000 600000],'XTickLabel',{'200000' '400000' '600000'},'NextPlot','add');

plot(h1,CyB,Z,'-sg','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',linewidth); hold on; 

plot(h1,APSB,Z,'-vm','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','m','LineWidth',linewidth)

set(h1,'FontSize',fontsize2,'XColor','k','YColor','k');

xlabel('Bacteria (cells mL^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

lgd = legend({'Cyanobacteria' 'APSB'},'Location','northeast','Color','none','FontSize',fontsize3);

title(lgd,'Observation')

text(-20000,-2.55,'c','FontSize',fontsize1,'FontWeight','bold')

box on;

% Simulation 

axesPosition = [800 950 200 300]; % Axes position, in pixels 

h2 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse',...
    'YLim',yLimit,'XLim',[0 600000],'XTick',[200000 400000 600000],'XTickLabel',{'200000' '400000' '600000'},'NextPlot','add');

plot(h2,CyBzt(:,end)/1000000,zz,'g','LineWidth',linewidth); hold on; 

plot(h2,APSBzt(:,end)/1000000,zz,'m','LineWidth',linewidth)

set(h2,'FontSize',fontsize2,'XColor','k','YColor','k');

xlabel('Bacteria (cells mL^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

lgd = legend({'Cyanobacteria' 'APSB'},'Location','northeast','Color','none','FontSize',fontsize3);

title(lgd,'Simulation')

text(-20000,-2.55,'d','FontSize',fontsize1,'FontWeight','bold')

box on;

% Non-metal variables 

% Observation 

% plot NH3 & H2S 

axesPosition = [50 500 200 300]; % Axes position, in pixels 

h1 = axes('Units','pixels','Position',axesPosition,'Color','w','XColor','g','YColor','k',...
    'YDir','reverse','XLim',[0 150],'XTick',[0:30:150],'YLim',yLimit,'NextPlot','add');

plot(h1,NH3,Z,'-sm','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','m','Linewidth',linewidth); hold on; 

plot(h1,H2S,Z,'-dc','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','c','LineWidth',linewidth); hold on; 

box on;

set(h1,'FontSize',fontsize2,'XColor','k','YColor','k');

xlabel(h1,'NH_{3}, H_{2}S (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel(h1,'Depth (m)','FontSize',fontsize1,'FontWeight','bold')

% plot NO3- 

h2 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 1000],'YLim',yLimit,'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h2,NO3,Z,'-^g','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',linewidth); hold on; 

set(h2,'FontSize',fontsize2,'XColor','k','YColor','k');

xlabel(h2,'NO_{3}^{-} (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

% plot PO43- & SO42-

h3 = axes('Units','pixels','Position',axesPosition+yWidth*[0 -2.5 0 2.5],'YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 5],'YLim',yLimit+[1.1*yOffset 0],'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h3,[-1 0],[-2 -1],'-sm','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','m','Linewidth',linewidth)

plot(h3,[-1 0],[-2 -1],'-dc','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','c','Linewidth',linewidth)

plot(h3,[-1 0],[-2 -1],'-^g','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',linewidth)

plot(h3,PO4,Z/1.03-2.2,'-or','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','r','Linewidth',linewidth)

plot(h3,SO4,Z/1.03-2.2,'-vb','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','b','Linewidth',linewidth)

set(h3,'FontSize',fontsize2,'XColor','k','YColor','none');

xlabel(h3,'PO_{4}^{3-} (µmol L^{-1}), SO_{4}^{2-} (mmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

lgd = legend(h3,{'NH_{3}','H_{2}S','NO_{3}^{-}','PO_{4}^{3-}','SO_{4}^{2-}'},'Location','northeast','Color','none','FontSize',fontsize3);

title(lgd,'Observation')

text(-0.1,-4.7,'e','FontSize',fontsize1,'FontWeight','bold')

% plot simulation 

axesPosition = [300 500 200 300]; % Axes position, in pixels 

% plot NH3 & H2S 

h5 = axes('Units','pixels','Position',axesPosition,'Color','w','XColor','g','YColor','k',...
    'YDir','reverse','XLim',[0 150],'XTick',[0:30:150],'YLim',yLimit,'NextPlot','add');

plot(h5,NH3zt(:,end)*1000,zz,'-m','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','m','Linewidth',linewidth); hold on; 

plot(h5,H2Szt(:,end)*1000,zz,'-c','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','c','LineWidth',linewidth); hold on; 

box on;

set(h5,'FontSize',fontsize2,'XColor','k','YColor','k');

xlabel(h5,'NH_{3}, H_{2}S (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel(h5,'Depth (m)','FontSize',fontsize1,'FontWeight','bold')

% plot NO3- 

h5 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 1000],'YLim',yLimit,'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h5,NO3zt(:,end)*1000,zz,'-g','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',linewidth); hold on; 

set(h5,'FontSize',fontsize2,'XColor','k','YColor','k');

xlabel(h5,'NO_{3}^{-} (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

% plot PO43- & SO42-

h6 = axes('Units','pixels','Position',axesPosition+yWidth*[0 -2.5 0 2.5],'YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 5],'YLim',yLimit+[1.1*yOffset 0],'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h6,[-1 -2],[-2 -1],'-m','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','m','Linewidth',linewidth)

plot(h6,[-1 -2],[-2 -1],'-c','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','c','Linewidth',linewidth)

plot(h6,[-1 -2],[-2 -1],'-g','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',linewidth)

plot(h6,PO4zt(:,end)*1000,zz/1.08-2.2,'-r','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','r','Linewidth',linewidth)

plot(h6,SO4zt(:,end),zz/1.08-2.2,'-b','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','b','Linewidth',linewidth)

set(h6,'FontSize',fontsize2,'XColor','k','YColor','none');

xlabel(h6,'PO_{4}^{3-} (µmol L^{-1}), SO_{4}^{2-} (mmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

lgd = legend(h6,{'NH_{3}','H_{2}S','NO_{3}^{-}','PO_{4}^{3-}','SO_{4}^{2-}'},'Location','northeast','Color','none','FontSize',fontsize3);

title(lgd,'Simulation')

text(-0.1,-4.7,'f','FontSize',fontsize1,'FontWeight','bold')

% Iron

% Observation 

axesPosition = [550 500 200 300]; % Axes position, in pixels 

h1 = axes('Units','pixels','Position',axesPosition,'FontSize',fontsize2,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 1.5],'NextPlot','add');

plot(h1,DFe,Z,'-^g','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',linewidth); hold on; 

plot(h1,PFe,Z,'-sr','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',linewidth)

xlabel('Fe (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

lgd = legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',fontsize3);

title(lgd,'Observation')

text(-0.05,-2.55,'g','FontSize',fontsize1,'FontWeight','bold')

box on;

% Simulation 

axesPosition = [800 500 200 300]; % Axes position, in pixels 

h2 = axes('Units','pixels','Position',axesPosition,'FontSize',fontsize2,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 1.5],'NextPlot','add');

plot(h2,DFezt(:,end)*1000,zz,'g','LineWidth',linewidth); hold on; 

plot(h2,PFezt(:,end)*1000,zz,'r','LineWidth',linewidth)

xlabel('Fe (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

lgd = legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',fontsize3);

title(lgd,'Simulation')

text(-0.05,-2.55,'h','FontSize',fontsize1,'FontWeight','bold')

box on;

% Manganese 

% Observation 

axesPosition = [50 50 200 300]; % Axes position, in pixels 

h1 = axes('Units','pixels','Position',axesPosition,'FontSize',fontsize2,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 3],'NextPlot','add');

plot(h1,DMn,Z,'-^g','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',linewidth); hold on; 

plot(h1,PMn,Z,'-sr','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',linewidth)

xlabel('Mn (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

lgd = legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',fontsize3);

title(lgd,'Observation')

text(-0.1,-2.55,'i','FontSize',fontsize1,'FontWeight','bold')

box on;

% Simulation 

axesPosition = [300 50 200 300]; % Axes position, in pixels 

h2 = axes('Units','pixels','Position',axesPosition,'FontSize',fontsize2,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 3],'NextPlot','add');

plot(h2,DMnzt(:,end)*1000,zz,'g','LineWidth',linewidth); hold on; 

plot(h2,PMnzt(:,end)*1000,zz,'r','LineWidth',linewidth)

xlabel('Mn (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

lgd = legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',fontsize3);

title(lgd,'Simulation')

text(-0.1,-2.55,'j','FontSize',fontsize1,'FontWeight','bold')

box on;

exportgraphics(gcf,'Results\Result_figure1.jpg','Resolution',600)

%% Residual (bacteria) 

% Filter the simulated data in the same depth with measured samples 

k=1;

for i=1:length(zz)

    for j=1:length(Z)

        if zz(i)==Z(j)

            z_r(k,1)=i;

            k=k+1;

        end

    end

end

figure(2)

set(gcf,'Units','pixels','position',[50 50 525 400])

% Cyanobacteria 

axesPosition = [50 50 200 300]; % Axes position, in pixels 

h1 = axes('Units','pixels','Position',axesPosition,'FontSize',fontsize2,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[-50000 50000],...
    'XTick',[-50000 -25000 0 25000 50000],'XTickLabel',{'-50000' '-25000' '0' '25000' '50000'},'NextPlot','add');

plot(h1,CyBzt(z_r,end)/1e6-CyB,Z,'-g','LineWidth',linewidth)

xline(0,'--')

xlabel('Residual (cells mL^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

legend({'Cyanobacteria'},'Location','northeast','Color','none','FontSize',fontsize3) 

text(-52000,-2.55,'a','FontSize',fontsize1,'FontWeight','bold')

box on;

% APSB 

axesPosition = [300 50 200 300]; % Axes position, in pixels 

h2 = axes('Units','pixels','Position',axesPosition,'FontSize',fontsize2,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[-250000 250000],...
    'XTick',[-250000 -125000 0 125000 250000],'XTickLabel',{'-250000' '-125000' '0' '125000' '250000'},'NextPlot','add');

plot(h2,APSBzt(z_r,end)/1e6-APSB,Z,'-m','LineWidth',linewidth)

xline(0,'--')

xlabel('Residual (cells mL^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

legend({'APSB'},'Location','northeast','Color','none','FontSize',fontsize3)

text(-210000,-2.55,'b','FontSize',fontsize1,'FontWeight','bold')

box on;

exportgraphics(gcf,'Results\Residual_bacteria.jpg','Resolution',600)

%% Iron and manganese discussion plot 

figure(3)

set(gcf,'Units','pixels','position',[50 50 525 750]) 

% Observation (Fe)

axesPosition = [50 400 200 300]; % Axes position, in pixels 

h1 = axes('Units','pixels','Position',axesPosition,'FontSize',fontsize2,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 1.5],'NextPlot','add');

plot(DFe,Z,'-^g','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',linewidth); hold on; 

plot(PFe,Z,'-sr','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',linewidth)

yline([9.5 12 13.5],'--')

xlabel('Fe (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',fontsize3) 

text(-0.05,-2.55,'a','FontSize',fontsize1,'FontWeight','bold')

box on;

% Observation (Mn) 

axesPosition = [300 400 200 300]; % Axes position, in pixels 

h2 = axes('Units','pixels','Position',axesPosition,'FontSize',fontsize2,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 3],'NextPlot','add');

plot(DMn,Z,'-^g','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',linewidth); hold on; 

plot(PMn,Z,'-sr','MarkerSize',markersize,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',linewidth)

yline([9.5 12 13.5],'--')

xlabel('Mn (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

text(-0.1,-2.55,'b','FontSize',fontsize1,'FontWeight','bold')

box on;

% Observation 

axesPosition = [50 50 200 300]; % Axes position, in pixels 

h3 = axes('Units','pixels','Position',axesPosition,'FontSize',fontsize2,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 1.5],'NextPlot','add');

plot(DFezt(:,end)*1000,zz,'g','LineWidth',linewidth); hold on; 

plot(PFezt(:,end)*1000,zz,'r','LineWidth',linewidth) 

yline([10 12],'--')

xlabel('Fe (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',fontsize3) 

text(-0.05,-2.55,'c','FontSize',fontsize1,'FontWeight','bold')

box on;

% Simulation 

axesPosition = [300 50 200 300]; % Axes position, in pixels 

h4 = axes('Units','pixels','Position',axesPosition,'FontSize',fontsize2,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 3],'NextPlot','add');

plot(DMnzt(:,end)*1000,zz,'g','LineWidth',linewidth); hold on; 

plot(PMnzt(:,end)*1000,zz,'r','LineWidth',linewidth)

yline([10 12],'--')

xlabel('Mn (µmol L^{-1})','FontSize',fontsize1,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize1,'FontWeight','bold')

text(-0.1,-2.55,'d','FontSize',fontsize1,'FontWeight','bold')

box on;

exportgraphics(gcf,'Results\Fe_Mn_discussion.jpg','Resolution',600)
