% This is to run the MyLake model using real data and plot the results 
clear all;clc;

%% Run the model for a specific day 

Datetime = "27-Jun-2015 00:00:00";
Initfile = "Init_data.xlsx";
Initsheet = 1;
Inputfile = "Input_data.xlsx";
Inputsheet = 1;
Parafile = "Para_data.xlsx";
Parasheet = 1;

[zz,Az,Vz,tt,Qst,Kzt,Tzt,rhozt,DOzt,pHzt,PO4zt,NO3zt,NH3zt,CyBzt,APSBzt,H2Szt,SO4zt,DFezt,PFezt,DMnzt,PMnzt,...
    Qzt_sed,PARzt,Wt,G_cybzt,G_apsbzt] = odlake_v6(Datetime,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet);

%% Write output file 
save 'OdLake_result_base.mat' zz Az Vz tt Qst Kzt Tzt rhozt DOzt pHzt PO4zt NO3zt NH3zt CyBzt APSBzt H2Szt SO4zt DFezt PFezt DMnzt PMnzt ...
    Qzt_sed PARzt Wt G_cybzt G_apsbzt;

%% Input raw dataset 
[Z,T,PAR,DO,pH,PO4,NO3,NH3,CyB,APSB,H2S,SO4,DFe,PFe,DMn,PMn]=rawdatainputs('LakeCadagnoData.xlsx',1);

%% Physical variables 

fontsize = 18;

axesPosition = [100 200 400 700]; % Axes position, in pixels 

yWidth = 40; % y axes spacing 

yLimit = [0 20];

yOffset = -yWidth*diff(yLimit)/axesPosition(3);

figure(1)

set(gcf,'Units','pixels','position',[100 100 1050 1000]) 

% plot observation 

% plot PAR

h1 = axes('Units','pixels','Position',axesPosition,'Color','w','XColor','g','YColor','k',...
    'YDir','reverse','XLim',[-3 3],'XTick',[-3 -2 -1 0 1 2 3],'XTickLabel',{'0.001' '0.01' '0.1' '1' '10' '100' '1000'},...
    'YLim',yLimit,'NextPlot','add');

plot(h1,log10(PAR(PAR>0)),Z(PAR>0),'-sg','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',2); hold on; 

set(h1,'FontSize',14,'XColor','k','YColor','k');

xlabel(h1,'PAR (µmol s^{-1} m^{-2})','FontSize',fontsize,'FontWeight','bold')

ylabel(h1,'Depth (m)','FontSize',fontsize,'FontWeight','bold')

% plot dissolved oxygen 

h2 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 500],'YLim',yLimit,'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h2,DO,Z,'-^b','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','b','Linewidth',2); hold on; 

set(h2,'FontSize',14,'XColor','k','YColor','k');

xlabel(h2,'Dissolved Oxygen (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

% plot temperature

h3 = axes('Units','pixels','Position',axesPosition+yWidth*[0 -2.5 0 2.5],'YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[4 14],'YLim',yLimit+[1.1*yOffset 0],'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h3,[-1 0],[-2 -1],'-sg','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',2)

plot(h3,[-1 0],[-2 -1],'-^b','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','b','Linewidth',2)

plot(h3,T,Z/1.03-2.2,'-or','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','r','Linewidth',2); hold on; 

set(h3,'FontSize',14,'XColor','k','YColor','none');

xlabel(h3,'Temperature (^{o}C)','FontSize',fontsize,'FontWeight','bold')

legend(h3,{'PAR','Dissolved Oxygen','Temperature'},'Location','northwest','Color','none','FontSize',12)

text(3.5,-4.7,'A','FontSize',fontsize,'FontWeight','bold')

h4 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','none','YColor','none','YLim',yLimit,'XTick',[],'XTickLabel',[],'YTickLabel',[],'NextPlot','add');

set(h4,'XColor','k','YColor','k')

box on;

% plot simulation 

axesPosition = [600 200 400 700]; 

% plot PAR 

h5 = axes('Units','pixels','Position',axesPosition,'Color','w','XColor','g','YColor','k',...
    'YDir','reverse','XLim',[-3 3],'XTick',[-3 -2 -1 0 1 2 3],'XTickLabel',{'0.001' '0.01' '0.1' '1' '10' '100' '1000'},...
    'YLim',yLimit,'NextPlot','add');

plot(h5,log10(PARzt(:,end)*1e6),zz,'g','Linewidth',2); hold on; 


set(h5,'FontSize',14,'XColor','k','YColor','k');

xlabel(h5,'PAR (µmol s^{-1} m^{-2})','FontSize',fontsize,'FontWeight','bold')

ylabel(h5,'Depth (m)','FontSize',fontsize,'FontWeight','bold')

% plot dissolved oxygen 

h6 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 500],'YLim',yLimit,'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h6,[-1 0],[-2 -1],'g','Linewidth',2)

plot(h6,DOzt(:,end)*1000,zz,'b','Linewidth',2); hold on; 

set(h6,'FontSize',14,'XColor','k','YColor','k');

xlabel(h6,'Dissolved Oxygen (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

legend(h6,{'PAR','Dissolved Oxygen'},'Location','northwest','Color','none','FontSize',12)

text(-30,-2.55,'B','FontSize',fontsize,'FontWeight','bold')

h7 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','none','YColor','none','YLim',yLimit,'XTick',[],'XTickLabel',[],'YTickLabel',[],'NextPlot','add');

set(h7,'XColor','k','YColor','k')

box on;

saveas(gcf,'Results\Phys_result.jpg')

%% Biological variables 

figure(2)

set(gcf,'Units','pixels','position',[100 200 1050 900]) 

% Observation 

subplot(121)

plot(CyB,Z,'-sg','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',2); hold on; 

plot(APSB,Z,'-vm','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','m','LineWidth',2)

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 600000],'XTick',[200000 400000 600000],'XTickLabel',{'200000' '400000' '600000'})

xlabel('Bacteria (cells mL^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

legend({'Cyanobacteria' 'APSB'},'Location','northeast','Color','none','FontSize',12) 

text(-20000,-1.5,'A','FontSize',fontsize,'FontWeight','bold')

% Simulation 

subplot(122)

plot(CyBzt(:,end)/1000000,zz,'g','LineWidth',2); hold on; 

plot(APSBzt(:,end)/1000000,zz,'m','LineWidth',2)

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 600000],'XTick',[200000 400000 600000],'XTickLabel',{'200000' '400000' '600000'})

xlabel('Bacteria (cells mL^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

legend({'Cyanobacteria' 'APSB'},'Location','northeast','Color','none','FontSize',12) 

text(-20000,-1.5,'B','FontSize',fontsize,'FontWeight','bold')

saveas(gcf,'Results\Bac_results.jpg')

%% Non-metal variables 

axesPosition = [100 200 400 700]; % Axes position, in pixels 

figure(3)

set(gcf,'Units','pixels','position',[100 100 1050 1000]) 

% Observation 

% plot NH3 & H2S 

h1 = axes('Units','pixels','Position',axesPosition,'Color','w','XColor','g','YColor','k',...
    'YDir','reverse','XLim',[0 150],'XTick',[0:30:150],'YLim',yLimit,'NextPlot','add');

plot(h1,NH3,Z,'-sm','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','m','Linewidth',2); hold on; 

plot(h1,H2S,Z,'-dc','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','c','LineWidth',2); hold on; 

box on;

set(h1,'FontSize',14,'XColor','k','YColor','k');

xlabel(h1,'NH_{3}, H_{2}S (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel(h1,'Depth (m)','FontSize',fontsize,'FontWeight','bold')

% plot NO3- 

h2 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 1000],'YLim',yLimit,'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h2,NO3,Z,'-^g','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',2); hold on; 

set(h2,'FontSize',14,'XColor','k','YColor','k');

xlabel(h2,'NO_{3}^{-} (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

% plot PO43- & SO42-

h3 = axes('Units','pixels','Position',axesPosition+yWidth*[0 -2.5 0 2.5],'YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 5],'YLim',yLimit+[1.1*yOffset 0],'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h3,[-1 0],[-2 -1],'-sm','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','m','Linewidth',2)

plot(h3,[-1 0],[-2 -1],'-dc','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','c','Linewidth',2)

plot(h3,[-1 0],[-2 -1],'-^g','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',2)

plot(h3,PO4,Z/1.03-2.2,'-or','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','r','Linewidth',2)

plot(h3,SO4,Z/1.03-2.2,'-vb','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','b','Linewidth',2)

set(h3,'FontSize',14,'XColor','k','YColor','none');

xlabel(h3,'PO_{4}^{3-} (µmol L^{-1}), SO_{4}^{2-} (mmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

legend(h3,{'NH_{3}','H_{2}S','NO_{3}^{-}','PO_{4}^{3-}','SO_{4}^{2-}'},'Location','northeast','Color','none','FontSize',12)

text(-0.1,-4.7,'A','FontSize',fontsize,'FontWeight','bold')

% plot simulation 

axesPosition = [600 200 400 700]; 

% plot NH3 & H2S 

h5 = axes('Units','pixels','Position',axesPosition,'Color','w','XColor','g','YColor','k',...
    'YDir','reverse','XLim',[0 150],'XTick',[0:30:150],'YLim',yLimit,'NextPlot','add');

plot(h5,NH3zt(:,end)*1000,zz,'-m','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','m','Linewidth',2); hold on; 

plot(h5,H2Szt(:,end)*1000,zz,'-c','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','c','LineWidth',2); hold on; 

box on;

set(h5,'FontSize',14,'XColor','k','YColor','k');

xlabel(h5,'NH_{3}, H_{2}S (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel(h5,'Depth (m)','FontSize',fontsize,'FontWeight','bold')

% plot NO3- 

h5 = axes('Units','pixels','Position',axesPosition,'XAxisLocation','top','YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 1000],'YLim',yLimit,'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h5,NO3zt(:,end)*1000,zz,'-g','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',2); hold on; 

set(h5,'FontSize',14,'XColor','k','YColor','k');

xlabel(h5,'NO_{3}^{-} (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

% plot PO43- & SO42-

h6 = axes('Units','pixels','Position',axesPosition+yWidth*[0 -2.5 0 2.5],'YDir','reverse','Color','none',...
    'XColor','b','YColor','none','XLim',[0 5],'YLim',yLimit+[1.1*yOffset 0],'YTick',[],'YTickLabel',[],'NextPlot','add');

plot(h6,[-1 -2],[-2 -1],'-m','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','m','Linewidth',2)

plot(h6,[-1 -2],[-2 -1],'-c','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','c','Linewidth',2)

plot(h6,[-1 -2],[-2 -1],'-g','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g','Linewidth',2)

plot(h6,PO4zt(:,end)*1000,zz/1.055-2.2,'-r','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','r','Linewidth',2)

plot(h6,SO4zt(:,end),zz/1.055-2.2,'-b','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','b','Linewidth',2)

set(h6,'FontSize',14,'XColor','k','YColor','none');

xlabel(h6,'PO_{4}^{3-} (µmol L^{-1}), SO_{4}^{2-} (mmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

legend(h6,{'NH_{3}','H_{2}S','NO_{3}^{-}','PO_{4}^{3-}','SO_{4}^{2-}'},'Location','northeast','Color','none','FontSize',12)

text(-0.1,-4.7,'B','FontSize',fontsize,'FontWeight','bold')

saveas(gcf,'Results\Nutn_result.jpg')


%% Iron

figure(4)

set(gcf,'Units','pixels','position',[100 200 1050 900]) 

% Observation 

subplot(121)

plot(DFe,Z,'-^g','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',2); hold on; 

plot(PFe,Z,'-sr','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',2)

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 1.5])

xlabel('Fe (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',12) 

text(-0.05,-1.5,'A','FontSize',fontsize,'FontWeight','bold')

% Simulation 

subplot(122)

plot(DFezt(:,end)*1000,zz,'g','LineWidth',2); hold on; 

plot(PFezt(:,end)*1000,zz,'r','LineWidth',2)

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 1.5])

xlabel('Fe (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',12) 

text(-0.05,-1.5,'B','FontSize',fontsize,'FontWeight','bold')

saveas(gcf,'Results\Fe_results.jpg')

%% Manganese 

figure(5)

set(gcf,'Units','pixels','position',[100 200 1050 900]) 

% Observation 

subplot(121)

plot(DMn,Z,'-^g','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',2); hold on; 

plot(PMn,Z,'-sr','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',2)

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 3])

xlabel('Mn (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',12) 

text(-0.1,-1.5,'A','FontSize',fontsize,'FontWeight','bold')

% Simulation 

subplot(122)

plot(DMnzt(:,end)*1000,zz,'g','LineWidth',2); hold on; 

plot(PMnzt(:,end)*1000,zz,'r','LineWidth',2)

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 3])

xlabel('Mn (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',12) 

text(-0.1,-1.5,'B','FontSize',fontsize,'FontWeight','bold')

saveas(gcf,'Results\Mn_results.jpg')

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

figure(6)

set(gcf,'Units','pixels','position',[100 200 1250 900])

% Cyanobacteria 

subplot(121)

plot(CyBzt(z_r,end)/1e6-CyB,Z,'-g','LineWidth',2)

xline(0,'--')

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[-50000 50000],...
    'XTick',[-50000 -25000 0 25000 50000],'XTickLabel',{'-50000' '-25000' '0' '25000' '50000'})

xlabel('Residual (cells mL^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

legend({'Cyanobacteria'},'Location','northeast','Color','none','FontSize',12) 

text(-52000,-1.5,'A','FontSize',fontsize,'FontWeight','bold')

subplot(122)

plot(APSBzt(z_r,end)/1e6-APSB,Z,'-m','LineWidth',2)

xline(0,'--')

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[-200000 200000],...
    'XTick',[-200000 -100000 0 100000 200000],'XTickLabel',{'-200000' '-100000' '0' '100000' '200000'})

xlabel('Residual (cells mL^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

legend({'APSB'},'Location','northeast','Color','none','FontSize',12)

text(-210000,-1.5,'B','FontSize',fontsize,'FontWeight','bold')

saveas(gcf,'Results\Residual_bacteria.jpg')

%% Iron and manganese discussion plot 

figure(7)

set(gcf,'Units','pixels','position',[100 0 1050 1400]) 

% Observation (Fe)

subplot(221)

plot(DFe,Z,'-^g','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',2); hold on; 

plot(PFe,Z,'-sr','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',2)

yline([9.5 12 13.5],'--')

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 1.5])

xlabel('Fe (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',12) 

text(-0.05,-3,'A','FontSize',fontsize,'FontWeight','bold')

% Observation (Mn) 

subplot(222)

plot(DMn,Z,'-^g','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',2); hold on; 

plot(PMn,Z,'-sr','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',2)

yline([9.5 12 13.5],'--')

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 3])

xlabel('Mn (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

text(-0.1,-3,'B','FontSize',fontsize,'FontWeight','bold')

% Observation 

subplot(223)

plot(DFezt(:,end)*1000,zz,'g','LineWidth',2); hold on; 

plot(PFezt(:,end)*1000,zz,'r','LineWidth',2) 

yline([10 12],'--')

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 1.5])

xlabel('Fe (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

legend({'Dissolved' 'Particulate'},'Location','northeast','Color','none','FontSize',12) 

text(-0.05,-3,'C','FontSize',fontsize,'FontWeight','bold')

% Simulation 

subplot(224)

plot(DMnzt(:,end)*1000,zz,'g','LineWidth',2); hold on; 

plot(PMnzt(:,end)*1000,zz,'r','LineWidth',2)

yline([10 12],'--')

set(gca,'FontSize',14,'XAxisLocation','top','YDir','reverse','YLim',yLimit,'XLim',[0 3])

xlabel('Mn (µmol L^{-1})','FontSize',fontsize,'FontWeight','bold')

ylabel('Depth (m)','FontSize',fontsize,'FontWeight','bold')

text(-0.1,-3,'D','FontSize',fontsize,'FontWeight','bold')

saveas(gcf,'Results\Fe_Mn_discussion.jpg')
