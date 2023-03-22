% This is to run the OdLake model v0.1 using real data and plot the results 
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

%% Input raw dataset 
[Z,T,PAR,DO,pH,PO4,NO3,NH3,CyB,APSB,H2S,SO4,DFe,PFe,DMn,PMn]=rawdatainputs('LakeCadagnoData.xlsx',1);

%% Physical components profiles 
figure(1)
subplot(1,3,1)
plot(T,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('Temperature')
ylabel('Depth (m)')
ylim([0 20])

subplot(1,3,2)
plot(pH,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('pH')
ylabel('Depth (m)')
ylim([0 20])

subplot(1,3,3)
plot(log10(PARzt(:,2)*1000000),zz,'b',log10(PAR(PAR>0)),Z(PAR>0),'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('PAR (µmol s-1 m-2)')
ylabel('Depth (m)')
ylim([0 20])
xticks([-3 -2 -1 0 1 2 3])
xticklabels({'0.001' '0.01' '0.1' '1' '10' '100' '1000'})
legend('Simulation', 'Real data')


%% Biological components 
figure(2)
subplot(1,3,1)
plot(DOzt(:,end)*1000,zz,'b',DO,Z,'r')
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('Dissolved Oxygen (µmol/L)')
ylabel('Depth (m)')
ylim([0 20])
legend('Simulation', 'Real data')

subplot(1,3,2)
plot(CyBzt(:,end)/1000000,zz,'b',CyB,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('Cyanobacteria (cells/mL)')
ylabel('Depth (m)')
ylim([0 20])

subplot(1,3,3)
plot(APSBzt(:,end)/1000000,zz,'b',APSB,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('APSB (cells/mL)')
ylabel('Depth (m)')
ylim([0 20])

%% Nutrient profiles 
figure(3)
subplot(1,5,1)
plot(PO4zt(:,end)*1000,zz,'b',PO4,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('PO4 (µmol/L)')
ylabel('Depth (m)')
ylim([0 20])
legend('Simulation', 'Real data')

subplot(1,5,2)
plot(NO3zt(:,end)*1000,zz,'b',NO3,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('NO3 (µmol/L)')
ylabel('Depth (m)')
ylim([0 20])

subplot(1,5,3)
plot(NH3zt(:,end)*1000,zz,'b',NH3,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('NH4 (µmol/L)')
ylabel('Depth (m)')
ylim([0 20])

subplot(1,5,4)
plot(H2Szt(:,end)*1000,zz,'b',H2S,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('H2S (µmol/L)')
ylabel('Depth (m)')
ylim([0 20])

subplot(1,5,5)
plot(SO4zt(:,end),zz,'b',SO4,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('SO4 (mmol/L)')
ylabel('Depth (m)')
ylim([0 20])

%% Metal 
figure(4)
subplot(1,4,1)
plot(DFezt(:,end)*1000,zz,'b',DFe,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('DFe (µmol/L)')
ylabel('Depth (m)')
ylim([0 20])
legend('Simulation', 'Real data')

subplot(1,4,2)
plot(PFezt(:,end)*1000,zz,'b',PFe,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('PFe (µmol/L)')
ylabel('Depth (m)')
ylim([0 20])

subplot(1,4,3)
plot(DMnzt(:,end)*1000,zz,'b',DMn,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('DMn (µmol/L)')
ylabel('Depth (m)')
ylim([0 20])

subplot(1,4,4)
plot(PMnzt(:,end)*1000,zz,'b',PMn,Z,'r');
ax=gca;
ax.YDir='reverse';
ax.XAxisLocation='top';
xlabel('PMn (µmol/L)')
ylabel('Depth (m)')
ylim([0 20])

%% Routine color map
dt=0.1;
zlim = [0 max(zz)];
tlim = [0 length(tt)*dt];
tt_std = [0:dt:(length(tt)-1)*dt]';

% thermocline depth
figure(5)
subplot(131)
contourf(tt_std,zz,DOzt*1000,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 400]);
colorbar;
ylabel(colorbar,'Dissolved Oxygen (µmol/L)')
set(gca,'fontsize',9);
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out')

subplot(132)
contourf(tt_std,zz,CyBzt/1e6,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 10e4]);
colorbar;
ylabel(colorbar,'Cyanobacteria (cells/mL)')
set(gca,'fontsize',9);
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out')

subplot(133)
contourf(tt_std,zz,APSBzt/1e6,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 5e5]);
colorbar;
ylabel(colorbar,'APSB (cells/mL)')
set(gca,'fontsize',9);
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out') 

figure(6)
subplot(151)
contourf(tt_std,zz,PO4zt*1000,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 4]);
colorbar;
ylabel(colorbar,'PO4 (µmol/L)')
set(gca,'fontsize',9); 
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out')

subplot(152)
contourf(tt_std,zz,NO3zt*1000,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 1000]);
colorbar;
ylabel(colorbar,'NO3 (µmol/L)')
set(gca,'fontsize',9);
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out')

subplot(153)
contourf(tt_std,zz,NH3zt*1000,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 60]);
colorbar;
ylabel(colorbar,'NH3 (µmol/L)')
set(gca,'fontsize',9);
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out')

subplot(154)
contourf(tt_std,zz,H2Szt*1000,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 160]);
colorbar;
ylabel(colorbar,'H2S (µmol/L)')
set(gca,'fontsize',9);
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out')

subplot(155)
contourf(tt_std,zz,SO4zt,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 5]);
colorbar;
ylabel(colorbar,'SO4 (mmol/L)')
set(gca,'fontsize',9);
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out')

figure(7)
subplot(141)
contourf(tt_std,zz,DFezt*1000,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 1.3]);
colorbar;
ylabel(colorbar,'DFe (µmol/L)')
set(gca,'fontsize',9);
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out')

subplot(142)
contourf(tt_std,zz,PFezt*1000,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 0.5]);
colorbar;
ylabel(colorbar,'PFe (µmol/L)')
set(gca,'fontsize',9);
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out')

subplot(143)
contourf(tt_std,zz,DMnzt*1000,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 2.2]);
colorbar;
ylabel(colorbar,'DMn (µmol/L)')
set(gca,'fontsize',9);
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out')

subplot(144)
contourf(tt_std,zz,PMnzt*1000,15,'linestyle','none')
shading interp
axis ij
set(gca,'xlim',tlim);
set(gca,'ylim',zlim);
caxis([0 1.5]);
colorbar;
ylabel(colorbar,'PMn (µmol/L)')
set(gca,'fontsize',9);
xlabel('Days')
ylabel('Depth (m)')
set(gca,'TickDir','out')

