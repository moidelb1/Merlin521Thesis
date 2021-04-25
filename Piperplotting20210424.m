%Benjamin Moidel
%April 24, 2021
%Pilot Test Result Plotting - Piper Arrow III

clc;clear;close all

%Load pilot test data
load Piperdata

%% Setup
%Pilot Selection
pilots = 6:11;      %only v10 pilots

%"Ideal" performance parameters: [min max]
cr = [700 1200];    %climb rate, clean [ft/min]
crpoh = 840;        %climb rate, clean from POH [ft/min]
sr = [400 600];     %sink rate, clean [ft/min]
sr1 = [400 600];    %sink rate, flaps 1 [ft/min]
sr2 = [400 600];    %sink rate, flaps 2 [ft/min]
sr3 = [400 600];    %sink rate, flaps 3 [ft/min]
Vapp = [75 90];     %approach speed [KIAS]
Vappf = [75 75];    %approach speed w/ flaps [KIAS]
Vclimb = [90 90];   %climb speed [KIAS]
Vs = 60;            %stall speed, clean [kn]
Vs0 = 55;           %stall speed, flaps [kn]
Vfe = 103;           %maximum flaps extended speed [kn]
Va = 118;           %maximum maneuvering speed, max weight [kn]
Vy = 90;            %best ROC speed, sea level [kn]
% Vto = 75;           %takeoff speed [kn]
Vto = 71;           %takeoff speed from POH [kn]
% Vtosf = 60;         %short field takeoff speed [kn]
Vtosf = 59;         %short field takeoff speed from POH [kn]
dto = 1800;         %takeoff distance [ft]
dtosf = 1000;       %short field takeoff distance [ft]
V3 = 115;           %speed for cond 3 S-Turn [kn]
V4 = 110;           %speed for cond 4 360 turn [kn]
V8 = 90;            %speed for cond 8 S-Turn w/ flaps [kn]
whead4 = 235;       %wind heading for cond 4 [degrees]
wspeed4 = 23;       %wind speed for cond 4 [kn]
whead5 = 130;       %wind heading for cond 5 [degrees]
wspeed5 = 23;       %wind speed for cond 5 [degrees]

%Legend switch
showlgd = 0;    %1 = show legend, 0 = no legend

%Get info about data
[np,nc,nt] = size(s);   %number of pilots, conditions, and trials

cp = interp1(1:64,hsv(64),linspace(1,60,length(pilots)));   %assign colors
                                                            %to pilots
ct = linspace(0.25,1,nt);       %brightness multiplier for trial #

lgdcell{length(pilots)*nt} = [];    %preallocate
for k = 1:nt    %build cell array for legend entries
    for i = pilots
        lgdcell{i-pilots(1)+1+(k-1)*length(pilots)} = sprintf('C%02i Trial %i',i,k);
    end
end

%Condition Names
condstr = {'Takeoff & Climb';
           'Descent & Landing';
           'S-Turn';
           '360 Turn (About a Point)';
           'Crosswind Landing';
           'Short Field Takeoff';
           'Landing w/ Flaps';
           'S-Turn w/ Flaps'};

%% Condition #1: Takeoff & Climb
close all
cond = 1;   %condition number

figure
hold on
ax1 = gca;  %get axes for altitude plot
figure
hold on
ax2 = gca;  %get axes for vertical speed plot (climb)
figure
hold on
ax3 = gca;  %get axes for elevator plot

figure
hold on
ax4 = gca;  %get axes for takeoff distance plot

figure
hold on
ax5 = gca;  %get axes for vertical speed plot (takeoff)

figure
hold on
ax6 = gca;  %get axes for altitude plot (takeoff)

figure
hold on
ax7 = gca;  %get axes for vertical force plot (takeoff)

figure
hold on
ax8 = gca;  %get axes for speed plot (takeoff)

tclimb = NaN(length(pilots),nt);     %preallocate
VTO = NaN(length(pilots),nt);
dTO = NaN(length(pilots),nt);
for k = 1:nt
    for i = pilots
        if isempty(s(i,cond,k).time)    %check for missing data
            %Plot dummy point
            plot(0,0,'k--')
        else
            alt = s(i,cond,k).alt/.3048;    %altitude [ft]
            t = s(i,cond,k).time;   %time [seconds]
            hdot = s(i,cond,k).hdot/.3048*60;   %climb rate [ft/min]
            de = s(i,cond,k).eta;   %elevator deflection [-1 to 1]
            N = s(i,cond,k).north/.3048;    %north position [ft]
            E = s(i,cond,k).east/.3048;     %east position [ft]
            V = s(i,cond,k).VIAS;       %indicated airspeed [kn]
            Z = s(i,cond,k).Z;          %vertical force [N]
            
%             indTO = find(alt>5,1); %first altitude over 10ft, broke ground
%             indTO = find(hdot>30,1); %broke ground when hdot > 30 fpm
            indTO = find(alt<4,1,'last');
            tTO = t(1:ceil(indTO*1.2));
            hdotTO = hdot(1:ceil(indTO*1.2));
            altTO = alt(1:ceil(indTO*1.2));
            ZTO = Z(1:ceil(indTO*1.2));
%             NTO = N(1:indTO);
%             ETO = E(1:indTO);
%             VTO = V(1:indTO);
            
            inds = and(alt<1350,alt>150);   %indices of interest to plot
            alt = alt(inds);
            t = t(inds);
            hdot = hdot(inds);
            de = de(inds);
            
            %Align data with origin of plot
            t = t-t(1);
            
            %Calculating performance metrics
            j = i-pilots(1)+1;      %alternate index variable
            tclimb(j,k) = t(end);   %saving last time for climb rate calcs
            VTO(j,k) = V(indTO);    %takeoff speed
            dTO(j,k) = norm([N(indTO),E(indTO)]);   %takeoff distance
            
            plot(ax1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax2,t,hdot,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax3,t,de,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax4,dTO(j,k),VTO(j,k),'*','Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax5,tTO,hdotTO,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax5,tTO(indTO),hdotTO(indTO),'*','Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax6,tTO,altTO,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax6,tTO(indTO),altTO(indTO),'*','Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax7,tTO,ZTO,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax7,tTO(indTO),ZTO(indTO),'*','Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax8,V(1:ceil(indTO*1.2)),ZTO,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax8,VTO(j,k),ZTO(indTO),'*','Color',cp(i-pilots(1)+1,:)*ct(k))
        end
    end
end
t1 = linspace(0,1200/cr(2)*60,100); %time vector [seconds]
t2 = linspace(0,1200/cr(1)*60,100); %time vector [seconds]
t3 = linspace(0,1200/crpoh*60,100); %time vector [seconds]
altmax = t1/60*cr(2)+150;           %altitude at maximum climb rate [ft]
altmin = t2/60*cr(1)+150;           %altitude at minimum climb rate [ft]
altpoh = t3/60*crpoh+150;           %altitude at POH climb rate [ft]

plot(ax1,t1,altmax,'k--','LineWidth',2)
plot(ax1,t2,altmin,'k--','LineWidth',2)
plot(ax1,t3,altpoh,'b--','LineWidth',2)
plot(ax2,t1,ones(size(t1))*cr(2),'k--','LineWidth',2)
plot(ax2,t2,ones(size(t2))*cr(1),'k--','LineWidth',2)
plot(ax2,t3,ones(size(t3))*crpoh,'k--','LineWidth',2)
plot(ax4,dto,Vto,'ko','MarkerSize',10,'LineWidth',2)

cravg = mean((1350-150)./tclimb*60,'all','omitnan'); %average climb rate
fprintf('Average climb rate for Condition #%i:\t%.1f fpm\n',cond,cravg)

%Plot Formatting
titstr = sprintf('Condition #%i: %s',cond,condstr{cond});

grid(ax1,'on')
xlabel(ax1,'Time (seconds)')
ylabel(ax1,'Altitude (ft)')
if showlgd
    legend(ax1,lgdcell,'NumColumns',nt,'Location','southeast')
end
title(ax1,titstr)

grid(ax2,'on')
xlabel(ax2,'Time (seconds)')
ylabel(ax2,'Vertical Speed (ft/min)')
if showlgd
    legend(ax2,lgdcell,'NumColumns',nt)
end
title(ax2,titstr)

grid(ax3,'on')
xlabel(ax3,'Time (seconds)')
ylabel(ax3,'Elevator Deflection (%)')
if showlgd
    legend(ax3,lgdcell,'NumColumns',nt)
end
title(ax3,titstr)

grid(ax4,'on')
xlabel(ax4,'Takeoff Distance (ft)')
ylabel(ax4,'Takeoff Speed (KIAS)')
if showlgd
    legend(ax4,lgdcell,'NumColumns',nt)
end
title(ax4,titstr)

grid(ax5,'on')
xlabel(ax5,'Time (seconds)')
ylabel(ax5,'Vertical Speed (ft/min)')
if showlgd
    legend(ax5,lgdcell,'NumColumns',nt)
end
title(ax5,titstr)

grid(ax6,'on')
xlabel(ax6,'Time (seconds)')
ylabel(ax6,'Altitude (ft)')
if showlgd
    legend(ax6,lgdcell,'NumColumns',nt)
end
title(ax6,titstr)

grid(ax7,'on')
xlabel(ax7,'Time (seconds)')
ylabel(ax7,'Z Body Force (N)')
if showlgd
    legend(ax7,lgdcell,'NumColumns',nt)
end
title(ax7,titstr)

grid(ax8,'on')
xlabel(ax8,'Indicated Airspeed (kn)')
ylabel(ax8,'Z Body Force (N)')
if showlgd
    legend(ax8,lgdcell,'NumColumns',nt)
end
title(ax8,titstr)

%% Condition #2: Descent & Landing
close all
cond = 2;   %condition number

figure
hold on
ax1 = gca;  %get axes for altitude plot

figure
hold on
ax2 = gca;  %get axes for vertical speed plot

figure
hold on
ax3 = gca;  %get axes for elevator plot

figure
hold on
ax4 = gca;  %get axes for airspeed plot

for k = 1:nt
    for i = pilots
        if isempty(s(i,cond,k).time)    %check for missing data
            %Plot dummy point
            plot(0,0,'k--')
        else
            alt = s(i,cond,k).alt/.3048;    %altitude [ft]
            t = s(i,cond,k).time;   %time [seconds]
            hdot = s(i,cond,k).hdot/.3048*60;   %sink rate [ft/min]
            de = s(i,cond,k).eta;   %elevator deflection [-1 to 1]
            V = s(i,cond,k).VIAS;   %indicated airspeed [kn]
            
            plot(ax4,alt,V,'Color',cp(i-pilots(1)+1,:)*ct(k))
            
            inds = and(alt<900,alt>100);   %indices of interest to plot
            alt = alt(inds);
            t = t(inds);
            hdot = hdot(inds);
            de = de(inds);
            V = V(inds);
            
            %Align data with origin of plot
            t = t-t(1);
            
            plot(ax1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax2,t,hdot,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax3,t,de,'Color',cp(i-pilots(1)+1,:)*ct(k))
        end
    end
end
t1 = linspace(0,800/sr(2)*60,100);  %time vector [seconds]
t2 = linspace(0,800/sr(1)*60,100);  %time vector [seconds]
altmax = 900-t1/60*sr(2);           %altitude at maximum sink rate [ft]
altmin = 900-t2/60*sr(1);           %altitude at minimum sink rate [ft]

plot(ax1,t1,altmax,'k--','LineWidth',2)
plot(ax1,t2,altmin,'k--','LineWidth',2)
plot(ax2,t1,-ones(size(t1))*sr(2),'k--','LineWidth',2)
plot(ax2,t2,-ones(size(t2))*sr(1),'k--','LineWidth',2)
plot(ax4,[0 1000],ones(1,2)*Vapp(1),'k--','LineWidth',2)
plot(ax4,[0 1000],ones(1,2)*Vapp(2),'k--','LineWidth',2)

%Plot Formatting
titstr = sprintf('Condition #%i: %s',cond,condstr{cond});

grid(ax1,'on')
xlabel(ax1,'Time (seconds)')
ylabel(ax1,'Altitude (ft)')
if showlgd
    legend(ax1,lgdcell,'NumColumns',nt,'Location','southeast')
end
title(ax1,titstr)

grid(ax2,'on')
xlabel(ax2,'Time (seconds)')
ylabel(ax2,'Vertical Speed (ft/min)')
if showlgd
    legend(ax2,lgdcell,'NumColumns',nt)
end
title(ax2,titstr)

grid(ax3,'on')
xlabel(ax3,'Time (seconds)')
ylabel(ax3,'Elevator Deflection (%)')
if showlgd
    legend(ax3,lgdcell,'NumColumns',nt)
end
title(ax3,titstr)

ax4.XDir = 'reverse';
ax4.XLim = [0 1000];
grid(ax4,'on')
xlabel(ax4,'Altitude (ft)')
ylabel(ax4,'IAS (kn)')
if showlgd
    legend(ax4,lgdcell,'NumColumns',nt)
end
title(ax4,titstr)

%% Condition #3: S-Turn
close all
cond = 3;   %condition number

figure
hold on
ax1 = gca;  %get axes for altitude plot

figure
hold on
ax2 = gca;  %get axes for vertical speed plot

figure
hold on
ax3 = gca;  %get axes for bank angle plot

figure
hold on
ax4 = gca;  %get axes for sideslip plot

figure
hold on
ax5 = gca;  %get axes for roll rate plot

figure
hold on
ax6 = gca;  %get axes for velocity plot

f = figure; %figure handle for controls plots
subplot(2,1,1)
hold on
ax7 = gca;  %get axes for aileron plot
subplot(2,1,2)
hold on
ax10 = gca; %get axes for rudder plot

figure
hold on
ax8 = gca;  %get axes for heading plot

figure
hold on
ax9 = gca;  %get axes for flight path plot

fc = figure;    %figure handle for input plots
subplot(4,1,1)
hold on
axc1 = gca; %get axes for stick-x
subplot(4,1,2)
hold on
axc2 = gca; %get axes for stick-y
subplot(4,1,3)
hold on
axc3 = gca; %get axes for pedals
subplot(4,1,4)
hold on
axc4 = gca; %get axes for throttle lever

fi = figure;    %figure handle for instruction plots
subplot(2,2,1)
hold on
axi1 = gca; %get axes for altitude
subplot(2,2,2)
hold on
axi2 = gca; %get axes for airspeed
subplot(2,2,3)
hold on
axi3 = gca; %get axes for sideslip
subplot(2,2,4)
hold on
axi4 = gca; %get axes for bank angle

for k = 1:nt
    for i = pilots
        if isempty(s(i,cond,k).time)    %check for missing data
            %Plot dummy point
            plot(0,0,'k--')
        else
            alt = s(i,cond,k).alt/.3048;    %altitude [ft]
            t = s(i,cond,k).time;   %time [seconds]
            hdot = s(i,cond,k).hdot/.3048*60;   %sink rate [ft/min]
            phi = s(i,cond,k).phi;      %bank angle [degrees]
            beta = s(i,cond,k).beta;    %sideslip angle [degrees]
            p = s(i,cond,k).p;          %roll rate [deg/s]
            V = s(i,cond,k).VIAS;       %indicated airspeed [kn]
            da = s(i,cond,k).xsi;       %aileron deflection [%]
            dr = s(i,cond,k).zeta;      %rudder deflection [%]
            psi = s(i,cond,k).psi;      %TRUE heading angle [degrees]
            E = s(i,cond,k).east/.3048; %east position [ft]
            N = s(i,cond,k).north/.3048;%north position [ft]
            stickx = s(i,cond,k).x_stick;
            sticky = s(i,cond,k).y_stick;
            pedal = s(i,cond,k).pedal;
            throt = s(i,cond,k).l_throt;
            
            ind = find(phi>=10,1);      %finding start of S-turn
            ind = ind-25;               %assume pilots bank 10deg/s
            if ind < 1
                ind = 1;
            end
            
            alt = alt(ind:end);
            t = t(ind:end);
            hdot = hdot(ind:end);
            phi = phi(ind:end);
            beta = beta(ind:end);
            p = p(ind:end);
            V = V(ind:end);
            da = da(ind:end);
            dr = dr(ind:end);
            psi = psi(ind:end);
            E = E(ind:end);
            N = N(ind:end);
            stickx = stickx(ind:end);
            sticky = sticky(ind:end);
            pedal = pedal(ind:end);
            throt = throt(ind:end);
            
            %Align data with origin of plot
            t = t-t(1);
            E = E-E(1);
            N = N-N(1);
            
            psi = mod(psi-90,360)+90; %shift to be between 0 and 360 degrees
            
            plot(ax1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax2,t,hdot,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax3,t,phi,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax4,t,beta,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax5,t,p,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax6,t,V,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax7,t,da,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax8,t,psi,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax9,E,N,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax10,t,dr,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axc1,t,stickx,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axc2,t,sticky,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axc3,t,pedal,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axc4,t,throt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi2,t,V,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi3,t,beta,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi4,t,phi,'Color',cp(i-pilots(1)+1,:)*ct(k))
        end
    end
end

%Plotting +/- 100 ft altitude limits (from instructions)
t1 = ax1.XLim;
altmax = 4100*ones(1,2);
altmin = 3900*ones(1,2);
plot(ax1,t1,altmax,'k--','LineWidth',2)
plot(ax1,t1,altmin,'k--','LineWidth',2)
plot(axi1,t1,altmax,'k--','LineWidth',2)
plot(axi1,t1,altmin,'k--','LineWidth',2)

%Plotting vertical speed limits
t2 = ax2.XLim;
hdotmax = 100*ones(1,2);    %100 ft/min vertical speed (arbitrary)
hdotmin = -100*ones(1,2);
plot(ax2,t2,hdotmax,'k--','LineWidth',2)
plot(ax2,t2,hdotmin,'k--','LineWidth',2)

%Plotting Ideal Banking
g = 32.2;   %gravitational acceleration [m/s^2]
pid = 10;   %roll rate for ideal flight [deg/s]
phiid = 20;     %ideal bank angle for S-turn [deg]
phimin = 15;    %minimum bank angle allowed for S-turn [deg]
phimax = 25;    %maximum bank angle allowed for S-turn [deg]
psi0 = 180;     %initial heading [deg]
rid = V3^2/g/tand(phiid); %ideal S-turn radius [ft]

sid = idealSTurn(psi0,V3,phiid,pid,4000);
s1 = idealSTurn(psi0,V3,[phimin,-phimax],pid,4000);
s2 = idealSTurn(psi0,V3,[phimax,-phimin],pid,4000);

plot(ax3,sid.time,sid.phi,'k--','LineWidth',2)
% plot(ax3,s1.time,s1.phi,'--','Color',[0.3,0.3,0.3],'LineWidth',2)
% plot(ax3,s2.time,s2.phi,'--','Color',[0.3,0.3,0.3],'LineWidth',2)
plot(ax8,sid.time,mod(sid.psi-90,360)+90,'k--','LineWidth',2)
plot(ax9,sid.east,sid.north,'k--','LineWidth',2)
% plot(ax9,s1.east,s1.north,'--','Color',[0.3,0.3,0.3],'LineWidth',2)
% plot(ax9,s2.east,s2.north,'--','Color',[0.3,0.3,0.3],'LineWidth',2)

plot(axi4,sid.time,sid.phi,'k','LineWidth',2)

%Plotting +/- Sideslip
t4 = ax4.XLim;
betamax = 1*ones(1,2);  %1 degree sideslip (arbitrary)
betamin = -1*ones(1,2);
plot(ax4,t4,betamax,'k--','LineWidth',2)
plot(ax4,t4,betamin,'k--','LineWidth',2)
plot(axi3,t4,betamax,'k--','LineWidth',2)
plot(axi3,t4,betamin,'k--','LineWidth',2)

%Plotting Airspeed Limits
t6 = ax6.XLim;
Vmax = Va*ones(1,2);        %maximum maneuvering speed [ft/s]
plot(ax6,t6,Vmax,'r--','LineWidth',2)
plot(ax6,t6,[V3,V3],'k--','LineWidth',2)
plot(axi2,t6,Vmax,'r--','LineWidth',2)
plot(axi2,t6,[V3,V3],'k--','LineWidth',2)

%Plot formatting
titstr = sprintf('Condition #%i: %s',cond,condstr{cond});

grid(ax1,'on')
xlabel(ax1,'Time (seconds)')
ylabel(ax1,'Altitude (ft)')
if showlgd
    legend(ax1,lgdcell,'NumColumns',nt,'Location','southeast')
end
title(ax1,titstr)

grid(ax2,'on')
xlabel(ax2,'Time (seconds)')
ylabel(ax2,'Vertical Speed (ft/min)')
if showlgd
    legend(ax2,lgdcell,'NumColumns',nt)
end
title(ax2,titstr)

grid(ax3,'on')
xlabel(ax3,'Time (seconds)')
ylabel(ax3,'Bank Angle (\circ)')
if showlgd
    legend(ax3,lgdcell,'NumColumns',nt)
end
title(ax3,titstr)

grid(ax4,'on')
xlabel(ax4,'Time (seconds)')
ylabel(ax4,'Sideslip Angle (\circ)')
if showlgd
    legend(ax4,lgdcell,'NumColumns',nt)
end
title(ax4,titstr)

grid(ax5,'on')
xlabel(ax5,'Time (seconds)')
ylabel(ax5,'Body-Axis Roll Rate (\circ/s)')
if showlgd
    legend(ax5,lgdcell,'NumColumns',nt)
end
title(ax5,titstr)

grid(ax6,'on')
xlabel(ax6,'Time (seconds)')
ylabel(ax6,'Indicated Airspeed (knots)')
if showlgd
    legend(ax6,lgdcell,'NumColumns',nt)
end
title(ax6,titstr)

grid(ax7,'on')
xlabel(ax7,'Time (seconds)')
ylabel(ax7,'Aileron Deflection (%)')
if showlgd
    legend(ax7,lgdcell,'NumColumns',nt)
end
sgtitle(f,titstr)

grid(ax8,'on')
xlabel(ax8,'Time (seconds)')
ylabel(ax8,'Aircraft Heading (\circ)')
if showlgd
    legend(ax8,lgdcell,'NumColumns',nt)
end
title(ax8,titstr)

grid(ax9,'on')
xlabel(ax9,'East Position (ft)')
ylabel(ax9,'North Position (ft)')
if showlgd
    legend(ax9,lgdcell,'NumColumns',nt)
end
title(ax9,titstr)
axis(ax9,'equal')

grid(ax10,'on')
xlabel(ax10,'Time (seconds)')
ylabel(ax10,'Rudder Deflection (%)')
if showlgd
    legend(ax10,lgdcell,'NumColumns',nt)
end

sgtitle(fc,titstr)
grid(axc1,'on')
grid(axc2,'on')
grid(axc3,'on')
grid(axc4,'on')
xlabel(axc1,'Time (seconds)')
xlabel(axc2,'Time (seconds)')
xlabel(axc3,'Time (seconds)')
xlabel(axc4,'Time (seconds)')
ylabel(axc1,'Stick-X')
ylabel(axc2,'Stick-Y')
ylabel(axc3,'Pedals')
ylabel(axc4,'Throttle')
if showlgd
    legend(axc4,lgdcell,'NumColumns',nt)
end

sgtitle(fi,titstr)
grid(axi1,'on')
grid(axi2,'on')
grid(axi3,'on')
grid(axi4,'on')
xlabel(axi1,'Time (seconds)')
xlabel(axi2,'Time (seconds)')
xlabel(axi3,'Time (seconds)')
xlabel(axi4,'Time (seconds)')
ylabel(axi1,'Altitude (ft)')
ylabel(axi2,'IAS (kn)')
ylabel(axi3,'Sideslip (\circ)')
ylabel(axi4,'Bank Angle (\circ)')

%% Condition #4: 360 Turn
close all
cond = 4;   %condition number

figure
hold on
ax1 = gca;  %get axes for altitude plot

figure
hold on
ax2 = gca;  %get axes for vertical speed plot

figure
hold on
ax3 = gca;  %get axes for bank angle plot

figure
hold on
ax4 = gca;  %get axes for sideslip plot

figure
hold on
ax5 = gca;  %get axes for roll rate plot

figure
hold on
ax6 = gca;  %get axes for velocity plot

f = figure; %figure handle for controls plots
subplot(2,1,1)
hold on
ax7 = gca;  %get axes for aileron plot
subplot(2,1,2)
hold on
ax10 = gca; %get axes for rudder plot

figure
hold on
ax8 = gca;  %get axes for heading plot

figure
hold on
ax9 = gca;  %get axes for flight path plot

fi = figure;    %figure handle for instruction plots
subplot(2,2,1)
hold on
axi1 = gca; %get axes for altitude
subplot(2,2,2)
hold on
axi2 = gca; %get axes for airspeed
subplot(2,2,3)
hold on
axi3 = gca; %get axes for sideslip
subplot(2,2,4)
hold on
axi4 = gca; %get axes for bank angle

for k = 1:nt
    for i = pilots
        if isempty(s(i,cond,k).time)    %check for missing data
            %Plot dummy point
            plot(0,0,'k--')
        else
            alt = s(i,cond,k).alt/.3048;    %altitude [ft]
            t = s(i,cond,k).time;   %time [seconds]
            hdot = s(i,cond,k).hdot/.3048*60;   %sink rate [ft/min]
            phi = s(i,cond,k).phi;      %bank angle [degrees]
            beta = s(i,cond,k).beta;    %sideslip angle [degrees]
            p = s(i,cond,k).p;          %roll rate [deg/s]
            V = s(i,cond,k).VIAS;       %indicated airspeed [kn]
            da = s(i,cond,k).xsi;       %aileron deflection [%]
            dr = s(i,cond,k).zeta;      %rudder deflection [%]
            psi = s(i,cond,k).psi;      %TRUE heading angle [degrees]
            E = s(i,cond,k).east/.3048; %east position [ft]
            N = s(i,cond,k).north/.3048;%north position [ft]
            
%             ind = find(phi>=10,1);      %finding start of S-turn
%             ind = ind-25;               %assume pilots bank 10deg/s
%             if ind < 1
%                 ind = 1;
%             end
            ind = 1;

            alt = alt(ind:end);
            t = t(ind:end);
            hdot = hdot(ind:end);
            phi = phi(ind:end);
            beta = beta(ind:end);
            p = p(ind:end);
            V = V(ind:end);
            da = da(ind:end);
            dr = dr(ind:end);
            psi = psi(ind:end);
            E = E(ind:end);
            N = N(ind:end);
            
            if t(end) > 300 %shorten data if longer than 2.5 minutes
                ind2 = find(t>=150,1);
                
                alt = alt(1:ind2);
                t = t(1:ind2);
                hdot = hdot(1:ind2);
                phi = phi(1:ind2);
                beta = beta(1:ind2);
                p = p(1:ind2);
                V = V(1:ind2);
                da = da(1:ind2);
                dr = dr(1:ind2);
                psi = psi(1:ind2);
                E = E(1:ind2);
                N = N(1:ind2);
            end
            
            %Align data with origin of plot
            t = t-t(1);
            E = E-E(1);
            N = N-N(1);
            
            psi = mod(psi,360); %shift to be between 0 and 360 degrees
            
            plot(ax1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax2,t,hdot,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax3,t,phi,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax4,t,beta,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax5,t,p,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax6,t,V,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax7,t,da,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax8,t,psi,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax9,E,N,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax10,t,dr,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi2,t,V,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi3,t,beta,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi4,t,phi,'Color',cp(i-pilots(1)+1,:)*ct(k))
        end
    end
end

%Plotting +/- 100 ft altitude limits (from instructions)
t1 = ax1.XLim;
altmax = 3100*ones(1,2);
altmin = 2900*ones(1,2);
plot(ax1,t1,altmax,'k--','LineWidth',2)
plot(ax1,t1,altmin,'k--','LineWidth',2)
plot(axi1,t1,altmax,'k--','LineWidth',2)
plot(axi1,t1,altmin,'k--','LineWidth',2)

%Plotting vertical speed limits
t2 = ax2.XLim;
hdotmax = 100*ones(1,2);    %100 ft/min vertical speed (arbitrary)
hdotmin = -100*ones(1,2);
plot(ax2,t2,hdotmax,'k--','LineWidth',2)
plot(ax2,t2,hdotmin,'k--','LineWidth',2)

%Plotting Ideal Banking
g = 32.2;   %gravitational acceleration [m/s^2]
pid = 10;   %roll rate for ideal flight [deg/s]
phiid = 20;     %ideal bank angle for S-turn [deg]
phimin = 15;    %minimum bank angle allowed for S-turn [deg]
phimax = 25;    %maximum bank angle allowed for S-turn [deg]
psi0 = 55;     %initial heading [deg]
rid = V4^2/g/tand(phiid); %ideal S-turn radius [ft]

sid = idealSTurn(psi0,V4,[phiid,phiid],pid,4000);
s1 = idealSTurn(psi0,V4,[phimin,phimin],pid,4000);
s2 = idealSTurn(psi0,V4,[phimax,phimax],pid,4000);

% plot(ax3,sid.time,sid.phi,'k--','LineWidth',2)
% plot(ax3,s1.time,s1.phi,'--','Color',[0.3,0.3,0.3],'LineWidth',2)
% plot(ax3,s2.time,s2.phi,'--','Color',[0.3,0.3,0.3],'LineWidth',2)
% plot(ax8,sid.time,sid.psi,'k--','LineWidth',2)
% plot(ax9,sid.east,sid.north,'k--','LineWidth',2)
% plot(ax9,s1.east,s1.north,'--','Color',[0.3,0.3,0.3],'LineWidth',2)
% plot(ax9,s2.east,s2.north,'--','Color',[0.3,0.3,0.3],'LineWidth',2)

% plot(axi4,sid.time,sid.phi,'k--','LineWidth',2)
plot(axi4,axi4.XLim,ones(1,2)*phiid,'k--','LineWidth',2)

%Plot wind direction + speed
x0 = 0.35;      %arrow start x-coord [0-1]
y0 = 0.65;      %arrow start y-coord [0-1]
r = 0.1;        %arrow length

x = x0 + r*[0,-sind(whead4)];
y = y0 + r*[0,-cosd(whead4)];
str = sprintf('%i kn from %i\\circ',wspeed4,whead4);
annotation(ax9.Parent,'textarrow',x,y,'String',str)

%Plotting +/- Sideslip
t4 = ax4.XLim;
betamax = 1*ones(1,2);  %1 degree sideslip (arbitrary)
betamin = -1*ones(1,2);
plot(ax4,t4,betamax,'k--','LineWidth',2)
plot(ax4,t4,betamin,'k--','LineWidth',2)
plot(axi3,t4,betamax,'k--','LineWidth',2)
plot(axi3,t4,betamin,'k--','LineWidth',2)

%Plotting Airspeed Limits
t6 = ax6.XLim;
Vmax = Va*ones(1,2);    %maximum maneuvering speed [ft/s]
plot(ax6,t6,Vmax,'r--','LineWidth',2)
plot(ax6,t6,[V4,V4],'k--','LineWidth',2)
plot(axi2,t6,Vmax,'r--','LineWidth',2)
plot(axi2,t6,[V4,V4],'k--','LineWidth',2)

%Plot formatting
titstr = sprintf('Condition #%i: %s',cond,condstr{cond});

grid(ax1,'on')
xlabel(ax1,'Time (seconds)')
ylabel(ax1,'Altitude (ft)')
if showlgd
    legend(ax1,lgdcell,'NumColumns',nt,'Location','southeast')
end
title(ax1,titstr)

grid(ax2,'on')
xlabel(ax2,'Time (seconds)')
ylabel(ax2,'Vertical Speed (ft/min)')
if showlgd
    legend(ax2,lgdcell,'NumColumns',nt)
end
title(ax2,titstr)

grid(ax3,'on')
xlabel(ax3,'Time (seconds)')
ylabel(ax3,'Bank Angle (\circ)')
if showlgd
    legend(ax3,lgdcell,'NumColumns',nt)
end
title(ax3,titstr)

grid(ax4,'on')
xlabel(ax4,'Time (seconds)')
ylabel(ax4,'Sideslip Angle (\circ)')
if showlgd
    legend(ax4,lgdcell,'NumColumns',nt)
end
title(ax4,titstr)

grid(ax5,'on')
xlabel(ax5,'Time (seconds)')
ylabel(ax5,'Body-Axis Roll Rate (\circ/s)')
if showlgd
    legend(ax5,lgdcell,'NumColumns',nt)
end
title(ax5,titstr)

grid(ax6,'on')
xlabel(ax6,'Time (seconds)')
ylabel(ax6,'Indicated Airspeed (knots)')
if showlgd
    legend(ax6,lgdcell,'NumColumns',nt)
end
title(ax6,titstr)

grid(ax7,'on')
xlabel(ax7,'Time (seconds)')
ylabel(ax7,'Aileron Deflection (%)')
if showlgd
    legend(ax7,lgdcell,'NumColumns',nt)
end
sgtitle(f,titstr)

grid(ax8,'on')
xlabel(ax8,'Time (seconds)')
ylabel(ax8,'Aircraft Heading (\circ)')
if showlgd
    legend(ax8,lgdcell,'NumColumns',nt)
end
title(ax8,titstr)

grid(ax9,'on')
xlabel(ax9,'East Position (ft)')
ylabel(ax9,'North Position (ft)')
if showlgd
    legend(ax9,lgdcell,'NumColumns',nt)
end
title(ax9,titstr)
axis(ax9,'equal')

grid(ax10,'on')
xlabel(ax10,'Time (seconds)')
ylabel(ax10,'Rudder Deflection (%)')
if showlgd
    legend(ax10,lgdcell,'NumColumns',nt)
end

sgtitle(fi,titstr)
grid(axi1,'on')
grid(axi2,'on')
grid(axi3,'on')
grid(axi4,'on')
xlabel(axi1,'Time (seconds)')
xlabel(axi2,'Time (seconds)')
xlabel(axi3,'Time (seconds)')
xlabel(axi4,'Time (seconds)')
ylabel(axi1,'Altitude (ft)')
ylabel(axi2,'IAS (kn)')
ylabel(axi3,'Sideslip (\circ)')
ylabel(axi4,'Bank Angle (\circ)')

%% Condition #5: Crosswind Landing
close all
cond = 5;   %condition number

figure
hold on
ax1 = gca;  %get axes for altitude plot

figure
hold on
ax2 = gca;  %get axes for vertical speed plot

f = figure; %get handle for controls plots figure
subplot(2,1,1)
hold on
ax3a = gca; %get axes for elevator plot
subplot(2,1,2)
hold on
ax3b = gca; %get axes for rudder plot

figure
hold on
ax4 = gca;  %get axes for airspeed plot

figure
hold on
ax5 = gca;  %get axes for sideslip plot

fi = figure;    %figure handle for instruction plots
subplot(2,2,1)
hold on
axi1 = gca;     %altitude plot
subplot(2,2,2)
hold on
axi2 = gca;     %heading plot
subplot(2,2,3)
hold on
axi3 = gca;     %aileron plot
subplot(2,2,4)
hold on
axi4 = gca;     %rudder plot

for k = 1:nt
    for i = pilots
        if isempty(s(i,cond,k).time)    %check for missing data
            %Plot dummy point
            plot(0,0,'k--')
        else
            alt = s(i,cond,k).alt/.3048;    %altitude [ft]
            t = s(i,cond,k).time;   %time [seconds]
            hdot = s(i,cond,k).hdot/.3048*60;   %sink rate [ft/min]
            de = s(i,cond,k).eta;   %elevator deflection [-1 to 1]
            dr = s(i,cond,k).zeta;  %rudder deflection [-1 to 1]
            da = s(i,cond,k).xsi;   %aileron deflection [-1 to 1]
            psi = s(i,cond,k).psi;  %heading angle [degrees]
            V = s(i,cond,k).VIAS;   %indicated airspeed [kn]
            beta = s(i,cond,k).beta;    %sideslip angle [degrees]
            
            plot(axi1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi2,t,psi,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi3,t,da,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi4,t,dr,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax5,t,beta,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax4,t,V,'Color',cp(i-pilots(1)+1,:)*ct(k))

            inds = and(alt<900,alt>100);   %indices of interest to plot
            alt = alt(inds);
            t = t(inds);
            hdot = hdot(inds);
            de = de(inds);
            dr = dr(inds);
            V = V(inds);
            beta = beta(inds);
            
            %Align data with origin of plot
            t = t-t(1);
            
            plot(ax1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax2,t,hdot,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax3a,t,de,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax3b,t,dr,'Color',cp(i-pilots(1)+1,:)*ct(k))
        end
    end
end
t1 = linspace(0,800/sr(2)*60,100);  %time vector [seconds]
t2 = linspace(0,800/sr(1)*60,100);  %time vector [seconds]
altmax = 900-t1/60*sr(2);           %altitude at maximum sink rate [ft]
altmin = 900-t2/60*sr(1);           %altitude at minimum sink rate [ft]

plot(ax1,t1,altmax,'k--','LineWidth',2)
plot(ax1,t2,altmin,'k--','LineWidth',2)
plot(ax2,t1,-ones(size(t1))*sr(2),'k--','LineWidth',2)
plot(ax2,t2,-ones(size(t2))*sr(1),'k--','LineWidth',2)

plot(ax4,ax4.XLim,ones(1,2)*Vapp(1),'k--','LineWidth',2)
plot(ax4,ax4.XLim,ones(1,2)*Vapp(2),'k--','LineWidth',2)

%Plot Formatting
titstr = sprintf('Condition #%i: %s',cond,condstr{cond});

grid(ax1,'on')
xlabel(ax1,'Time (seconds)')
ylabel(ax1,'Altitude (ft)')
if showlgd
    legend(ax1,lgdcell,'NumColumns',nt,'Location','southeast')
end
title(ax1,titstr)

grid(ax2,'on')
xlabel(ax2,'Time (seconds)')
ylabel(ax2,'Vertical Speed (ft/min)')
if showlgd
    legend(ax2,lgdcell,'NumColumns',nt)
end
title(ax2,titstr)

grid(ax3a,'on')
xlabel(ax3a,'Time (seconds)')
ylabel(ax3a,'Elevator Deflection (%)')
if showlgd
    legend(ax3a,lgdcell,'NumColumns',nt)
end
sgtitle(f,titstr)

grid(ax3b,'on')
xlabel(ax3b,'Time (seconds)')
ylabel(ax3b,'Rudder Deflection (%)')
if showlgd
    legend(ax3b,lgdcell,'NumColumns',nt)
end

grid(ax4,'on')
xlabel(ax4,'Time (seconds)')
ylabel(ax4,'IAS (kn)')
if showlgd
    legend(ax4,lgdcell,'NumColumns',nt)
end
title(ax4,titstr)

grid(ax5,'on')
xlabel(ax5,'Time (seconds)')
ylabel(ax5,'Sideslip (\circ)')
if showlgd
    legend(ax5,lgdcell,'NumColumns',nt)
end
title(ax5,titstr)

sgtitle(fi,titstr)
grid(axi1,'on')
grid(axi2,'on')
grid(axi3,'on')
grid(axi4,'on')
xlabel(axi1,'Time (seconds)')
xlabel(axi2,'Time (seconds)')
xlabel(axi3,'Time (seconds)')
xlabel(axi4,'Time (seconds)')
ylabel(axi1,'Altitude (ft)')
ylabel(axi2,'Heading (\circ)')
ylabel(axi3,'Aileron (%)')
ylabel(axi4,'Rudder (%)')

%% Condition #6: Short Field Takeoff
% fuelpc = 0.5;                   %gas tank fullness when initialized
% Wto = 7429 - (1-fuelpc)*655;    %weight when initialized

close all
cond = 6;   %condition number

figure
hold on
ax1 = gca;  %get axes for altitude plot

figure
hold on
ax2 = gca;  %get axes for vertical speed plot

figure
hold on
ax3 = gca;  %get axes for elevator plot

figure
hold on
ax4 = gca;  %get axes for climb rate plot

figure
hold on
ax5 = gca;  %get axes for takeoff performance plot

figure
hold on
ax6 = gca;  %get axes for takeoff performance plot

crf0 = zeros(1,nt*length(pilots));  %preallocation
crf1 = zeros(1,nt*length(pilots));
for k = 1:nt
    for i = pilots
        if isempty(s(i,cond,k).time)    %check for missing data
            %Plot dummy point
            plot(0,0,'k--')
        else
            alt = s(i,cond,k).alt/.3048;    %altitude [ft]
            t = s(i,cond,k).time;   %time [seconds]
            hdot = s(i,cond,k).hdot/.3048*60;   %climb rate [ft/min]
            de = s(i,cond,k).eta;   %elevator deflection [-1 to 1]
            df = s(i,cond,k).d_flap;    %flap stage #
            N = s(i,cond,k).north/.3048;    %north position [ft]
            E = s(i,cond,k).east/.3048;     %east position [ft]
            V = s(i,cond,k).VIAS;       %indicated airspeed [kn]
            nz = s(i,cond,k).nz;        %load factor (g)
            L = -s(i,cond,k).Z;         %lift force [N]
            
            %Trim data to just the climb
            ind1 = find(alt>150,1);     %first altitude > 150
            ind2 = find(alt>1350,1)-1;  %last altitude < 1350 (climbing)
            
            %Estimate average climb rate with and without flaps
            indto1 = find(alt>10,1);
            indto = find(hdot>30,1);        %find takeoff from runway
            indf = find(df>0,1,'last');    %find when flaps are retracted
            
            j = i-3 + (k-1)*length(pilots); %index used for ROC bar plot
            
            %Rate-of-climb for flaps 1 and flaps 0 (retracted)
            crf1(j) = (alt(indf)-alt(indto))/(t(indf)-t(indto))*60;
            crf0(j) = (alt(ind2)-alt(indf))/(t(ind2)-t(indf))*60;
            cr6(j) = (alt(ind2)-alt(ind1))/(t(ind2)-t(ind1))*60;
            altf(j) = alt(indf);    %altitude when flaps retracted [ft]
            
            tTO = t;    %save unaltered time vector for takeoff plot
            hdotTO = hdot;
%             fprintf('%.2f s\t%.2f ft\t%.2f KIAS\t%.0f ft\n',t(indto1)-t(indto),...
%                 alt(indto1)-alt(indto),V(indto1)-V(indto),...
%                 norm([N(indto1),E(indto1)])-norm([N(indto),E(indto)]))
            fprintf('%.2f s\t%.2f ft\t%.2f KIAS\t%.0f ft\n',t(indto),...
                alt(indto),V(indto),...
                norm([N(indto),E(indto)]))
            
            alt = alt(ind1:ind2);
            t = t(ind1:ind2);
            hdot = hdot(ind1:ind2);
            de = de(ind1:ind2);
            df = df(ind1:ind2);
            
            %Takeoff performance
            VTO = V(indto);     %takeoff speed [KIAS]
            dTO = norm([N(indto),E(indto)]);    %takeoff distance [ft]
            
            %Align data with origin of plot
            t = t-t(1);
            
            plot(ax1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax2,t,hdot,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax3,t,de,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax5,dTO,VTO,'*','Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax6,tTO,hdotTO,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax6,tTO(indto),hdotTO(indto),'*','Color',cp(i-pilots(1)+1,:)*ct(k))
        end
    end
end

crf1avg = mean(crf1);   %average flaps 1 climb rate [ft/min]
crf0avg = mean(crf0);   %average no flaps climb rate [ft/min]
fprintf('Average climb rate, clean:  \t%.1f fpm\n',crf0avg)
fprintf('Average climb rate, flaps 1:\t%.1f fpm\n',crf1avg)
altfavg = mean(altf);   %average altitude
fprintf('Average altitude flaps retracted:\t%.0f ft\n',altfavg)

%Create bar plot
X = categorical({'Flaps Retracted','Flaps 1'});
X = reordercats(X,{'Flaps Retracted','Flaps 1'});
b = bar(ax4,X,[crf0;crf1]);
for k = 1:nt    %setting bar colors
    for i = pilots
        j = i-3 + (k-1)*length(pilots);
        b(j).FaceColor = 'flat';
        b(j).CData = cp(i-pilots(1)+1,:)*ct(k);
    end
end
% plot(ax4,ax4.XLim,[cr(1),cr(1)],'k--','LineWidth',2)

%Plotting Expected Maximum/Minimum Climb Rates
t1 = linspace(0,1200/cr(2)*60,100); %time vector [seconds]
t2 = linspace(0,1200/cr(1)*60,100); %time vector [seconds]
altmax = t1/60*cr(2)+150;           %altitude at maximum climb rate [ft]
altmin = t2/60*cr(1)+150;           %altitude at minimum climb rate [ft]

plot(ax1,t1,altmax,'k--','LineWidth',2)
plot(ax1,t2,altmin,'k--','LineWidth',2)
plot(ax2,t1,ones(size(t1))*cr(2),'k--','LineWidth',2)
plot(ax2,t2,ones(size(t2))*cr(1),'k--','LineWidth',2)

plot(ax5,dtosf,Vtosf,'ko','MarkerSize',10,'LineWidth',2)

%Plot Formatting
titstr = sprintf('Condition #%i: %s',cond,condstr{cond});

grid(ax1,'on')
xlabel(ax1,'Time (seconds)')
ylabel(ax1,'Altitude (ft)')
if showlgd
    legend(ax1,lgdcell,'NumColumns',nt,'Location','southeast')
end
title(ax1,titstr)

grid(ax2,'on')
xlabel(ax2,'Time (seconds)')
ylabel(ax2,'Vertical Speed (ft/min)')
if showlgd
    legend(ax2,lgdcell,'NumColumns',nt)
end
title(ax2,titstr)

grid(ax3,'on')
xlabel(ax3,'Time (seconds)')
ylabel(ax3,'Elevator Deflection (%)')
if showlgd
    legend(ax3,lgdcell,'NumColumns',nt)
end
title(ax3,titstr)

grid(ax4,'on')
ylabel(ax4,'Average Climb Rate (ft/min)')
title(ax4,titstr)

grid(ax5,'on')
xlabel(ax5,'Takeoff Distance (ft)')
ylabel(ax5,'Takeoff Speed (KIAS)')
if showlgd
    legend(ax5,lgdcell,'NumColumns',nt)
end
title(ax5,titstr)

%% Condition #7: Landing w/ Flaps
close all
cond = 7;   %condition number

figure
hold on
ax1 = gca;  %get axes for altitude plot

figure
hold on
ax2 = gca;  %get axes for vertical speed plot

f1 = figure;    %get handle for controls plot figure
subplot(2,1,1)
hold on
ax3a = gca;     %get axes for elevator plot
subplot(2,1,2)
hold on
ax3b = gca;     %get axes for flap plot

f2 = figure;    %get handle for altitude per flap stage figure
subplot(3,1,1)
hold on
ax4a = gca;     %get axes for flaps 1
subplot(3,1,2)
hold on
ax4b = gca;     %get axes for flaps 2
subplot(3,1,3)
hold on
ax4c = gca;     %get axes for flaps 3

figure
hold on
ax5 = gca;  %get axes for flap vs altitude plot
figure
hold on
ax6 = gca;  %get axes for velocity vs altitude plot

for k = 1:nt
    for i = pilots
        if isempty(s(i,cond,k).time)    %check for missing data
            %Plot dummy point
            plot(0,0,'k--')
        else
            alt = s(i,cond,k).alt/.3048;    %altitude [ft]
            t = s(i,cond,k).time;   %time [seconds]
            hdot = s(i,cond,k).hdot/.3048*60;   %sink rate [ft/min]
            de = s(i,cond,k).eta;   %elevator deflection [-1 to 1]
            df = s(i,cond,k).d_flap;    %flap deployment [stage #]
            V = s(i,cond,k).VIAS;   %airspeed [kn]
            
            %Separating altitude vs time by flap stage
            ind0 = find(df==0,1,'last');    %find index of last flaps 0
            ind1 = find(df==1,1,'last');    %find index of last flaps 1
            ind2 = find(df==2,1,'last');    %find index of last flaps 2
            ind3 = find(alt>100,1,'last');  %find last time altitude > 100
            
            alt1 = alt(ind0:ind1)-alt(ind0);
            alt2 = alt(ind1:ind2)-alt(ind1);
            alt3 = alt(ind2:ind3)-alt(ind2);
            t1 = t(ind0:ind1)-t(ind0);
            t2 = t(ind1:ind2)-t(ind1);
            t3 = t(ind2:ind3)-t(ind2);
            
            %Finding data between 2 altitudes
            inds = and(alt<900,alt>100);    %indices of interest to plot
            alt = alt(inds);
            t = t(inds);
            hdot = hdot(inds);
            de = de(inds);
            df = df(inds);
            V = V(inds);
            
            %Align data with origin of plot
            t = t-t(1);
            
            plot(ax1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax2,t,hdot,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax3a,t,de,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax3b,t,df,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax4a,t1,alt1,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax4b,t2,alt2,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax4c,t3,alt3,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax5,alt,df,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax6,alt,V,'Color',cp(i-pilots(1)+1,:)*ct(k))
        end
    end
end

%Plotting Expected Minimum/Maximum Sink Rate
t1 = linspace(0,800/sr(2)*60,100);  %time vector [seconds]
t2 = linspace(0,800/sr(1)*60,100);  %time vector [seconds]
altmax = 900-t1/60*sr(2);           %altitude at maximum sink rate [ft]
altmin = 900-t2/60*sr(1);           %altitude at minimum sink rate [ft]

plot(ax1,t1,altmax,'k--','LineWidth',2)
plot(ax1,t2,altmin,'k--','LineWidth',2)
plot(ax2,t1,-ones(size(t1))*sr(2),'k--','LineWidth',2)
plot(ax2,t2,-ones(size(t2))*sr(1),'k--','LineWidth',2)

%Plotting Minimum/Maximum Sink Rates by Flap Stage
t1max = linspace(0,200/sr1(2)*60,100);  %time, flaps 1, max sink [seconds]
t1min = linspace(0,200/sr1(1)*60,100);  %time, flaps 1, min sink [seconds]
t2max = linspace(0,400/sr2(2)*60,100);  %time, flaps 2, max sink [seconds]
t2min = linspace(0,400/sr2(1)*60,100);  %time, flaps 2, min sink [seconds]
t3max = linspace(0,400/sr3(2)*60,100);  %time, flaps 3, max sink [seconds]
t3min = linspace(0,400/sr3(1)*60,100);  %time, flaps 3, min sink [seconds]

alt1max = -t1max/60*sr1(2); %altitude, flaps 1, max sink [ft]
alt1min = -t1min/60*sr1(1); %altitude, flaps 1, min sink [ft]
alt2max = -t2max/60*sr2(2); %altitude, flaps 2, max sink [ft]
alt2min = -t2min/60*sr2(1); %altitude, flaps 2, min sink [ft]
alt3max = -t3max/60*sr3(2); %altitude, flaps 3, max sink [ft]
alt3min = -t3min/60*sr3(1); %altitude, flaps 3, min sink [ft]

plot(ax4a,t1max,alt1max,'k--','LineWidth',2)
plot(ax4a,t1min,alt1min,'k--','LineWidth',2)
plot(ax4b,t2max,alt2max,'k--','LineWidth',2)
plot(ax4b,t2min,alt2min,'k--','LineWidth',2)
plot(ax4c,t3max,alt3max,'k--','LineWidth',2)
plot(ax4c,t3min,alt3min,'k--','LineWidth',2)

%Plotting Instruction Performance Metrics
altid = [900,800,800,400,400,100];
dfid = [1,1,2,2,3,3];
plot(ax5,altid,dfid,'k--','LineWidth',2)

%Plot Formatting
titstr = sprintf('Condition #%i: %s',cond,condstr{cond});

grid(ax1,'on')
xlabel(ax1,'Time (seconds)')
ylabel(ax1,'Altitude (ft)')
if showlgd
    legend(ax1,lgdcell,'NumColumns',nt,'Location','southeast')
end
title(ax1,titstr)

grid(ax2,'on')
xlabel(ax2,'Time (seconds)')
ylabel(ax2,'Vertical Speed (ft/min)')
if showlgd
    legend(ax2,lgdcell,'NumColumns',nt)
end
title(ax2,titstr)

grid(ax3a,'on')
xlabel(ax3a,'Time (seconds)')
ylabel(ax3a,'Elevator Deflection (%)')
if showlgd
    legend(ax3a,lgdcell,'NumColumns',nt)
end
sgtitle(f1,titstr)

grid(ax3b,'on')
xlabel(ax3b,'Time (seconds)')
ylabel(ax3b,'Flap Stage #')
if showlgd
    legend(ax3b,lgdcell,'NumColumns',nt)
end

grid(ax4a,'on')
xlabel(ax4a,'Time (seconds)')
ylabel(ax4a,'Altitude (ft)')
if showlgd
    legend(ax4a,lgdcell,'NumColumns',nt)
end
sgtitle(f2,titstr)
title(ax4a,'Flaps 1')

grid(ax4b,'on')
xlabel(ax4b,'Time (seconds)')
ylabel(ax4b,'Altitude (ft)')
if showlgd
    legend(ax4b,lgdcell,'NumColumns',nt)
end
title(ax4b,'Flaps 2')

grid(ax4c,'on')
xlabel(ax4c,'Time (seconds)')
ylabel(ax4c,'Altitude (ft)')
if showlgd
    legend(ax4c,lgdcell,'NumColumns',nt)
end
title(ax4c,'Flaps 3')

grid(ax5,'on')
xlabel(ax5,'Altitude (ft)')
ylabel(ax5,'Flap Stage #')
ax5.XDir = 'reverse';
if showlgd
    legend(ax5,lgdcell,'NumColumns',nt)
end
title(ax5,titstr)

grid(ax6,'on')
xlabel(ax6,'Altitude (ft)')
ylabel(ax6,'Indicated Airspeed (knots)')
ax6.XDir = 'reverse';
if showlgd
    legend(ax6,lgdcell,'NumColumns',nt)
end
title(ax6,titstr)


%% Condition #8: S-Turn w/ Flaps
close all
cond = 8;   %condition number

figure
hold on
ax1 = gca;  %get axes for altitude plot
figure
hold on
ax2 = gca;  %get axes for vertical speed plot
figure
hold on
ax3 = gca;  %get axes for bank angle plot
figure
hold on
ax4 = gca;  %get axes for sideslip plot
figure
hold on
ax5 = gca;  %get axes for roll rate plot
figure
hold on
ax6 = gca;  %get axes for velocity plot
f = figure; %figure handle for controls plots
subplot(3,1,1)
hold on
ax7 = gca;  %get axes for aileron plot
subplot(3,1,2)
hold on
ax10 = gca; %get axes for rudder plot
subplot(3,1,3)
hold on
ax11 = gca; %get axes for flap plot
figure
hold on
ax8 = gca;  %get axes for heading plot
figure
hold on
ax9 = gca;  %get axes for flight path plot

fi = figure;    %figure handle for instruction plots
subplot(2,2,1)
hold on
axi1 = gca; %get axes for altitude
subplot(2,2,2)
hold on
axi2 = gca; %get axes for airspeed
subplot(2,2,3)
hold on
axi3 = gca; %get axes for sideslip
subplot(2,2,4)
hold on
axi4 = gca; %get axes for bank angle

for k = 1:nt
    for i = pilots
        if isempty(s(i,cond,k).time)    %check for missing data
            %Plot dummy point
            plot(0,0,'k--')
        else
            alt = s(i,cond,k).alt/.3048;    %altitude [ft]
            t = s(i,cond,k).time;   %time [seconds]
            hdot = s(i,cond,k).hdot/.3048*60;   %sink rate [ft/min]
            phi = s(i,cond,k).phi;      %bank angle [degrees]
            beta = s(i,cond,k).beta;    %sideslip angle [degrees]
            p = s(i,cond,k).p;          %roll rate [deg/s]
            V = s(i,cond,k).VIAS;       %indicated airspeed [kn]
            da = s(i,cond,k).xsi;       %aileron deflection [%]
            dr = s(i,cond,k).zeta;      %rudder deflection [%]
            df = s(i,cond,k).d_flap;    %flap setting [stage #]
            psi = s(i,cond,k).psi;      %TRUE heading angle [degrees]
            E = s(i,cond,k).east/.3048; %east position [ft]
            N = s(i,cond,k).north/.3048;%north position [ft]
            
            ind = find(phi>=10,1);      %finding start of S-turn
            ind = ind-25;               %assume pilots bank 10deg/s
            if ind < 1
                ind = 1;
            end
            
            ind2 = find(df>=3,1,'last');    %find end of S-turn
            
            alt = alt(ind:ind2);
            t = t(ind:ind2);
            hdot = hdot(ind:ind2);
            phi = phi(ind:ind2);
            beta = beta(ind:ind2);
            p = p(ind:ind2);
            V = V(ind:ind2);
            da = da(ind:ind2);
            dr = dr(ind:ind2);
            df = df(ind:ind2);
            psi = psi(ind:ind2);
            E = E(ind:ind2);
            N = N(ind:ind2);
            
            %Align data with origin of plot
            t = t-t(1);
            E = E-E(1);
            N = N-N(1);
            
            psi = mod(psi,360); %shift to be between 0 and 360 degrees
            
            plot(ax1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax2,t,hdot,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax3,t,phi,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax4,t,beta,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax5,t,p,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax6,t,V,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax7,t,da,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax8,t,psi,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax9,E,N,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax10,t,dr,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(ax11,t,df,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi1,t,alt,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi2,t,V,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi3,t,beta,'Color',cp(i-pilots(1)+1,:)*ct(k))
            plot(axi4,t,phi,'Color',cp(i-pilots(1)+1,:)*ct(k))
        end
    end
end

%Plotting +/- 100 ft altitude limits (from instructions)
t1 = ax1.XLim;
altmax = 4100*ones(1,2);
altmin = 3900*ones(1,2);
plot(ax1,t1,altmax,'k--','LineWidth',2)
plot(ax1,t1,altmin,'k--','LineWidth',2)
plot(axi1,t1,altmax,'k--','LineWidth',2)
plot(axi1,t1,altmin,'k--','LineWidth',2)

%Plotting vertical speed limits
t2 = ax2.XLim;
hdotmax = 100*ones(1,2);    %100 ft/min vertical speed (arbitrary)
hdotmin = -100*ones(1,2);
plot(ax2,t2,hdotmax,'k--','LineWidth',2)
plot(ax2,t2,hdotmin,'k--','LineWidth',2)

%Plotting Ideal Banking
g = 32.2;       %gravitational acceleration [m/s^2]
pid = 10;       %roll rate for ideal flight [deg/s]
phiid = 20;     %ideal bank angle for S-turn [deg]
phimin = 15;    %minimum bank angle allowed for S-turn [deg]
phimax = 25;    %maximum bank angle allowed for S-turn [deg]
psi0 = 180;     %initial heading [deg]
rid = V8^2/g/tand(phiid); %ideal S-turn radius [ft]

sid = idealSTurn(psi0,V8,phiid,pid,4000);
s1 = idealSTurn(psi0,V8,[phimin,-phimax],pid,4000);
s2 = idealSTurn(psi0,V8,[phimax,-phimin],pid,4000);

plot(ax3,sid.time,sid.phi,'k--','LineWidth',2)
% plot(ax3,s1.time,s1.phi,'--','Color',[0.3,0.3,0.3],'LineWidth',2)
% plot(ax3,s2.time,s2.phi,'--','Color',[0.3,0.3,0.3],'LineWidth',2)
plot(ax8,sid.time,sid.psi,'k--','LineWidth',2)
plot(ax9,sid.east,sid.north,'k--','LineWidth',2)
% plot(ax9,s1.east,s1.north,'--','Color',[0.3,0.3,0.3],'LineWidth',2)
% plot(ax9,s2.east,s2.north,'--','Color',[0.3,0.3,0.3],'LineWidth',2)

plot(axi4,sid.time,sid.phi,'k--','LineWidth',2)

%Plotting +/- Sideslip
t4 = ax4.XLim;
betamax = 1*ones(1,2);  %1 degree sideslip (arbitrary)
betamin = -1*ones(1,2);
plot(ax4,t4,betamax,'k--','LineWidth',2)
plot(ax4,t4,betamin,'k--','LineWidth',2)
plot(axi3,t4,betamax,'k--','LineWidth',2)
plot(axi3,t4,betamin,'k--','LineWidth',2)

%Plotting Airspeed Limits
t6 = ax6.XLim;
Vmax = Vfe*ones(1,2);       %maximum maneuvering speed [ft/s]
plot(ax6,t6,Vmax,'r--','LineWidth',2)
plot(ax6,t6,[V8,V8],'k--','LineWidth',2)
plot(axi2,t6,Vmax,'r--','LineWidth',2)
plot(axi2,t6,[V8,V8],'k--','LineWidth',2)

%Plot formatting
titstr = sprintf('Condition #%i: %s',cond,condstr{cond});

grid(ax1,'on')
xlabel(ax1,'Time (seconds)')
ylabel(ax1,'Altitude (ft)')
if showlgd
    legend(ax1,lgdcell,'NumColumns',nt,'Location','southeast')
end
title(ax1,titstr)

grid(ax2,'on')
xlabel(ax2,'Time (seconds)')
ylabel(ax2,'Vertical Speed (ft/min)')
if showlgd
    legend(ax2,lgdcell,'NumColumns',nt)
end
title(ax2,titstr)

grid(ax3,'on')
xlabel(ax3,'Time (seconds)')
ylabel(ax3,'Bank Angle (\circ)')
if showlgd
    legend(ax3,lgdcell,'NumColumns',nt)
end
title(ax3,titstr)

grid(ax4,'on')
xlabel(ax4,'Time (seconds)')
ylabel(ax4,'Sideslip Angle (\circ)')
if showlgd
    legend(ax4,lgdcell,'NumColumns',nt)
end
title(ax4,titstr)

grid(ax5,'on')
xlabel(ax5,'Time (seconds)')
ylabel(ax5,'Body-Axis Roll Rate (\circ/s)')
if showlgd
    legend(ax5,lgdcell,'NumColumns',nt)
end
title(ax5,titstr)

grid(ax6,'on')
xlabel(ax6,'Time (seconds)')
ylabel(ax6,'Indicated Airspeed (knots)')
if showlgd
    legend(ax6,lgdcell,'NumColumns',nt)
end
title(ax6,titstr)

grid(ax7,'on')
xlabel(ax7,'Time (seconds)')
ylabel(ax7,'Aileron (%)')
if showlgd
    legend(ax7,lgdcell,'NumColumns',nt)
end
sgtitle(f,titstr)

grid(ax8,'on')
xlabel(ax8,'Time (seconds)')
ylabel(ax8,'Aircraft Heading (\circ)')
if showlgd
    legend(ax8,lgdcell,'NumColumns',nt)
end
title(ax8,titstr)

grid(ax9,'on')
xlabel(ax9,'East Position (ft)')
ylabel(ax9,'North Position (ft)')
if showlgd
    legend(ax9,lgdcell,'NumColumns',nt)
end
title(ax9,titstr)
axis(ax9,'equal')

grid(ax10,'on')
xlabel(ax10,'Time (seconds)')
ylabel(ax10,'Rudder (%)')
if showlgd
    legend(ax10,lgdcell,'NumColumns',nt)
end

grid(ax11,'on')
xlabel(ax11,'Time (seconds)')
ylabel(ax11,'Flap Stage')
if showlgd
    legend(ax11,lgdcell,'NumColumns',nt)
end

sgtitle(fi,titstr)
grid(axi1,'on')
grid(axi2,'on')
grid(axi3,'on')
grid(axi4,'on')
xlabel(axi1,'Time (seconds)')
xlabel(axi2,'Time (seconds)')
xlabel(axi3,'Time (seconds)')
xlabel(axi4,'Time (seconds)')
ylabel(axi1,'Altitude (ft)')
ylabel(axi2,'IAS (kn)')
ylabel(axi3,'Sideslip (\circ)')
ylabel(axi4,'Bank Angle (\circ)')
