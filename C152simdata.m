%Benjamin Moidel
%November 23, 2020
%Pilot Test Result Plotting - Cessna 152

clc;clear;close all

%Filepath to pilot datalog .mat files
fpath = ['D:\Storage\Documents\School Stuff\Grad School',...
         '\Research\Pilot Test Results\MATLAB Datalogs'];

%Names of each pilot folder (add to this list as more are created)
% -UNUSED-
% pilotpath = {'C01 20201020';
%              'C02 20201022';
%              'C03 20201026';
%              'C04 20201028';
%              'C05 20201102';
%              'C06 20201103';
%              'C07 20201105';
%              'C08 20201109'};
% -UNUSED-

%Loading Datalog Parameters
pilot = 1:9;    %vector list of pilots to load (1st dim)
task = 1:8;   %vector list of tasks to load (2nd dim)
trial = 1:3;    %vector list of trials to load (3rd dim)

imax = length(pilot);   %size of 1st dim of data structure
jmax = length(task);    %size of 2nd dim of data structure
kmax = length(trial);   %size of 3rd dim of data structure
nfiles = imax*jmax*kmax;    %number of datalog files to load

fprintf('Loading data...\n')
fprintf('Number of files:\t%i\n',nfiles)    %display file count
prog = round(linspace(1,nfiles,21));    %progress updates at 5% intervals

tic     %start load timer
s(imax,jmax,kmax) = struct(datadummy());    %preallocate structure
f = 0;  %file counter
pbar = '';   %preallocate text progress bar
for i = 1:imax
    for j = 1:jmax
        for k = 1:kmax
            %Progress Update
            c = k+(j-1)*kmax+(i-1)*jmax*kmax;   %iteration counter
            if any(prog==c)     %when iteration is a multiple of 5%...
                for z=1:20
                    if prog(z) < c  %...print "-" until that %...
                        pbar(z) = '-';
                    else            %...then print " " until end of bar
                        pbar(z) = ' ';
                    end
                end
                home    %align next line at top of command window
                fprintf('|%s%3i%%%s|\n',pbar(1:10),find(prog==c,1)*5-5,...
                    pbar(11:20))    %print progress bar
            end
            
            %Creating filename
            if pilot(i) <= 3
                modver = 18;    %C01-C03 used v18
            else
                modver = 19;    %C04-C08 used v19
            end
            fname = sprintf('Cessna_152_v%02i_C%02i-%i-%i.mat',modver,...
                pilot(i),task(j),trial(k));
            fload = fullfile(fpath,fname);  %merge filename and filepath
            if isfile(fload)
                s(i,j,k) = load(fload);
                f=f+1;
            else
                warning('File cannot be found: %s\n',fname)
            end
        end
    end
end
t = toc;    %end load timer
fprintf('Load complete. %i of %i files found and loaded in %.3f seconds.\n',...
    f,i*j*k,t);

save C152data s imax jmax kmax pilot task trial

%% Plotting
close all
clc;clear
load C152data

%Selecting Data
pilots = 4:9;
conds = [1 6];
trials = 1:3;

imax = length(pilots);
jmax = length(conds);
kmax = length(trials);

%Color Selection
cpi = interp1(1:64,hsv(64),linspace(1,60,imax));    %pilot colors
ctr = linspace(0.25,1,kmax);    %brightness multiplier for trial #

%Plotting flight path
lgdcell{imax*kmax} = '';
for j = 1:jmax  %iterate through conditions
    figure
    hold on
    grid on
    view(3)     %change to 3D view
    c = conds(j);   %abbreviated condition index
    for k = 1:kmax  %iterate through trials
        t = trials(k);  %abbreviated trial index
        for i = 1:imax  %iterate through pilots
            p = pilots(i);  %abbreviated pilot index
            if(isempty(s(p,c,t).time))  %check for absent datalogs
                x0 = s(1,c,1).east(1);
                y0 = s(1,c,1).north(1);
                z0 = 0;
                plot3(x0,y0,z0,'k--')   %plot dummy point
            else
                plot3(s(p,c,t).east,s(p,c,t).north,s(p,c,t).alt,...
                      'Color',cpi(i,:)*ctr(k))
            end
            lgdcell{i+(k-1)*imax} = sprintf('C%02i Trial %i',...
                                            p,t);
        end
    end
    daspect([1 1 1])
    zlim([0,inf])
    axis vis3d
    xlabel('East Position (m)')
    ylabel('North Position (m)')
    zlabel('Altitude (m)')
    legend(lgdcell,'NumColumns',kmax)
    title(sprintf('Condition #%i',c))
end

%Plotting Control Inputs
for j = 1:jmax
    figure
    c = conds(j);   %abbreviated condition index
    for k = 1:kmax
        t = trials(k);  %abbreviated trial index
        for i = 1:imax
            p = pilots(i);  %abbreviated pilot index
            if(isempty(s(p,c,t).time))  %check for absent datalogs
                subplot(4,1,1)
                hold on
                plot(0,0,'k--')     %plot dummy point
                subplot(4,1,2)
                hold on
                plot(0,0,'k--')     %plot dummy point
                subplot(4,1,3)
                hold on
                plot(0,0,'k--')     %plot dummy point
                subplot(4,1,4)
                hold on
                plot(0,0,'k--')     %plot dummy point
            else
                subplot(4,1,1)
                hold on
                plot(s(p,c,t).time,s(p,c,t).x_stick,...
                     'Color',cpi(i,:)*ctr(k))
                
                subplot(4,1,2)
                hold on
                plot(s(p,c,t).time,s(p,c,t).y_stick,...
                     'Color',cpi(i,:)*ctr(k))
                
                subplot(4,1,3)
                hold on
                plot(s(p,c,t).time,s(p,c,t).pedal,...
                     'Color',cpi(i,:)*ctr(k))
                
                subplot(4,1,4)
                hold on
                plot(s(p,c,t).time,s(p,c,t).l_throt,...
                     'Color',cpi(i,:)*ctr(k))
            end
        end
    end
    subplot(4,1,1)  %stick x-deflection (roll)
    grid on
%     xlabel('Time (in seconds)')
    ylabel('Stick X-Deflection, Roll')
    title(sprintf('Condition #%i',c))
    
    subplot(4,1,2)  %stick y-deflection (pitch)
    grid on
%     xlabel('Time (in seconds)')
    ylabel('Stick Y-Deflection, Pitch')
    
    subplot(4,1,3)  %pedal deflection (yaw)
    grid on
%     xlabel('Time (in seconds)')
    ylabel('Pedal Deflection')
    
    subplot(4,1,4)  %throttle setting
    grid on
    xlabel('Time (in seconds)')
    ylabel('Throttle Setting')
    legend(lgdcell,'NumColumns',kmax)
end

%Plotting Attitude Time History
for j = 1:jmax
    figure
    c = conds(j);   %abbreviated condition index
    for k = 1:kmax
        t = trials(k);  %abbreviated trial index
        for i = 1:imax
            p = pilots(i);  %abbreviated pilot index
            if(isempty(s(p,c,t).time))  %check for absent datalogs
                subplot(4,1,1)
                hold on
                plot(0,0,'k--')     %plot dummy point
                subplot(4,1,2)
                hold on
                plot(0,0,'k--')     %plot dummy point
                subplot(4,1,3)
                hold on
                plot(0,0,'k--')     %plot dummy point
                subplot(4,1,4)
                hold on
                plot(0,0,'k--')     %plot dummy point
            else
                subplot(3,1,1)
                hold on
                plot(s(p,c,t).time,s(p,c,t).phi,...
                     'Color',cpi(i,:)*ctr(k))
                
                subplot(3,1,2)
                hold on
                plot(s(p,c,t).time,s(p,c,t).theta,...
                     'Color',cpi(i,:)*ctr(k))
                
                subplot(3,1,3)
                hold on
                plot(s(p,c,t).time,s(p,c,t).psi,...
                     'Color',cpi(i,:)*ctr(k))
            end
        end
    end
    subplot(3,1,1)  %bank angle
    grid on
%     xlabel('Time (in seconds)')
    ylabel('Bank (\circ)')
    title(sprintf('Condition #%i',c))
    
    subplot(3,1,2)  %pitch angle
    grid on
%     xlabel('Time (in seconds)')
    ylabel('Pitch (\circ)')
    
    subplot(3,1,3)  %heading angle
    grid on
    xlabel('Time (in seconds)')
    ylabel('Heading (\circ)')
    legend(lgdcell,'NumColumns',kmax)
end

%Plotting Flap Deflection
for j = 1:jmax
    figure
    c = conds(j);   %abbreviated condition index
    for k = 1:kmax
        t = trials(k);  %abbreviated trial index
        for i = 1:imax
            p = pilots(i);  %abbreviated pilot index
            if(isempty(s(p,c,t).time))
                subplot(5,1,1)
                hold on
                plot(0,s(1,1,1).X(1),'k--')     %plot dummy point
                subplot(5,1,2)
                hold on
                plot(0,s(1,1,1).Z(1),'k--')     %plot dummy point
                subplot(5,1,3)
                hold on
                plot(0,0,'k--')     %plot dummy point
                subplot(5,1,4)
                hold on
                plot(0,s(1,1,1).VIAS(1),'k--')  %plot dummy point
                subplot(5,1,5)
                hold on
                plot(0,0,'k--')     %plot dummy point
            else
                subplot(5,1,1)
                hold on
                plot(s(p,c,t).time,s(p,c,t).X,...
                     'Color',cpi(i,:)*ctr(k))
                
                subplot(5,1,2)
                hold on
                plot(s(p,c,t).time,s(p,c,t).Z,...
                     'Color',cpi(i,:)*ctr(k))
                
                subplot(5,1,3)
                hold on
                plot(s(p,c,t).time,s(p,c,t).alpha,...
                     'Color',cpi(i,:)*ctr(k))
                
                subplot(5,1,4)
                hold on
                plot(s(p,c,t).time,s(p,c,t).VIAS,...
                     'Color',cpi(i,:)*ctr(k))
                
                subplot(5,1,5)
                hold on
                plot(s(p,c,t).time,s(p,c,t).d_flap,...
                     'Color',cpi(i,:)*ctr(k))
            end
        end
    end
    subplot(5,1,1)
    grid on
    ylabel('X-Body Force (in N)')
    title(sprintf('Condition #%i',c))
    
    subplot(5,1,2)
    grid on
    ylabel('Z-Body Force (in N)')
    
    subplot(5,1,3)
    grid on
    ylabel('Angle of Attack (in \circ)')
    
    subplot(5,1,4)
    grid on
    ylabel('Indicated Airspeed (in knots)')
    
    subplot(5,1,5)
    grid on
    xlabel('Time (in seconds)')
    ylabel('Flap Setting (stage #)')
    legend(lgdcell,'NumColumns',kmax)
end

%% Functions
function stemp = datadummy()
%Returns the template structure for standard EIII datalog files
%Inputs:
%   none
%Outputs:
%   stemp = empty template structure for EIII datalog variables

stemp.time = [];    %Simulation elapsed time [seconds]
stemp.event = [];   %Event marker, zero for no event being recorded,
                        %incremental positive integer to flag events
stemp.north = [];   %North position (from simulation start) [m]
stemp.east = [];    %East position (from simulation start) [m]
stemp.alt = [];     %Altitude above sea level [m]
stemp.lambda = [];  %Latitude [deg]
stemp.mu = [];      %Longitude [deg]
stemp.u = [];       %Body axis forward speed [m/s]
stemp.v = [];       %Body axis lateral speed [m/s]
stemp.w = [];       %Body axis vertical (down) speed [m/s]
stemp.ax = [];      %Body axis forward acceleration [m/s^2]
stemp.ay = [];      %Body axis lateral acceleration [m/s^2]
stemp.az = [];      %Body axis normal acceleration [m/s^2]
stemp.phi = [];     %Euler roll attitude [deg]
stemp.theta = [];   %Euler pitch attitude [deg]
stemp.psi = [];     %Euler heading angle [deg]
stemp.p = [];       %Body axis roll rate [deg/s]
stemp.q = [];       %Body axis pitch rate [deg/s]
stemp.r = [];       %Body axis yaw rate [deg/s]
stemp.pdot = [];    %Body axis roll acceleration [deg/s^2]
stemp.qdot = [];    %Body axis pitch acceleration [deg/s^2]
stemp.rdot = [];    %Body axis yaw acceleration [deg/s^2]
stemp.alpha = [];   %Incidence angle [deg]
stemp.beta = [];    %Sideslip angle [deg]
stemp.gamma = [];   %Flight path angle [deg]
stemp.V = [];       %True (inertial) speed [kn]
stemp.VTAS = [];    %True air speed [kn]
stemp.VIAS = [];    %Indicated air speed [kn]
stemp.VEAS = [];    %Equivalent air speed [kn]
stemp.VGS = [];     %Ground speed [kn]
stemp.track = [];   %Ground track angle [deg]
stemp.hdot = [];    %Rate of change of altitude [m/s]
stemp.mach = [];    %Mach number
stemp.nz = [];      %Load factor (g)
stemp.h = [];       %Height above terrain (radio altitude) [m]
stemp.eta = [];     %Elevator deflection [-1 to 1]
stemp.xsi = [];     %Aileron deflection [-1 to 1]
stemp.zeta = [];    %Rudder deflection [-1 to 1]
stemp.d_ab = [];    %Airbrake (spoiler) deflection [0-1]
stemp.d_flap = [];  %Flap deflection [0-n stage]
stemp.d_gear = [];  %Undercarriage deflection [0-1]
stemp.d_spoiler = [];%Roll spoiler deflection [deg]
stemp.lth = [];     %Left engine thrust [N]
stemp.rth = [];     %Right engine thrust [N]
stemp.y_stick = []; %Longitudinal stick position [-1 to 1]
stemp.x_stick = []; %Lateral stick position [-1 to 1]
stemp.pedal = [];   %Pedal position [-1 to 1]
stemp.l_throt = []; %Left throttle position [0-1]
stemp.r_throt = []; %Right throttle position [0-1]
stemp.sw_flap = []; %Flap switch position [0-n stage]
stemp.sw_gear = []; %Gear switch position [0-1]
stemp.sw_ab = [];   %Airbrake (spoiler) switch position [0-1]
stemp.mass = [];    %Aircraft total mass [kg]
stemp.Ixx = [];     %Moment of inertia [kg*m^2]
stemp.Iyy = [];     %Moment of inertia [kg*m^2]
stemp.Izz = [];     %Moment of inertia [kg*m^2]
stemp.Ixz = [];     %Moment of inertia [kg*m^2]
stemp.fuel = [];    %Fuel state [full 0-1 empty]
stemp.cg_x = [];    %Current CG position [m]
stemp.cg_y = [];    %Current CG position [m]
stemp.cg_z = [];    %Current CG position [m]
stemp.l_act = [];   %Left motion actuator [0-1]
stemp.r_act = [];   %Right motion actuator [0-1]
stemp.T = [];       %Air temperature [K]
stemp.P = [];       %Air pressure [Pa]
stemp.rho = [];     %Air density [kg/m^3]
stemp.temp_ratio = [];  %Air temperature ratio
stemp.delta = [];   %Air pressure ratio
stemp.sigma = [];   %Air density ratio
stemp.X = [];       %Total Body Axis Force [N]
stemp.Y = [];       %Total Body Axis Force [N]
stemp.Z = [];       %Total Body Axis Force [N]
stemp.L = [];       %Total Moment [N*m]
stemp.M = [];       %Total Moment [N*m]
stemp.N = [];       %Total Moment [N*m]
% stemp.datacell = [];    %cell containing the previous data + headings
end