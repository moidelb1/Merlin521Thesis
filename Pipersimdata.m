%Benjamin Moidel
%March 16, 2021
%Pilot Test Results Structure Creation - Piper PA-28R-201 Arrow III

clc;clear;close all

%Filepath to pilot datalog .mat files
fpath = ['D:\Storage\Documents\School Stuff\Grad School',...
         '\Research\Pilot Test Results\MATLAB Datalogs'];

%Loading Datalog Parameters
pilot = 1:11;   %vector list of pilots to load (1st dim)
task = 1:8;     %vector list of tasks to load (2nd dim)
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
            switch pilot(i)
                case 1
                    modver = 4;     %P01 used v04
                case 2
                    modver = 5;     %P02 used v05
                case {3,4}
                    modver = 6;     %P03 & P04 used v06
                case 5
                    modver = 8;     %P05 used v08
                otherwise
                    modver = 10;    %P06+ use v10
            end
            fname = sprintf('Piper_PA-28R-201_v%02i_P%02i-%i-%i.mat',...
                modver,pilot(i),task(j),trial(k));
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

save Piperdata s imax jmax kmax pilot task trial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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