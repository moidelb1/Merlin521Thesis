%Convert .asc datalog to .mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   +----------------------------------+                  %
%                   |  MERLIN MP521 FLIGHT SIMULATOR   |                  %
%                   |     Datalog Conversion Code      |                  %
%                   |   Pilot Tests Batch Conversion   |                  %
%                   |                                  |                  %
%                   |   Written By: Benjamin Moidel    |                  %
%                   |    Last Updated: Feb 25, 2021    |                  %
%                   +----------------------------------+                  %
%                                                                         %
%                              INSTRUCTIONS                               %
%   1) Fill in the model name as "modver" in line 28 of the code          %
%   1) Run datalogConvertPilot.m                                          %
%   2) Select the folder containing the subject's testing datalog files.  %
%       This folder should contain subfolders numbered by condition       %
%       (1,2,...,8).                                                      %
%       Example: ...\Pilot Test Results\C01 20201020                      %
%   3) Choose a destination directory to save the .mat files.             %
%   4) Verify the newly created files have the proper variables as listed %
%       in "Doc UG06_3.0 Excal datalog format.pdf"                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all

%Setup
modver = 'Piper_PA-28R-201_v10';  %model version name

warning('on')
warning('off','backtrace')  %deactivate warning line numbers

mods.Interpreter = 'tex';   %error dialog box options
mods.WindowStyle = 'modal';

%Recall last filepaths used
if exist('datalogConvertPilotlastpath.mat','file') == 2
    load('datalogConvertPilotlastpath.mat')     %load saved filepaths
end

if ~exist('lastlp','var')   %check for last load filepath in workspace
    lastlp = 'C:\';     %default load filepath to C drive
end

if ~exist('lastsp','var')   %check for last save filepath in workspace
    lastsp = 'C:\';     %default save filepath to C drive
end

%Test Results Directory (Ben's PC)
% trdir = ['D:\Storage\Documents\School Stuff\Grad School\Research\',...
%     'Pilot Test Results\C01 20201020'];     %hard-coded
trdir = uigetdir(lastlp,'Select Test Results Directory');   %user selected

if trdir == 0   %if directory selection cancelled
    errordlg('\fontsize{12}No folder selected.',...
        'Select Test Results Directory',mods);  %creates error dialog box
    error('Error loading data. No folder selected.')  %ends script early
end

%Selecting directory to save new datalog files
% savepath = ['D:\Storage\Documents\School Stuff\Grad School\Research\',...
%     'Pilot Test Results\Named Results'];    %hard-coded
savepath = uigetdir(lastsp,'Select Destination Folder');    %user selected

if savepath == 0    %if directory selection cancelled
    errordlg('\fontsize{12}No folder selected.',...
        'Select Destination Folder',mods);  %creates error dialog box
    error('Error saving data. No folder selected.')     %ends script early
end

%Locating datalog files to load
loadpath{8} = [];       %preallocate 8 condition folders
loadname{4,8} = [];     %preallocate 24 filenames (8 conditions, 3 trials)
dltime = nan(4,8);      %preallocate datalog timestamp matrix
for j = 1:8     %8 conditions
    loadpath{j} = [trdir,sprintf('\\%i',j)];
    files = dir(loadpath{j});
    tn = 1;     %reset trial number
    
    %Pick out standard datalog files
    for i = 1:length(files)     %check all files in folder
        if contains(files(i).name,'datalog.asc')
        	loadname{tn,j} = files(i).name; %collect datalog filenames
            dltime(tn,j) = str2double(loadname{tn,j}(9:14)); %timestamps
            tn = tn + 1;        %increment trial number
        end
    end
    %Sort file names by trial # using datalog timestamp
    [~,newind] = sort(dltime(:,j));
    loadname(:,j) = loadname(newind,j);     %reorder filename columns
end

%Generating new names for datalog workspace files
ind = regexp(trdir,'\');    %find index of all \ in filepath
ind = ind(end);     %index of final \ in filepath, precedes subject code
subjstr = trdir(ind+1:ind+3); %subject code string
savename{4,8} = []; %preallocate save name cell array
for j = 1:8     %8 conditions
    for i = 1:4     %at most 4 trials, usually 2 or 3
        savename{i,j} = sprintf('%s_%s-%i-%i',modver,subjstr,j,i);
    end
end

%Loading, formatting, and saving each file
ss = zeros(size(savename));    %preallocation
for j = 1:8     %8 conditions
    for i = 1:4     %at most 4 trials, usually 2 or 3
        if ~isnan(dltime(i,j))  %check if trial datalog exists
            ss(i,j) = loadnsave(loadpath{j},loadname{i,j},savepath,...
                savename{i,j},1,0);     %returns 1 if save successful
            if i == 4   %if a condition had 4 trials
                warning(['There is a trial 4 of condition %i ',...
                    'for subject %s'],j,subjstr)
            end
        elseif i <= 3   %if less than 3 datalogs were found per condition
            warning(['No datalog found for trial %i ',...
                'of condition %i for subject %s'],i,j,subjstr)
        end
    end
end

%Confirm save was successful
ss = logical(ss); %convert to logical vector
msg = [{'The following files were created:'};
    reshape(savename(ss),[],1)];
msgbox(msg,'Save Successful','help')

%Save last used load/save directories
[~,ind] = regexp(trdir,'\');    %find all \ in selected load directory
lastlp = trdir(1:ind(end)-1);   %get directory above trdir
[~,ind] = regexp(savepath,'\'); %find all \ in selected save directory
lastsp = savepath(1:ind(end)-1);    %get directory above savepath
save('datalogConvertPilotlastpath','lastlp','lastsp')


function result = loadnsave(lpath,lname,spath,sname,stype,printtime)
%Loads a .asc datalog file, converts it to the specified file type, and
%saves the new file.
%Inputs:
%   lpath = load file path as a char vector or string
%   lname = load file name as a char vector or string (with file extension)
%   spath = save file path as a char vector or string
%   sname = save file name as a char vector or string (without file
%           extension)
%   stype = file type index for saving
%           1 = .mat
%           2 = .asc
%   `       3 = .txt
%           4 = .dat
%           5 = .xlsx
%           6 = .xlsm
%           7 = .xls
%           8 = .csv
%   printtime = logical (0 or 1); if 1, will print load and save times to
%               command window
%Outputs:
%   result = diagnostic integer
%           0 = save unsuccessful
%           1 = save successful

result = 0; %set as 0 in case function encounters error and ends

%Sets datalog type index based on predicted datalog type
tic     %start load timer
if contains(lname,'aeroelastic')     %datalog aeroelastic
    dltype = 1;
    sizeM = [32 inf];
elseif contains(lname,'controls')    %datalog controls
    dltype = 2;
    sizeM = [128 inf];
elseif contains(lname,'datalog')     %datalog
    dltype = 3;
    sizeM = [75 inf];
else
    dltype = 0;
end

%Read data from plain text datalog file
fileId = fopen(fullfile(lpath,lname));    %get file ID number
A = fscanf(fileId,'%f',sizeM)';     %read number data into matrix form

if stype == 1   %create structure variable containing data
    switch dltype
        case 1      %datalog aeroelastic
            data.time = A(:,1);     %Simulation elapsed time [seconds]
            data.event = A(:,2);    %Event marker, zero for no event being 
                                    %recorded, incremental positive integer
                                    %to flag events
            data.H01incid = A(:,3);     %HPanel01 incidence increase [deg]
            data.H01dihed = A(:,4);     %HPanel01 dihedral increase [deg]
            data.H01zdisp = A(:,5);     %HPanel01 downward displacement [m]
            
            data.H02incid = A(:,6);     %HPanel02 incidence increase [deg]
            data.H02dihed = A(:,7);     %HPanel02 dihedral increase [deg]
            data.H02zdisp = A(:,8);     %HPanel02 downward displacement [m]
            
            data.H03incid = A(:,9);     %HPanel03 incidence increase [deg]
            data.H03dihed = A(:,10);    %HPanel03 dihedral increase [deg]
            data.H03zdisp = A(:,11);    %HPanel03 downward displacement [m]
            
            data.H04incid = A(:,12);    %HPanel04 incidence increase [deg]
            data.H04dihed = A(:,13);    %HPanel04 dihedral increase [deg]
            data.H04zdisp = A(:,14);    %HPanel04 downward displacement [m]
            
            data.H05incid = A(:,15);    %HPanel05 incidence increase [deg]
            data.H05dihed = A(:,16);    %HPanel05 dihedral increase [deg]
            data.H05zdisp = A(:,17);    %HPanel05 downward displacement [m]
            
            data.H06incid = A(:,18);    %HPanel06 incidence increase [deg]
            data.H06dihed = A(:,19);    %HPanel06 dihedral increase [deg]
            data.H06zdisp = A(:,20);    %HPanel06 downward displacement [m]
            
            data.H07incid = A(:,21);    %HPanel07 incidence increase [deg]
            data.H07dihed = A(:,22);    %HPanel07 dihedral increase [deg]
            data.H07zdisp = A(:,23);    %HPanel07 downward displacement [m]
            
            data.H08incid = A(:,24);    %HPanel08 incidence increase [deg]
            data.H08dihed = A(:,25);    %HPanel08 dihedral increase [deg]
            data.H08zdisp = A(:,26);    %HPanel08 downward displacement [m]
            
            data.H09incid = A(:,27);    %HPanel09 incidence increase [deg]
            data.H09dihed = A(:,28);    %HPanel09 dihedral increase [deg]
            data.H09zdisp = A(:,29);    %HPanel09 downward displacement [m]
            
            data.H10incid = A(:,30);    %HPanel10 incidence increase [deg]
            data.H10dihed = A(:,31);    %HPanel10 dihedral increase [deg]
            data.H10zdisp = A(:,32);    %HPanel10 downward displacement [m]
        case 2      %datalog controls
            data.time = A(:,1);     %Simulation elapsed time [seconds]
            data.event = A(:,2);    %Event marker, zero for no event being 
            
            data.H01teds = A(:,3);      %HPanel01 starboard TED [-1 to 1]
            data.H01tedp = A(:,4);      %HPanel01 port TED [-1 to 1]
            data.H01tedds = A(:,5);     %HPanel01 starboard TED [deg]
            data.H01teddp = A(:,6);     %HPanel01 port TED [deg]
            data.H01leds = A(:,7);      %HPanel01 starboard LED [-1 to 1]
            data.H01ledp = A(:,8);      %HPanel01 port LED [-1 to 1]
            data.H01usds = A(:,9);      %HPanel01 starboard USD [-1 to 1]
            data.H01usdp = A(:,10);     %HPanel01 port USD [-1 to 1]
            data.H01alps = A(:,11);     %HPanel01 starboard alpha [-1 to 1]
            data.H01alpp = A(:,12);     %HPanel01 port alpha [-1 to 1]
            data.H01swp = A(:,13);      %HPanel01 sweep angle [deg]
            
            data.H02teds = A(:,14);     %HPanel02 starboard TED [-1 to 1]
            data.H02tedp = A(:,15);     %HPanel02 port TED [-1 to 1]
            data.H02tedds = A(:,16);    %HPanel02 starboard TED [deg]
            data.H02teddp = A(:,17);    %HPanel02 port TED [deg]
            data.H02leds = A(:,18);     %HPanel02 starboard LED [-1 to 1]
            data.H02ledp = A(:,19);     %HPanel02 port LED [-1 to 1]
            data.H02usds = A(:,20);     %HPanel02 starboard USD [-1 to 1]
            data.H02usdp = A(:,21);     %HPanel02 port USD [-1 to 1]
            data.H02alps = A(:,22);     %HPanel02 starboard alpha [-1 to 1]
            data.H02alpp = A(:,23);     %HPanel02 port alpha [-1 to 1]
            data.H02swp = A(:,24);      %HPanel02 sweep angle [deg]
            
            data.H03teds = A(:,25);     %HPanel03 starboard TED [-1 to 1]
            data.H03tedp = A(:,26);     %HPanel03 port TED [-1 to 1]
            data.H03tedds = A(:,27);    %HPanel03 starboard TED [deg]
            data.H03teddp = A(:,28);    %HPanel03 port TED [deg]
            data.H03leds = A(:,29);     %HPanel03 starboard LED [-1 to 1]
            data.H03ledp = A(:,30);     %HPanel03 port LED [-1 to 1]
            data.H03usds = A(:,31);     %HPanel03 starboard USD [-1 to 1]
            data.H03usdp = A(:,32);     %HPanel03 port USD [-1 to 1]
            data.H03alps = A(:,33);     %HPanel03 starboard alpha [-1 to 1]
            data.H03alpp = A(:,34);     %HPanel03 port alpha [-1 to 1]
            data.H03swp = A(:,35);      %HPanel03 sweep angle [deg]
            
            data.H04teds = A(:,36);     %HPanel04 starboard TED [-1 to 1]
            data.H04tedp = A(:,37);     %HPanel04 port TED [-1 to 1]
            data.H04tedds = A(:,38);    %HPanel04 starboard TED [deg]
            data.H04teddp = A(:,39);    %HPanel04 port TED [deg]
            data.H04leds = A(:,40);     %HPanel04 starboard LED [-1 to 1]
            data.H04ledp = A(:,41);     %HPanel04 port LED [-1 to 1]
            data.H04usds = A(:,42);     %HPanel04 starboard USD [-1 to 1]
            data.H04usdp = A(:,43);     %HPanel04 port USD [-1 to 1]
            data.H04alps = A(:,44);     %HPanel04 starboard alpha [-1 to 1]
            data.H04alpp = A(:,45);     %HPanel04 port alpha [-1 to 1]
            data.H04swp = A(:,46);      %HPanel04 sweep angle [deg]
            
            data.H05teds = A(:,47);     %HPanel05 starboard TED [-1 to 1]
            data.H05tedp = A(:,48);     %HPanel05 port TED [-1 to 1]
            data.H05tedds = A(:,49);    %HPanel05 starboard TED [deg]
            data.H05teddp = A(:,50);    %HPanel05 port TED [deg]
            data.H05leds = A(:,51);     %HPanel05 starboard LED [-1 to 1]
            data.H05ledp = A(:,52);     %HPanel05 port LED [-1 to 1]
            data.H05usds = A(:,53);     %HPanel05 starboard USD [-1 to 1]
            data.H05usdp = A(:,54);     %HPanel05 port USD [-1 to 1]
            data.H05alps = A(:,55);     %HPanel05 starboard alpha [-1 to 1]
            data.H05alpp = A(:,56);     %HPanel05 port alpha [-1 to 1]
            data.H05swp = A(:,57);      %HPanel05 sweep angle [deg]
            
            data.H06teds = A(:,58);     %HPanel06 starboard TED [-1 to 1]
            data.H06tedp = A(:,59);     %HPanel06 port TED [-1 to 1]
            data.H06tedds = A(:,60);    %HPanel06 starboard TED [deg]
            data.H06teddp = A(:,61);    %HPanel06 port TED [deg]
            data.H06leds = A(:,62);     %HPanel06 starboard LED [-1 to 1]
            data.H06ledp = A(:,63);     %HPanel06 port LED [-1 to 1]
            data.H06usds = A(:,64);     %HPanel06 starboard USD [-1 to 1]
            data.H06usdp = A(:,65);     %HPanel06 port USD [-1 to 1]
            data.H06alps = A(:,66);     %HPanel06 starboard alpha [-1 to 1]
            data.H06alpp = A(:,67);     %HPanel06 port alpha [-1 to 1]
            data.H06swp = A(:,68);      %HPanel06 sweep angle [deg]
            
            data.H07teds = A(:,69);     %HPanel07 starboard TED [-1 to 1]
            data.H07tedp = A(:,70);     %HPanel07 port TED [-1 to 1]
            data.H07tedds = A(:,71);    %HPanel07 starboard TED [deg]
            data.H07teddp = A(:,72);    %HPanel07 port TED [deg]
            data.H07leds = A(:,73);     %HPanel07 starboard LED [-1 to 1]
            data.H07ledp = A(:,74);     %HPanel07 port LED [-1 to 1]
            data.H07usds = A(:,75);     %HPanel07 starboard USD [-1 to 1]
            data.H07usdp = A(:,76);     %HPanel07 port USD [-1 to 1]
            data.H07alps = A(:,77);     %HPanel07 starboard alpha [-1 to 1]
            data.H07alpp = A(:,78);     %HPanel07 port alpha [-1 to 1]
            data.H07swp = A(:,79);      %HPanel07 sweep angle [deg]
            
            data.H08teds = A(:,80);     %HPanel08 starboard TED [-1 to 1]
            data.H08tedp = A(:,81);     %HPanel08 port TED [-1 to 1]
            data.H08tedds = A(:,82);    %HPanel08 starboard TED [deg]
            data.H08teddp = A(:,83);    %HPanel08 port TED [deg]
            data.H08leds = A(:,84);     %HPanel08 starboard LED [-1 to 1]
            data.H08ledp = A(:,85);     %HPanel08 port LED [-1 to 1]
            data.H08usds = A(:,86);     %HPanel08 starboard USD [-1 to 1]
            data.H08usdp = A(:,87);     %HPanel08 port USD [-1 to 1]
            data.H08alps = A(:,88);     %HPanel08 starboard alpha [-1 to 1]
            data.H08alpp = A(:,89);     %HPanel08 port alpha [-1 to 1]
            data.H08swp = A(:,90);      %HPanel08 sweep angle [deg]
            
            data.H09teds = A(:,91);     %HPanel09 starboard TED [-1 to 1]
            data.H09tedp = A(:,92);     %HPanel09 port TED [-1 to 1]
            data.H09tedds = A(:,93);    %HPanel09 starboard TED [deg]
            data.H09teddp = A(:,94);    %HPanel09 port TED [deg]
            data.H09leds = A(:,95);     %HPanel09 starboard LED [-1 to 1]
            data.H09ledp = A(:,96);     %HPanel09 port LED [-1 to 1]
            data.H09usds = A(:,97);     %HPanel09 starboard USD [-1 to 1]
            data.H09usdp = A(:,98);     %HPanel09 port USD [-1 to 1]
            data.H09alps = A(:,99);     %HPanel09 starboard alpha [-1 to 1]
            data.H09alpp = A(:,100);    %HPanel09 port alpha [-1 to 1]
            data.H09swp = A(:,101);     %HPanel09 sweep angle [deg]
            
            data.H10teds = A(:,102);    %HPanel10 starboard TED [-1 to 1]
            data.H10tedp = A(:,103);    %HPanel10 port TED [-1 to 1]
            data.H10tedds = A(:,104);   %HPanel10 starboard TED [deg]
            data.H10teddp = A(:,105);   %HPanel10 port TED [deg]
            data.H10leds = A(:,106);    %HPanel10 starboard LED [-1 to 1]
            data.H10ledp = A(:,107);    %HPanel10 port LED [-1 to 1]
            data.H10usds = A(:,108);    %HPanel10 starboard USD [-1 to 1]
            data.H10usdp = A(:,109);    %HPanel10 port USD [-1 to 1]
            data.H10alps = A(:,110);    %HPanel10 starboard alpha [-1 to 1]
            data.H10alpp = A(:,111);    %HPanel10 port alpha [-1 to 1]
            data.H10swp = A(:,112);     %HPanel10 sweep angle [deg]
            
            data.V01ted = A(:,113);     %VPanel01 TED deflection [-1 to 1]
            data.V01beta = A(:,114);    %VPanel01 beta deflection [-1 to 1]
            data.V01tedd = A(:,115);    %VPanel01 TED deflection [deg]
            
            data.V02ted = A(:,116);     %VPanel02 TED deflection [-1 to 1]
            data.V02beta = A(:,117);    %VPanel02 beta deflection [-1 to 1]
            data.V02tedd = A(:,118);    %VPanel02 TED deflection [deg]
            
            data.V03ted = A(:,119);     %VPanel03 TED deflection [-1 to 1]
            data.V03beta = A(:,120);    %VPanel03 beta deflection [-1 to 1]
            data.V03tedd = A(:,121);    %VPanel03 TED deflection [deg]
            
            data.V04ted = A(:,122);     %VPanel04 TED deflection [-1 to 1]
            data.V04beta = A(:,123);    %VPanel04 beta deflection [-1 to 1]
            data.V04tedd = A(:,124);    %VPanel04 TED deflection [deg]
            
            data.thrust1 = A(:,125);    %Engine 1 thrust [N]
            data.thrust2 = A(:,126);    %Engine 2 thrust [N]
            data.thrust3 = A(:,127);    %Engine 3 thrust [N]
            data.thrust4 = A(:,128);    %Engine 4 thrust [N]
        case 3      %datalog (standard)
            data.time = A(:,1);     %Simulation elapsed time [seconds]
            data.event = A(:,2);    %Event marker, zero for no event being 
                                    %recorded, incremental positive integer
                                    %to flag events
            data.north = A(:,3);    %North position (from simulation start)
                                    %[m]
            data.east = A(:,4);     %East position (from simulation start) 
                                    %[m]
            data.alt = A(:,5);      %Altitude above sea level [m]
            data.lambda = A(:,6);   %Latitude [deg]
            data.mu = A(:,7);       %Longitude [deg]
            data.u = A(:,8);        %Body axis forward speed [m/s]
            data.v = A(:,9);        %Body axis lateral speed [m/s]
            data.w = A(:,10);       %Body axis vertical (down) speed [m/s]
            data.ax = A(:,11);      %Body axis forward acceleration [m/s^2]
            data.ay = A(:,12);      %Body axis lateral acceleration [m/s^2]
            data.az = A(:,13);      %Body axis normal acceleration [m/s^2]
            data.phi = A(:,14);     %Euler roll attitude [deg]
            data.theta = A(:,15);   %Euler pitch attitude [deg]
            data.psi = A(:,16);     %Euler heading angle [deg]
            data.p = A(:,17);       %Body axis roll rate [deg/s]
            data.q = A(:,18);       %Body axis pitch rate [deg/s]
            data.r = A(:,19);       %Body axis yaw rate [deg/s]
            data.pdot = A(:,20);    %Body axis roll acceleration [deg/s^2]
            data.qdot = A(:,21);    %Body axis pitch acceleration [deg/s^2]
            data.rdot = A(:,22);    %Body axis yaw acceleration [deg/s^2]
            data.alpha = A(:,23);   %Incidence angle [deg]
            data.beta = A(:,24);    %Sideslip angle [deg]
            data.gamma = A(:,25);   %Flight path angle [deg]
            data.V = A(:,26);       %True (inertial) speed [kn]
            data.VTAS = A(:,27);    %True air speed [kn]
            data.VIAS = A(:,28);    %Indicated air speed [kn]
            data.VEAS = A(:,29);    %Equivalent air speed [kn]
            data.VGS = A(:,30);     %Ground speed [kn]
            data.track = A(:,31);   %Ground track angle [deg]
            data.hdot = A(:,32);    %Rate of change of altitude [m/s]
            data.mach = A(:,33);    %Mach number
            data.nz = A(:,34);      %Load factor (g)
            data.h = A(:,35);       %Height above terrain (radio altitude) 
                                    %[m]
            data.eta = A(:,36);     %Elevator deflection [-1 to 1]
            data.xsi = A(:,37);     %Aileron deflection [-1 to 1]
            data.zeta = A(:,38);    %Rudder deflection [-1 to 1]
            data.d_ab = A(:,39);    %Airbrake (spoiler) deflection [0-1]
            data.d_flap = A(:,40);  %Flap deflection [0-n stage]
            data.d_gear = A(:,41);  %Undercarriage deflection [0-1]
            data.d_spoiler = A(:,42);%Roll spoiler deflection [deg]
            data.lth = A(:,43);     %Left engine thrust [N]
            data.rth = A(:,44);     %Right engine thrust [N]
            data.y_stick = A(:,45); %Longitudinal stick position [-1 to 1]
            data.x_stick = A(:,46); %Lateral stick position [-1 to 1]
            data.pedal = A(:,47);   %Pedal position [-1 to 1]
            data.l_throt = A(:,48); %Left throttle position [0-1]
            data.r_throt = A(:,49); %Right throttle position [0-1]
            data.sw_flap = A(:,50); %Flap switch position [0-n stage]
            data.sw_gear = A(:,51); %Gear switch position [0-1]
            data.sw_ab = A(:,52);   %Airbrake (spoiler) switch position 
                                    %[0-1]
            data.mass = A(:,53);    %Aircraft total mass [kg]
            data.Ixx = A(:,54);     %Moment of inertia [kg*m^2]
            data.Iyy = A(:,55);     %Moment of inertia [kg*m^2]
            data.Izz = A(:,56);     %Moment of inertia [kg*m^2]
            data.Ixz = A(:,57);     %Moment of inertia [kg*m^2]
            data.fuel = A(:,58);    %Fuel state [full 0-1 empty]
            data.cg_x = A(:,59);    %Current CG position [m]
            data.cg_y = A(:,60);    %Current CG position [m]
            data.cg_z = A(:,61);    %Current CG position [m]
            data.l_act = A(:,62);   %Left motion actuator [0-1]
            data.r_act = A(:,63);   %Right motion actuator [0-1]
            data.T = A(:,64);       %Air temperature [K]
            data.P = A(:,65);       %Air pressure [Pa]
            data.rho = A(:,66);     %Air density [kg/m^3]
            data.temp_ratio = A(:,67);  %Air temperature ratio
            data.delta = A(:,68);   %Air pressure ratio
            data.sigma = A(:,69);   %Air density ratio
            data.X = A(:,70);       %Total Body Axis Force [N]
            data.Y = A(:,71);       %Total Body Axis Force [N]
            data.Z = A(:,72);       %Total Body Axis Force [N]
            data.L = A(:,73);       %Total Moment [N*m]
            data.M = A(:,74);       %Total Moment [N*m]
            data.N = A(:,75);       %Total Moment [N*m]
        otherwise
            msg = {['\fontsize{12}Invalid file selected.',...
                'Select a *.asc file ending'];
                'with "datalog aeroelastic", "datalog controls",';
                'or "datalog".'};
            f = errordlg(msg,'Select File',mods);
            error('Invalid file selected.')
    end
    
else    %create cell array containing data and headings
    %Determine size of data matrix and proper headings
    switch dltype
        case 1      %datalog aeroelastic
            headings = cell(32,1);
            headings(1:2) = {'Time (sec)';'Event'};
            for i = 0:9
                str = sprintf('HPanel%02i',i+1);
                headings(3+3*i) = {[str,' Alpha/Twist (',char(176),')']};
                headings(4+3*i) = {[str,' Dihedral (',char(176),')']};
                headings(5+3*i) = {[str,' Z-Displacement (m)']};
            end
        case 2      %datalog controls
            headings = cell(128,1);
            headings(1:2) = {'Time (sec)';'Event'};
            for i = 0:9
                str = sprintf('HPanel%02i',i+1);
                headings(3+11*i) = {[str,' TED Starboard/Right (0-1)']};
                headings(4+11*i) = {[str,' TED Port/Left (0-1)']};
                headings(5+11*i) = {[str,' TED Deflection Starboard (',...
                    char(176),')']};
                headings(6+11*i) = {[str,' TED Deflection Port (',...
                    char(176),')']};
                headings(7+11*i) = {[str,' LED Starboard (0-1)']};
                headings(8+11*i) = {[str,' LED Port (0-1)']};
                headings(9+11*i) = {[str,' USD Starboard (0-1)']};
                headings(10+11*i) = {[str,' USD Port (0-1)']};
                headings(11+11*i) = {[str,' Alpha Starboard (0-1)']};
                headings(12+11*i) = {[str,' Alpha Port (0-1)']};
                headings(13+11*i) = {[str,' Sweep (',char(176),')']};
            end
            for i = 0:3
                str = sprintf('VPanel%02i',i+1);
                headings(113+3*i) = {[str,' TED (0-1)']};
                headings(114+3*i) = {[str,' Beta (0-1)']};
                headings(115+3*i) = {[str,' TED Deflection (',char(176),...
                    ')']};
            end
            for i = 0:3
                str = sprintf('Engine %i',i+1);
                headings(125+i) = {[str,' Thrust (N)']};
            end
        case 3      %datalog (standard)
            headings = {'Time (sec)';
                'Event';
                'North Position (m)';
                'East Position (m)';
                'Altitude (m)';
                ['Latitude (',char(176),')'];
                ['Longitude (',char(176),')'];
                'Forward Speed (m/s)';
                'Lateral Speed (m/s)';
                'Vertical Speed (m/s)';
                ['Forward Accel. (m/s',char(178),')'];
                ['Lateral Accel. (m/s',char(178),')'];
                ['Normal Accel. (m/s',char(178),')'];
                ['Roll Attitude (',char(176),')'];
                ['Pitch Attitude (',char(176),')'];
                ['Heading Angle (',char(176),')'];
                ['Roll Rate (',char(176),'/s)'];
                ['Pitch Rate (',char(176),'/s)'];
                ['Yaw Rate (',char(176),'/s)'];
                ['Roll Accel. (',char(176),'/s',char(178),')'];
                ['Pitch Accel. (',char(176),'/s',char(178),')'];
                ['Yaw Accel. (',char(176),'/s',char(178),')'];
                ['AoA (',char(176),')'];
                ['Sideslip Angle (',char(176),')'];
                ['Flight Path Angle (',char(176),')'];
                'True Inertial Speed (kn)';
                'True Air Speed (kn)';
                'Indicated Air Speed (kn)';
                'Equivalent Air Speed (kn)';
                'Ground Speed (kn)';
                ['Ground Track Angle (',char(176),')'];
                'Altitude Rate of Change (m/s)';
                'Mach';
                'Load Factor (g)';
                'Height Above Terrain/Radio Altitude (m)';
                'Elevator Deflection (+/-1)';
                'Aileron Deflection (+/-1)';
                'Rudder Deflection (+/-1)';
                'Airbrake Deflection (0-1)';
                'Flap Deflection (0-stage #)';
                'Undercarriage/Gear Deflection (0-1)';
                ['Roll Spoiler Deflection (',char(176),')'];
                'Left Engine Thrust (N)';
                'Right Engine Thrust (N)';
                'Longitudinal Stick Position (+/-1)';
                'Lateral Stick Position (+/-1)';
                'Pedal Position (+/-1)';
                'Left Throttle Position (0-1)';
                'Right Throttle Position (0-1)';
                'Flap Switch Position (0-1)';
                'Gear Switch Position (0-1)';
                'Airbrake Switch Position (0-1)';
                'Aircraft Mass (kg)';
                ['Ixx (kg',char(183),'m',char(178),')'];
                ['Iyy (kg',char(183),'m',char(178),')'];
                ['Izz (kg',char(183),'m',char(178),')'];
                ['Ixz (kg',char(183),'m',char(178),')'];
                'Fuel State (0-1)';
                'CG x-position (m)';
                'CG y-position (m)';
                'CG z-position (m)';
                'Left Motion Actuator (0-1)';
                'Right Motion Actuator (0-1)';
                'Air Temperature (K)';
                'Air Pressure (Pa)';
                ['Air Density (kg/m',char(179),')'];
                'Air Temperature Ratio';
                'Air Pressure Ratio';
                'Air Density Ratio';
                'X Body Force (N)';
                'Y Body Force (N)';
                'Z Body Force (N)';
                'Roll Moment (Nm)';
                'Pitch Moment (Nm)';
                'Yaw Moment (Nm)'};
    end
    datacell = [headings';num2cell(A)]; %combine headings and data
end

if printtime
    t1 = toc;   %end load timer
    fprintf('Time taken to load data: %.4f seconds\n',t1) %print load timer
end

%Save file
tic     %start save timer
switch stype
    case 1  %save as workspace
        save(fullfile(spath,sname),'-struct','data')
    case 2  %save as ASCII .asc file
        writecell(datacell,fullfile(spath,sname),'FileType','text')
    case {3,4}  %save as delimited text or .dat file
        writecell(datacell,fullfile(spath,sname),'Delimiter','space')
    case {5,6,7}    %save as spreadsheet .xlsx, .xlsm, or .xls file
        writecell(datacell,fullfile(spath,sname),'Sheet','datalog')
    case 8  %save as .csv file
        writecell(datacell,fullfile(spath,sname))
end

if printtime
    t2 = toc;   %end save timer
    fprintf('Time taken to save data: %.4f seconds\n',t2) %print save timer
end

result = 1;     %indicates successful save
end