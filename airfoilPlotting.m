%Benjamin Moidel
%April 22, 2021
%C152/Piper Arrow III Airfoil Plots

%% Loading data from Excel
clc;clear;close all

fname = {'Cessna 152 Data.xlsx';
         'Piper PA-28R-201 Arrow Data.xlsx'};     %file names
fpath = 'D:\Storage\Documents\School Stuff\Grad School\Research';   %path
f = fullfile(fpath,fname);
sh = {'0012';'2412';'65-415'};  %sheet names

%Naming Airfoil Data
lgdC152title = 'Re = 3e6';
lgdC152 = {'NACA 0012';
           'NACA 0012 flap -18\circ';
           'NACA 0012 flap -9\circ';
           'NACA 0012 flap 11\circ';
           'NACA 0012 flap 12\circ';
           'NACA 0012 flap 23\circ';
           'NACA 0012 flap 25\circ';
           'NACA 2412';
           'NACA 2412 flap -19\circ';
           'NACA 2412 flap -10\circ';
           'NACA 2412 flap 8\circ';
           'NACA 2412 flap 10\circ';
           'NACA 2412 flap 16\circ';
           'NACA 2412 flap 20\circ';
           'NACA 2412 flap 30\circ'};
lgdPipertitle = 'Re = 4e6';
lgdPiper = {'NACA 0012';
            'NACA 0012 flap 14\circ';
            'NACA 0012 flap 28\circ';
            'NACA 65_2-415';
            'NACA 65_2-415 flap -25\circ';
            'NACA 65_2-415 flap -13\circ';
            'NACA 65_2-415 flap 6\circ';
            'NACA 65_2-415 flap 10\circ';
            'NACA 65_2-415 flap 12.5\circ';
            'NACA 65_2-415 flap 25\circ';
            'NACA 65_2-415 flap 40\circ'};


%Grab data for each condition run in XFLR5 (C152, Re = 3e6)
Mc = NaN(15,51,4);  %preallocate
Mc(1,:,:) = readmatrix(f{1},'Sheet',sh{1},'Range','A4:D54');    %NACA 0012
Mc(2,:,:) = readmatrix(f{1},'Sheet',sh{1},'Range','N4:Q54');    %flap -18
Mc(3,:,:) = readmatrix(f{1},'Sheet',sh{1},'Range','F4:I54');    %flap -9
Mc(4,:,:) = readmatrix(f{1},'Sheet',sh{1},'Range','V4:Y54');    %flap 11
Mc(5,:,:) = readmatrix(f{1},'Sheet',sh{1},'Range','AD4:AG54');  %flap 12
Mc(6,:,:) = readmatrix(f{1},'Sheet',sh{1},'Range','AL4:AO54');  %flap 23
Mc(7,:,:) = readmatrix(f{1},'Sheet',sh{1},'Range','AT4:AW54');  %flap 25
Mc(8,:,:) = readmatrix(f{1},'Sheet',sh{2},'Range','A4:D54');    %NACA 2412
Mc(9,:,:) = readmatrix(f{1},'Sheet',sh{2},'Range','F4:I54');    %flap -19
Mc(10,:,:) = readmatrix(f{1},'Sheet',sh{2},'Range','N4:Q54');   %flap -10
Mc(11,:,:) = readmatrix(f{1},'Sheet',sh{2},'Range','V4:Y54');   %flap 8
Mc(12,:,:) = readmatrix(f{1},'Sheet',sh{2},'Range','BB4:BE54'); %flap 10
Mc(13,:,:) = readmatrix(f{1},'Sheet',sh{2},'Range','AD4:AG54'); %flap 16
Mc(14,:,:) = readmatrix(f{1},'Sheet',sh{2},'Range','BR4:BU54'); %flap 20
Mc(15,:,:) = readmatrix(f{1},'Sheet',sh{2},'Range','BZ4:CC54'); %flap 30

%Grab data for each condition run in XFLR5 (Arrow III, Re = 4e6)
Mp = NaN(11,121,4); %preallocate
Mp(1,:,:) = readmatrix(f{2},'Sheet',sh{1},'Range','A4:D124');   %NACA 0012
Mp(2,:,:) = readmatrix(f{2},'Sheet',sh{1},'Range','AL4:AO124'); %flap 14
Mp(3,:,:) = readmatrix(f{2},'Sheet',sh{1},'Range','AT4:AW124'); %flap 28
Mp(4,:,:) = readmatrix(f{2},'Sheet',sh{3},'Range','A4:D124');   %NACA 65415
Mp(5,:,:) = readmatrix(f{2},'Sheet',sh{3},'Range','F4:I124');   %flap -25
Mp(6,:,:) = readmatrix(f{2},'Sheet',sh{3},'Range','N4:Q124');   %flap -13
Mp(7,:,:) = readmatrix(f{2},'Sheet',sh{3},'Range','V4:Y124');   %flap 6
Mp(8,:,:) = readmatrix(f{2},'Sheet',sh{3},'Range','AL4:AO124'); %flap 10
Mp(9,:,:) = readmatrix(f{2},'Sheet',sh{3},'Range','AD4:AG124'); %flap 12.5
Mp(10,:,:) = readmatrix(f{2},'Sheet',sh{3},'Range','AT4:AW124'); %flap 25
Mp(11,:,:) = readmatrix(f{2},'Sheet',sh{3},'Range','BB4:BE124'); %flap 40

%Angle of attack
AoAc = Mc(1,:,1);
AoAp = Mp(1,:,1);

%Lift coefficient
Clc = Mc(:,:,2);
Clp = Mp(:,:,2);

%Drag coefficient
Cdc = Mc(:,:,3);
Cdp = Mp(:,:,3);

%Pitching Moment coefficient
Cmc = Mc(:,:,4);
Cmp = Mp(:,:,4);

save airfoilData.mat AoAc AoAp Clc Clp Cdc Cdp Cmc Cmp lgdC152 lgdPiper...
     lgdC152title lgdPipertitle

%% Plotting
clc;clear;close all
load airfoilData.mat

n1 = 7; %number of NACA 0012 lines
n2 = length(lgdC152)-n1;    %number of NACA 2412 lines
c1 = hsv(n1-1)*0.9;     %assign line colors
c1 = [0 0 0;c1];        %make first line black
c2 = hsv(n2-1)*0.9;     %assign line colors
c2 = [0 0 0;c2];        %make first line black
c = [c1;c2];

figure
subplot(2,2,1)
hold on
grid on
for i = 1:n1
    plot(AoAc,Clc(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_l')
% xlabel('Angle of Attack (\circ)')
% ylabel('2-D Lift Coefficient')
% l = legend(lgdC152);
% title(l,lgdC152title)

subplot(2,2,2)
hold on
grid on
for i = 1:n1
    plot(AoAc,Cdc(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_d')
ylim([0 0.1])

subplot(2,2,3)
hold on
grid on
for i = 1:n1
    plot(AoAc,Cmc(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_{m_{c/4}}')

l = legend(lgdC152(1:n1));
title(l,lgdC152title)
pos = l.Position;
newx = 0.5+(0.5-pos(3))/2;
newy = (0.5-pos(4))/2;
l.Position = [newx,newy,pos(3:4)];

figure
subplot(2,2,1)
hold on
grid on
for i = n1+1:n1+n2
    plot(AoAc,Clc(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_l')
% xlabel('Angle of Attack (\circ)')
% ylabel('2-D Lift Coefficient')
% l = legend(lgdC152);
% title(l,lgdC152title)

subplot(2,2,2)
hold on
grid on
for i = n1+1:n1+n2
    plot(AoAc,Cdc(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_d')
ylim([0 0.1])

subplot(2,2,3)
hold on
grid on
for i = n1+1:n1+n2
    plot(AoAc,Cmc(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_{m_{c/4}}')

l = legend(lgdC152(n1+1:end));
title(l,lgdC152title)
pos = l.Position;
newx = 0.5+(0.5-pos(3))/2;
newy = (0.5-pos(4))/2;
l.Position = [newx,newy,pos(3:4)];


n3 = 3; %number of NACA 0012 lines
n4 = length(lgdPiper)-n3;    %number of NACA 2412 lines
c3 = hsv(n3-1)*0.9;     %assign line colors
c3 = [0 0 0;c3];        %make first line black
c4 = hsv(n4-1)*0.9;     %assign line colors
c4 = [0 0 0;c4];        %make first line black
c = [c3;c4];

figure
subplot(2,2,1)
hold on
grid on
for i = 1:n3
    plot(AoAp,Clp(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_l')

subplot(2,2,2)
hold on
grid on
for i = 1:n3
    plot(AoAp,Cdp(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_d')
ylim([0 0.1])

subplot(2,2,3)
hold on
grid on
for i = 1:n3
    plot(AoAp,Cmp(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_{m_{c/4}}')

l = legend(lgdPiper(1:n3));
title(l,lgdPipertitle)
pos = l.Position;
newx = 0.5+(0.5-pos(3))/2;
newy = (0.5-pos(4))/2;
l.Position = [newx,newy,pos(3:4)];

figure
subplot(2,2,1)
hold on
grid on
for i = n3+1:n3+n4
    plot(AoAp,Clp(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_l')

subplot(2,2,2)
hold on
grid on
for i = n3+1:n3+n4
    plot(AoAp,Cdp(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_d')
ylim([0 0.1])

subplot(2,2,3)
hold on
grid on
for i = n3+1:n3+n4
    plot(AoAp,Cmp(i,:),'Color',c(i,:))
end
xlabel('\alpha (\circ)')
ylabel('C_{m_{c/4}}')

l = legend(lgdPiper(n3+1:end));
title(l,lgdPipertitle)
pos = l.Position;
newx = 0.5+(0.5-pos(3))/2;
newy = (0.5-pos(4))/2;
l.Position = [newx,newy,pos(3:4)];