%Benjamin Moidel
%October 20, 2020
%Pilot Condition Randomization Matrix

clc;clear;close all

numtask = 8;    %number of tasks
numpilots = 10; %number of pilots

M = zeros(numtask,numpilots);
for i = 1:numpilots
    M(:,i) = randperm(numtask);     %randomize task order for each pilot
end

%Plotting average condition # for each testing position
Mavg = mean(M,2);
subplot(2,1,1);
bar(Mavg)
hold on
xlabel('Testing Position/Order')
ylabel('Average Condition #')

%Plotting average testing position for each condition #
tavg = zeros(numtask,1);
for i = 1:numtask
    [row,~] = find(M==i);
    tavg(i) = mean(row);
end

subplot(2,1,2)
bar(tavg)
xlabel('Condition #')
ylabel('Average Testing Position')