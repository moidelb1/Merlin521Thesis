%Benjamin Moidel
%April 20, 2020
%Pilot Task Evaluation Responses

%% Loading data from Excel
clc;clear;close all

fname = 'Task Evaluation Survey PA-28R-201 (Responses) 20210421.xlsx';
fpath = 'D:\Storage\Documents\School Stuff\Grad School\Research';
f = fullfile(fpath,fname);

%aircraft characteristics responses
Mac = readmatrix(f,'Sheet','Sorted Form Responses','Range','C31:M54');
%pilot effort responses
Mpe = readmatrix(f,'Sheet','Sorted Form Responses','Range','C58:M81');
%overall rating responses
Mor = readmatrix(f,'Sheet','Sorted Form Responses','Range','C85:M108');

save PiperResponses.mat Mac Mpe Mor

%% Analyzing and Plotting data
clc;clear;close all

load PiperResponses.mat

%Data dimensions
numPilots = 11;     %total # of subjects that flew models
numTasks = 8;       %total # of tasks/conditions each subject flew
numTrials = 3;      %total # of trials per task
p0 = 1:5;           %subject #'s that flew PA-28R-201 v9 and before
p10 = 6:numPilots;  %subjects that flew PA-28R-201 v10
maxAC = 7;  %maximum aircraft characteristics rating
maxPE = 9;  %maximum pilot effort rating
maxOR = 10; %maximum overall rating

alpha = 0.05;   %statistical significance (two-tailed)
N = 10000;      %number of Monte-Carlo cases for power calculations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% medAC = (maxAC+1)/2;    %median aircraft characteristics rating
% medPE = (maxPE+1)/2;    %median pilot effort rating
% medOR = (maxOR+1)/2;    %median overall rating
medAC = 4.5;    %5 is perceptibly realistic, 4 is not
medPE = 7.5;    %8 is perceptibly realistic, 7 is not
medOR = 7.5;    %8 is perceptibly realistic, 7 is not

%Reshaping Data Matrices so task responses are all in same ROW
Mac0 = reshape(Mac(:,p0)',length(p0)*numTrials,numTasks)';
Mpe0 = reshape(Mpe(:,p0)',length(p0)*numTrials,numTasks)';
Mor0 = reshape(Mor(:,p0)',length(p0)*numTrials,numTasks)';
Mac10 = reshape(Mac(:,p10)',length(p10)*numTrials,numTasks)';
Mpe10 = reshape(Mpe(:,p10)',length(p10)*numTrials,numTasks)';
Mor10 = reshape(Mor(:,p10)',length(p10)*numTrials,numTasks)';

%Preallocation
z = zeros(numTasks,1);  %dummy zero vector with the correct size

pAC10sr = z;    %p-value for v10 aircraft characteristics signed rank
hAC10sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powAC10sr = z;  %power of test

pPE10sr = z;    %p-value for v10 pilot effort signed rank
hPE10sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powPE10sr = z;  %power of test

pOR10sr = z;    %p-value for v10 overall rating signed rank
hOR10sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powOR10sr = z;  %power of test

pAC0sr = z;    %p-value for early aircraft characteristics signed rank
hAC0sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powAC0sr = z;  %power of test

pPE0sr = z;    %p-value for early effort signed rank
hPE0sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powPE0sr = z;  %power of test

pOR0sr = z;    %p-value for early overall rating signed rank
hOR0sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powOR0sr = z;  %power of test

pACmw = z;      %p-value for aircraft characteristics Mann-Whitney U-test
hACmw = z;      %1 = reject null hypothesis, 0 = fail to reject
powACmw = z;    %power of test

pPEmw = z;      %p-value for pilot effort Mann-Whitney U-test
hPEmw = z;      %1 = reject null hypothesis, 0 = fail to reject
powPEmw = z;    %power of test

pORmw = z;      %p-value for overall rating Mann-Whitney U-test
powORmw = z;    %power of test
hORmw = z;      %1 = reject null hypothesis, 0 = fail to reject

%Wilcoxon Signed Rank Test for v10 (task # = row index)
for i = 1:numTasks
    [pAC10sr(i),hAC10sr(i)] = signrank(Mac10(i,:),medAC);
    [pPE10sr(i),hPE10sr(i)] = signrank(Mpe10(i,:),medPE);
    [pOR10sr(i),hOR10sr(i)] = signrank(Mor10(i,:),medOR);
    powAC10sr(i) = signRankPower(median(Mac10(i,:),'omitnan'),...
        std(Mac10(i,:),'omitnan'),length(p10)*numTrials,medAC,N);
    powPE10sr(i) = signRankPower(median(Mpe10(i,:),'omitnan'),...
        std(Mpe10(i,:),'omitnan'),length(p10)*numTrials,medPE,N);
    powOR10sr(i) = signRankPower(median(Mor10(i,:),'omitnan'),...
        std(Mor10(i,:),'omitnan'),length(p10)*numTrials,medOR,N);
    
    %Combining all tasks into one p-value per scale
    [pAC10,hAC10] = signrank(reshape(Mac10,1,[]),medAC);
    powAC10 = signRankPower(median(Mac10,'all','omitnan'),...
        std(Mac10,0,'all','omitnan'),numel(Mac10),medAC,N);
    [pPE10,hPE10] = signrank(reshape(Mpe10,1,[]),medPE);
    powPE10 = signRankPower(median(Mpe10,'all','omitnan'),...
        std(Mpe10,0,'all','omitnan'),numel(Mpe10),medPE,N);
    [pOR10,hOR10] = signrank(reshape(Mor10,1,[]),medOR);
    powOR10 = signRankPower(median(Mor10,'all','omitnan'),...
        std(Mor10,0,'all','omitnan'),numel(Mor10),medOR,N);
end
%Wilcoxon Signed Rank Test for early versions (task # = row index)
for i = 1:numTasks
    [pAC0sr(i),hAC0sr(i)] = signrank(Mac0(i,:),medAC);
    [pPE0sr(i),hPE0sr(i)] = signrank(Mpe0(i,:),medPE);
    [pOR0sr(i),hOR0sr(i)] = signrank(Mor0(i,:),medOR);
    powAC0sr(i) = signRankPower(median(Mac0(i,:),'omitnan'),...
        std(Mac0(i,:),'omitnan'),length(p0)*numTrials,medAC,N);
    powPE0sr(i) = signRankPower(median(Mpe0(i,:),'omitnan'),...
        std(Mpe0(i,:),'omitnan'),length(p0)*numTrials,medPE,N);
    powOR0sr(i) = signRankPower(median(Mor0(i,:),'omitnan'),...
        std(Mor0(i,:),'omitnan'),length(p0)*numTrials,medOR,N);
    
    %Combining all tasks into one p-value per scale
    [pAC0,hAC0] = signrank(reshape(Mac0,1,[]),medAC);
    powAC0 = signRankPower(median(Mac0,'all','omitnan'),...
        std(Mac0,0,'all','omitnan'),numel(Mac0),medAC,N);
    [pPE0,hPE0] = signrank(reshape(Mpe0,1,[]),medPE);
    powPE0 = signRankPower(median(Mpe0,'all','omitnan'),...
        std(Mpe0,0,'all','omitnan'),numel(Mpe0),medPE,N);
    [pOR0,hOR0] = signrank(reshape(Mor0,1,[]),medOR);
    powOR0 = signRankPower(median(Mor0,'all','omitnan'),...
        std(Mor0,0,'all','omitnan'),numel(Mor0),medOR,N);
end

powAC10sr = round(powAC10sr,2); %round to nearest 0.01
powPE10sr = round(powPE10sr,2);
powOR10sr = round(powOR10sr,2);
powAC10 = round(powAC10,2);
powPE10 = round(powPE10,2);
powOR10 = round(powOR10,2);

powAC0sr = round(powAC0sr,2);   %round to nearest 0.01
powPE0sr = round(powPE0sr,2);
powOR0sr = round(powOR0sr,2);
powAC0 = round(powAC0,2);
powPE0 = round(powPE0,2);
powOR0 = round(powOR0,2);

%Mann-Whitney U-Test/Wilcoxon Rank Sum Test comparing v10 to before
for i = 1:numTasks
    [pACmw(i),hACmw(i)] = ranksum(Mac10(i,:),Mac0(i,:));
    [pPEmw(i),hPEmw(i)] = ranksum(Mpe10(i,:),Mpe0(i,:));
    [pORmw(i),hORmw(i)] = ranksum(Mor10(i,:),Mor0(i,:));
    powACmw(i) = rankSumPower(median(Mac10(i,:),'omitnan'),...
        std(Mac10(i,:),'omitnan'),length(p10)*numTrials,...
        median(Mac0(i,:),'omitnan'),std(Mac0(i,:),'omitnan'),...
        length(p0)*numTrials,N);
    powPEmw(i) = rankSumPower(median(Mpe10(i,:),'omitnan'),...
        std(Mpe10(i,:),'omitnan'),length(p10)*numTrials,...
        median(Mpe0(i,:),'omitnan'),std(Mpe0(i,:),'omitnan'),...
        length(p0)*numTrials,N);
    powORmw(i) = rankSumPower(median(Mor10(i,:),'omitnan'),...
        std(Mor10(i,:),'omitnan'),length(p10)*numTrials,...
        median(Mor0(i,:),'omitnan'),std(Mor0(i,:),'omitnan'),...
        length(p0)*numTrials,N);
end

powACmw = round(powACmw,2); %round to nearest 0.01
powPEmw = round(powPEmw,2);
powORmw = round(powORmw,2);

n10 = numel(Mac10); %total # of responses for v10
n0 = numel(Mac0);   %total # of responses for v09 and before

[pAC,hAC] = ranksum(reshape(Mac10,1,n10),reshape(Mac0,1,n0));
[pPE,hPE] = ranksum(reshape(Mpe10,1,n10),reshape(Mpe0,1,n0));
[pOR,hOR] = ranksum(reshape(Mor10,1,n10),reshape(Mor0,1,n0));
powAC = rankSumPower(median(Mac10,'all','omitnan'),...
    std(Mac10,0,'all','omitnan'),numel(Mac10),...
    median(Mac0,'all','omitnan'),std(Mac0,0,'all','omitnan'),...
    numel(Mac0),N);
powPE = rankSumPower(median(Mpe10,'all','omitnan'),...
    std(Mpe10,0,'all','omitnan'),numel(Mpe10),...
    median(Mpe0,'all','omitnan'),std(Mpe0,0,'all','omitnan'),...
    numel(Mpe0),N);
powOR = rankSumPower(median(Mac10,'all','omitnan'),...
    std(Mor10,0,'all','omitnan'),numel(Mor10),...
    median(Mor0,'all','omitnan'),std(Mor0,0,'all','omitnan'),...
    numel(Mor0),N);

powAC = round(powAC,2); %round to nearest 0.01
powPE = round(powPE,2);
powOR = round(powOR,2);

%Plotting Wilcoxon Signed Rank Results v10
figure
hold on
b1 = bar([pAC10sr,pPE10sr,pOR10sr]);
xlim([0.5,numTasks+0.5])
grid on
[x1,x2,x3] = b1.XEndPoints; %get x-location of tops of bars
[y1,y2,y3] = b1.YEndPoints; %get y-location of tops of bars
xtips = [x1,x2,x3];
ytips = [y1,y2,y3];
labels = pad(string([powAC10sr,powPE10sr,powOR10sr]),8,'left');
text(xtips,ytips,labels,...     %add power values to plot
    'Rotation',90,...           %rotate text 90 degrees CCW
    'VerticalAlignment','middle')
if max(ytips)>alpha    %only plot significance threshold if one bar is over
    plot([0.5,numTasks+0.5],ones(1,2)*alpha,'k--') %plot alpha
end
xlabel('Task #')
ylabel('p-Value')
title('Wilcoxon Signed Rank Test PA-28R-201 v10')
legend('Aircraft Characteristics','Pilot Effort','Overall Rating',...
    '\alpha = 0.05')

%Plotting Wilcoxon Signed Rank Results v09 and earlier
figure
hold on
b2 = bar([pAC0sr,pPE0sr,pOR0sr]);
xlim([0.5,numTasks+0.5])
grid on
[x1,x2,x3] = b2.XEndPoints; %get x-location of tops of bars
[y1,y2,y3] = b2.YEndPoints; %get y-location of tops of bars
xtips = [x1,x2,x3];
ytips = [y1,y2,y3];
labels = pad(string([powAC0sr,powPE0sr,powOR0sr]),8,'left');
text(xtips,ytips,labels,...     %add power values to plot
    'Rotation',90,...           %rotate text 90 degrees CCW
    'VerticalAlignment','middle')
if max(ytips)>alpha    %only plot significance threshold if one bar is over
    plot([0.5,numTasks+0.5],ones(1,2)*alpha,'k--') %plot alpha
end
xlabel('Task #')
ylabel('p-Value')
title('Wilcoxon Signed Rank Test PA-28R-201 Early Models')
legend('Aircraft Characteristics','Pilot Effort','Overall Rating',...
    '\alpha = 0.05')

%Plotting Mann-Whitney U-Test results comparing v10 and others by task
figure
hold on
b3 = bar([pACmw,pPEmw,pORmw]);
xlim([0.5,numTasks+0.5])
grid on
[x1,x2,x3] = b3.XEndPoints; %get x-locations of tops of bars
[y1,y2,y3] = b3.YEndPoints; %get y-locations of tops of bars
xtips = [x1,x2,x3];
ytips = [y1,y2,y3];
labels = pad(string([powACmw,powPEmw,powORmw]),8,'left');
text(xtips,ytips,labels,...     %add power values to plot
    'Rotation',90,...           %rotate text 90 degrees CCW
    'VerticalAlignment','middle')
if max(ytips)>alpha    %only plot significance threshold if one bar is over
    plot([0.5,numTasks+0.5],ones(1,2)*alpha,'k--') %plot alpha
end
xlabel('Task #')
ylabel('p-Value')
title('Mann-Whitney U-Test/Wilcoxon Rank Sum Test')
legend('Aircraft Characteristics','Pilot Effort','Overall Rating',...
    '\alpha = 0.05')

%Plotting combined Wilcoxon Signed Rank results v10
figure
hold on
bars = categorical({'Aircraft Characteristics','Pilot Effort',...
       'Overall Rating'});
b4 = bar(bars,[pAC10,pPE10,pOR10]);
b4.FaceColor = 'flat';  %allow for setting color of bars
b4.CData = eye(3);  %make bars red, green, blue
grid on
ylabel('p-Value')
title('Wilcoxon Signed Rank Test PA-28R-201 v10')
xtips = b4.XEndPoints;  %get x-locations of tops of bars
ytips = b4.YEndPoints;  %get y-locations of tops of bars
labels = string([powAC10,powPE10,powOR10]); %convert Power to string array
text(xtips,ytips,labels,...             %add power values to plot
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

%Plotting combined Wilcoxon Signed Rank results, early versions
figure
hold on
bars = categorical({'Aircraft Characteristics','Pilot Effort',...
       'Overall Rating'});
b5 = bar(bars,[pAC0,pPE0,pOR0]);
b5.FaceColor = 'flat';  %allow for setting color of bars
b5.CData = eye(3);  %make bars red, green, blue
grid on
ylabel('p-Value')
title('Wilcoxon Signed Rank Test PA-28R-201 Early Models')
xtips = b5.XEndPoints;  %get x-locations of tops of bars
ytips = b5.YEndPoints;  %get y-locations of tops of bars
labels = string([powAC0,powPE0,powOR0]);    %convert Power to string array
text(xtips,ytips,labels,...             %add power values to plot
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

%Plotting Mann-Whitney U-Test results for all v10 and previous data
figure
hold on
bars = categorical({'Aircraft Characteristics','Pilot Effort',...
       'Overall Rating'});
b6 = bar(bars,[pAC,pPE,pOR]);
b6.FaceColor = 'flat';  %allow for setting color of bars
b6.CData = eye(3);  %make bars red, green, blue
grid on
ylabel('p-Value')
title('Mann-Whitney U-Test/Wilcoxon Rank Sum Test')
xtips = b6.XEndPoints;  %get x-locations of tops of bars
ytips = b6.YEndPoints;  %get y-locations of tops of bars
labels = string([powAC,powPE,powOR]);   %convert Power to string array
text(xtips,ytips,labels,...             %add power values to plot
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%---------------------+-----------------------------+---------------------%
%                     | Stacked Divergent Bar Plots |                     %
%---------------------+-----------------------------+---------------------%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure Sizing
fw = 6;     %figure width in inches
fh = 2.25;  %figure height in inches

%Aircraft Characteristics Ratings
fAC = figure;   %create figure and get handle
subplot(1,2,1)
axAC0 = gca; %get current axes object
[bAC0,statsAC0] = stackDivPlot(axAC0,Mac0,maxAC,medAC);%early model ratings
title('Early Models')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mac0(i,:))));
    text(axAC0.XLim(2),i,txt)
end
legend('off')   %use legend generated from second subplot
axAC0.OuterPosition = [0 0 0.35 1];   %shrink plotbox to fit legend & text
axAC0.XGrid = 'on';  %turn vertical gridlines on

subplot(1,2,2)
axAC = gca; %get current axes object
[bAC,statsAC] = stackDivPlot(axAC,Mac10,maxAC,medAC);   %v10 ratings
title('Piper\_PA-28R-201\_v10')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mac10(i,:))));
    text(axAC.XLim(2),i,txt)
end
l = legend;     %get legend object
oldpos = l.Position;    %get current legend position
l.Position = [1-oldpos(3)-0.01,oldpos(2:4)];    %shift legend right
axAC.OuterPosition = [0.45 0 0.35 1]; %shrink plotbox to fit legend & text
axAC.XGrid = 'on';  %turn vertical gridlines on

fAC.Units = 'inches';   %set units to inches
oldpos = fAC.Position;  %get current figure window position
fAC.Position = [oldpos(1:2),fw,fh]; %resize window, keep same location

%Pilot Effort Ratings
fPE = figure;   %create figure and get handle
subplot(1,2,1)
axPE0 = gca; %get current axes object
[bPE0,statsPE0] = stackDivPlot(axPE0,Mpe0,maxPE,medPE);
title('Early Models')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mpe0(i,:))));
    text(axPE0.XLim(2),i,txt)
end
legend('off')   %use legend generated from second subplot
axPE0.OuterPosition = [0 0 0.35 1];   %shrink plotbox to fit legend & text
axPE0.XGrid = 'on';  %turn vertical gridlines on

subplot(1,2,2)
axPE = gca; %get current axes object
[bPE,statsPE] = stackDivPlot(axPE,Mpe10,maxPE,medPE);
title('Piper\_PA-28R-201\_v10')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mpe10(i,:))));
    text(axPE.XLim(2),i,txt)
end
l = legend;     %get legend object
oldpos = l.Position;    %get current legend position
l.Position = [1-oldpos(3)-0.01,oldpos(2:4)];    %shift legend right
axPE.OuterPosition = [0.45 0 0.35 1]; %shrink plotbox to fit legend & text
axPE.XGrid = 'on';  %turn vertical gridlines on

fPE.Units = 'inches';   %set units to inches
oldpos = fPE.Position;  %get current figure window position
fPE.Position = [oldpos(1:2),fw,fh]; %resize window, keep same location

%Overall Ratings
fOR = figure;   %create figure and get handle
subplot(1,2,1)
axOR0 = gca; %get current axes object
[bOR0,statsOR0] = stackDivPlot(axOR0,Mor0,maxOR,medOR);
title('Early Models')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mor0(i,:))));
    text(axOR0.XLim(2),i,txt)
end
legend('off')   %use legend generated from second subplot
axOR0.OuterPosition = [0 0 0.35 1];   %shrink plotbox to fit legend & text
axOR0.XGrid = 'on';  %turn vertical gridlines on

subplot(1,2,2)
axOR = gca; %get current axes object
[bOR,statsOR] = stackDivPlot(axOR,Mor10,maxOR,medOR);
title('Piper\_PA-28R-201\_v10')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mor10(i,:))));
    text(axOR.XLim(2),i,txt)
end
l = legend;     %get legend object
oldpos = l.Position;    %get current legend position
l.Position = [1-oldpos(3)-0.01,oldpos(2:4)];    %shift legend right
axOR.OuterPosition = [0.45 0 0.35 1]; %shrink plotbox to fit legend & text
axOR.XGrid = 'on';  %turn vertical gridlines on

fOR.Units = 'inches';   %set units to inches
oldpos = fOR.Position;  %get current figure window position
fOR.Position = [oldpos(1:2),fw,fh]; %resize window, keep same location


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pow = signRankPower(sampmean,sampstd,sampsize,hypmean,N)
%Calculates the power of a Wilcoxon Signed Rank Test using a Monte-Carlo
%Analysis of normally distributed samples with the same mean and variance
%as the sample of survey responses.
%Inputs:
%   sampmean = mean of the sample
%   sampstd = standard deviation of the sample
%   sampsize = number of data points in the sample
%   hypmean = hypothesized mean (default = 0)
%   N = number of Monte-Carlo samples to generate (default = 1000)
%Outputs:
%   pow = power of the test = 1-P(Type II Error)

if nargin<4
    hypmean = 0;    %default hypothesized mean
end
if nargin<5
    N = 1000;   %default number of Monte-Carlo cases
end

h = zeros(N,1); %preallocate logical array for rejecting null hypothesis
for i = 1:N
    dat = normrnd(sampmean,sampstd,[sampsize,1]);   %normally distributed
    [~,h(i)] = signrank(dat,hypmean);   %apply Wilcoxon Signed Rank Test
end

pow = sum(h)/N; %probability that we reject null hyp (we know it's false
                    %since data was generated with a mean =/= hypmean)
end

function pow = rankSumPower(m1,s1,n1,m2,s2,n2,N)
%Computes the power of a Wilcoxon Rank Sum Test using a Monte-Carlo
%Analysis of normally distributed samples with the same mean and variance
%as the samples of survey responses.
%Inputs:
%   m1 = sample mean of group 1
%   s1 = sample standard deviation of group 1
%   n1 = size of group 1
%   m2 = sample mean of group 2
%   s2 = sample standard deviation of group 2
%   n2 = size of group 2
%   N = number of Monte-Carlo samples to generate (default = 1000)
%Outputs:
%   pow = power of the test = 1-P(Type II Error)

if nargin<7
    N = 1000;   %default number of Monte-Carlo cases
end

h = zeros(N,1); %preallocate logical array for rejecting null hypothesis
for i = 1:N
    dat1 = normrnd(m1,s1,[n1,1]);   %normally distributed
    dat2 = normrnd(m2,s2,[n2,1]);   %normally distributed
    [~,h(i)] = ranksum(dat1,dat2);  %apply Wilcoxon Rank Sum Test
end

pow = sum(h)/N; %probability that we reject null hyp (we know it's false as
                    %long as m1 =/= m2)
end



