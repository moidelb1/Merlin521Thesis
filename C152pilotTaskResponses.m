%Benjamin Moidel
%December 23, 2020
%Pilot Task Evaluation Responses

%% Loading data from Excel
clc;clear;close all

fname = 'Task Evaluation Survey C152 (Responses) 20201223.xlsx';
fpath = 'D:\Storage\Documents\School Stuff\Grad School\Research';
f = fullfile(fpath,fname);

%aircraft characteristics responses
Mac = readmatrix(f,'Sheet','Sorted Form Responses','Range','C31:K54');
%pilot effort responses
Mpe = readmatrix(f,'Sheet','Sorted Form Responses','Range','C58:K81');
%overall rating responses
Mor = readmatrix(f,'Sheet','Sorted Form Responses','Range','C85:K108');

save C152Responses.mat Mac Mpe Mor

%% Analyzing and Plotting data
clc;clear;close all

load C152Responses.mat

%Data dimensions
numPilots = 9;      %total # of subjects that flew models
numTasks = 8;       %total # of tasks/conditions each subject flew
numTrials = 3;      %total # of trials per task
p18 = 1:3;          %subject #'s that flew C152 v18
p19 = 4:numPilots;  %subjects that flew C152 v19
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
Mac18 = reshape(Mac(:,p18)',length(p18)*numTrials,numTasks)';
Mpe18 = reshape(Mpe(:,p18)',length(p18)*numTrials,numTasks)';
Mor18 = reshape(Mor(:,p18)',length(p18)*numTrials,numTasks)';
Mac19 = reshape(Mac(:,p19)',length(p19)*numTrials,numTasks)';
Mpe19 = reshape(Mpe(:,p19)',length(p19)*numTrials,numTasks)';
Mor19 = reshape(Mor(:,p19)',length(p19)*numTrials,numTasks)';

%Preallocation
z = zeros(numTasks,1);  %dummy zero vector with the correct size

pAC19sr = z;    %p-value for v19 aircraft characteristics signed rank
hAC19sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powAC19sr = z;  %power of test

pPE19sr = z;    %p-value for v19 pilot effort signed rank
hPE19sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powPE19sr = z;  %power of test

pOR19sr = z;    %p-value for v19 overall rating signed rank
hOR19sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powOR19sr = z;  %power of test

pAC18sr = z;    %p-value for v18 aircraft characteristics signed rank
hAC18sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powAC18sr = z;  %power of test

pPE18sr = z;    %p-value for v18 pilot effort signed rank
hPE18sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powPE18sr = z;  %power of test

pOR18sr = z;    %p-value for v18 overall rating signed rank
hOR18sr = z;    %1 = reject null hypothesis, 0 = fail to reject
powOR18sr = z;  %power of test

pACmw = z;      %p-value for aircraft characteristics Mann-Whitney U-test
hACmw = z;      %1 = reject null hypothesis, 0 = fail to reject
powACmw = z;    %power of test

pPEmw = z;      %p-value for pilot effort Mann-Whitney U-test
hPEmw = z;      %1 = reject null hypothesis, 0 = fail to reject
powPEmw = z;    %power of test

pORmw = z;      %p-value for overall rating Mann-Whitney U-test
powORmw = z;    %power of test
hORmw = z;      %1 = reject null hypothesis, 0 = fail to reject

%Wilcoxon Signed Rank Test for v19 (task # = row index)
for i = 1:numTasks
    [pAC19sr(i),hAC19sr(i)] = signrank(Mac19(i,:),medAC);
    [pPE19sr(i),hPE19sr(i)] = signrank(Mpe19(i,:),medPE);
    [pOR19sr(i),hOR19sr(i)] = signrank(Mor19(i,:),medOR);
    powAC19sr(i) = signRankPower(median(Mac19(i,:),'omitnan'),...
        std(Mac19(i,:),'omitnan'),length(p19)*numTrials,medAC,N);
    powPE19sr(i) = signRankPower(median(Mpe19(i,:),'omitnan'),...
        std(Mpe19(i,:),'omitnan'),length(p19)*numTrials,medPE,N);
    powOR19sr(i) = signRankPower(median(Mor19(i,:),'omitnan'),...
        std(Mor19(i,:),'omitnan'),length(p19)*numTrials,medOR,N);
    
    %Combining all tasks into one p-value per scale
    [pAC19,hAC19] = signrank(reshape(Mac19,1,[]),medAC);
    powAC19 = signRankPower(median(Mac19,'all','omitnan'),...
        std(Mac19,0,'all','omitnan'),numel(Mac19),medAC,N);
    [pPE19,hPE19] = signrank(reshape(Mpe19,1,[]),medPE);
    powPE19 = signRankPower(median(Mpe19,'all','omitnan'),...
        std(Mpe19,0,'all','omitnan'),numel(Mpe19),medPE,N);
    [pOR19,hOR19] = signrank(reshape(Mor19,1,[]),medOR);
    powOR19 = signRankPower(median(Mor19,'all','omitnan'),...
        std(Mor19,0,'all','omitnan'),numel(Mor19),medOR,N);
end
%Wilcoxon Signed Rank Test for v18 (task # = row index)
for i = 1:numTasks
    [pAC18sr(i),hAC18sr(i)] = signrank(Mac18(i,:),medAC);
    [pPE18sr(i),hPE18sr(i)] = signrank(Mpe18(i,:),medPE);
    [pOR18sr(i),hOR18sr(i)] = signrank(Mor18(i,:),medOR);
    powAC18sr(i) = signRankPower(median(Mac18(i,:),'omitnan'),...
        std(Mac18(i,:),'omitnan'),length(p18)*numTrials,medAC,N);
    powPE18sr(i) = signRankPower(median(Mpe18(i,:),'omitnan'),...
        std(Mpe18(i,:),'omitnan'),length(p18)*numTrials,medPE,N);
    powOR18sr(i) = signRankPower(median(Mor18(i,:),'omitnan'),...
        std(Mor18(i,:),'omitnan'),length(p18)*numTrials,medOR,N);
    
    %Combining all tasks into one p-value per scale
    [pAC18,hAC18] = signrank(reshape(Mac18,1,[]),medAC);
    powAC18 = signRankPower(median(Mac18,'all','omitnan'),...
        std(Mac18,0,'all','omitnan'),numel(Mac18),medAC,N);
    [pPE18,hPE18] = signrank(reshape(Mpe18,1,[]),medPE);
    powPE18 = signRankPower(median(Mpe18,'all','omitnan'),...
        std(Mpe18,0,'all','omitnan'),numel(Mpe18),medPE,N);
    [pOR18,hOR18] = signrank(reshape(Mor18,1,[]),medOR);
    powOR18 = signRankPower(median(Mor18,'all','omitnan'),...
        std(Mor18,0,'all','omitnan'),numel(Mor18),medOR,N);
end

powAC19sr = round(powAC19sr,2); %round to nearest 0.01
powPE19sr = round(powPE19sr,2);
powOR19sr = round(powOR19sr,2);
powAC19 = round(powAC19,2);
powPE19 = round(powPE19,2);
powOR19 = round(powOR19,2);

powAC18sr = round(powAC18sr,2); %round to nearest 0.01
powPE18sr = round(powPE18sr,2);
powOR18sr = round(powOR18sr,2);
powAC18 = round(powAC18,2);
powPE18 = round(powPE18,2);
powOR18 = round(powOR18,2);

%Mann-Whitney U-Test/Wilcoxon Rank Sum Test comparing v18 and v19
for i = 1:numTasks
    [pACmw(i),hACmw(i)] = ranksum(Mac19(i,:),Mac18(i,:));
    [pPEmw(i),hPEmw(i)] = ranksum(Mpe19(i,:),Mpe18(i,:));
    [pORmw(i),hORmw(i)] = ranksum(Mor19(i,:),Mor18(i,:));
    powACmw(i) = rankSumPower(median(Mac19(i,:),'omitnan'),...
        std(Mac19(i,:),'omitnan'),length(p19)*numTrials,...
        median(Mac18(i,:),'omitnan'),std(Mac18(i,:),'omitnan'),...
        length(p18)*numTrials,N);
    powPEmw(i) = rankSumPower(median(Mpe19(i,:),'omitnan'),...
        std(Mpe19(i,:),'omitnan'),length(p19)*numTrials,...
        median(Mpe18(i,:),'omitnan'),std(Mpe18(i,:),'omitnan'),...
        length(p18)*numTrials,N);
    powORmw(i) = rankSumPower(median(Mor19(i,:),'omitnan'),...
        std(Mor19(i,:),'omitnan'),length(p19)*numTrials,...
        median(Mor18(i,:),'omitnan'),std(Mor18(i,:),'omitnan'),...
        length(p18)*numTrials,N);
end

powACmw = round(powACmw,2); %round to nearest 0.01
powPEmw = round(powPEmw,2);
powORmw = round(powORmw,2);

n19 = numel(Mac19); %total # of responses for v19
n18 = numel(Mac18); %total # of responses for v18

[pAC,hAC] = ranksum(reshape(Mac19,1,n19),reshape(Mac18,1,n18));
[pPE,hPE] = ranksum(reshape(Mpe19,1,n19),reshape(Mpe18,1,n18));
[pOR,hOR] = ranksum(reshape(Mor19,1,n19),reshape(Mor18,1,n18));
powAC = rankSumPower(median(Mac19,'all','omitnan'),...
    std(Mac19,0,'all','omitnan'),numel(Mac19),...
    median(Mac18,'all','omitnan'),std(Mac18,0,'all','omitnan'),...
    numel(Mac18),N);
powPE = rankSumPower(median(Mpe19,'all','omitnan'),...
    std(Mpe19,0,'all','omitnan'),numel(Mpe19),...
    median(Mpe18,'all','omitnan'),std(Mpe18,0,'all','omitnan'),...
    numel(Mpe18),N);
powOR = rankSumPower(median(Mac19,'all','omitnan'),...
    std(Mor19,0,'all','omitnan'),numel(Mor19),...
    median(Mor18,'all','omitnan'),std(Mor18,0,'all','omitnan'),...
    numel(Mor18),N);

powAC = round(powAC,2); %round to nearest 0.01
powPE = round(powPE,2);
powOR = round(powOR,2);

%Plotting Wilcoxon Signed Rank Results v19
figure
hold on
b1 = bar([pAC19sr,pPE19sr,pOR19sr]);
% plot([0.5,numTasks+0.5],ones(1,2)*alpha,'k--') %plot alpha
xlim([0.5,numTasks+0.5])
grid on
xlabel('Task #')
ylabel('p-Value')
title('Wilcoxon Signed Rank Test C152 v19')
legend('Aircraft Characteristics','Pilot Effort','Overall Rating')
[x1,x2,x3] = b1.XEndPoints; %get x-location of tops of bars
[y1,y2,y3] = b1.YEndPoints; %get y-location of tops of bars
xtips = [x1,x2,x3];
ytips = [y1,y2,y3];
labels = pad(string([powAC19sr,powPE19sr,powOR19sr]),8,'left');
text(xtips,ytips,labels,...     %add power values to plot
    'Rotation',90,...           %rotate text 90 degrees CCW
    'VerticalAlignment','middle')

%Plotting Wilcoxon Signed Rank Results v18
figure
hold on
b2 = bar([pAC18sr,pPE18sr,pOR18sr]);
plot([0.5,numTasks+0.5],ones(1,2)*alpha,'k--') %plot alpha
xlim([0.5,numTasks+0.5])
grid on
xlabel('Task #')
ylabel('p-Value')
title('Wilcoxon Signed Rank Test C152 v18')
legend('Aircraft Characteristics','Pilot Effort','Overall Rating')
[x1,x2,x3] = b2.XEndPoints; %get x-location of tops of bars
[y1,y2,y3] = b2.YEndPoints; %get y-location of tops of bars
xtips = [x1,x2,x3];
ytips = [y1,y2,y3];
labels = pad(string([powAC18sr,powPE18sr,powOR18sr]),8,'left');
text(xtips,ytips,labels,...     %add power values to plot
    'Rotation',90,...           %rotate text 90 degrees CCW
    'VerticalAlignment','middle')

%Plotting Mann-Whitney U-Test results comparing v18 & v19 by task
figure
hold on
b3 = bar([pACmw,pPEmw,pORmw]);
plot([0.5,numTasks+0.5],ones(1,2)*alpha,'k--') %plot alpha
xlim([0.5,numTasks+0.5])
grid on
xlabel('Task #')
ylabel('p-Value')
title('Mann-Whitney U-Test/Wilcoxon Rank Sum Test v19 vs v18')
legend('Aircraft Characteristics','Pilot Effort','Overall Rating',...
    '\alpha = 0.05')
[x1,x2,x3] = b3.XEndPoints; %get x-locations of tops of bars
[y1,y2,y3] = b3.YEndPoints; %get y-locations of tops of bars
xtips = [x1,x2,x3];
ytips = [y1,y2,y3];
labels = pad(string([powACmw,powPEmw,powORmw]),8,'left');
text(xtips,ytips,labels,...     %add power values to plot
    'Rotation',90,...           %rotate text 90 degrees CCW
    'VerticalAlignment','middle')

%Plotting combined Wilcoxon Signed Rank results v19
figure
hold on
bars = categorical({'Aircraft Characteristics','Pilot Effort',...
       'Overall Rating'});
b4 = bar(bars,[pAC19,pPE19,pOR19]);
b4.FaceColor = 'flat';  %allow for setting color of bars
b4.CData = eye(3);  %make bars red, green, blue
grid on
ylabel('p-Value')
title('Wilcoxon Signed Rank Test C152 v19')
xtips = b4.XEndPoints; %get x-location of tops of bars
ytips = b4.YEndPoints; %get y-location of tops of bars
labels = string([powAC19,powPE19,powOR19]); %convert Power to string array
text(xtips,ytips,labels,...     %add power values to plot
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

%Plotting combined Wilcoxon Signed Rank results v18
figure
hold on
bars = categorical({'Aircraft Characteristics','Pilot Effort',...
       'Overall Rating'});
b5 = bar(bars,[pAC18,pPE18,pOR18]);
b5.FaceColor = 'flat';  %allow for setting color of bars
b5.CData = eye(3);  %make bars red, green, blue
grid on
ylabel('p-Value')
title('Wilcoxon Signed Rank Test C152 v18')
xtips = b5.XEndPoints; %get x-location of tops of bars
ytips = b5.YEndPoints; %get y-location of tops of bars
labels = string([powAC18,powPE18,powOR18]); %convert Power to string array
text(xtips,ytips,labels,...     %add power values to plot
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

%Plotting Mann-Whitney U-Test results for all v18 & v19 data
figure
hold on
bars = categorical({'Aircraft Characteristics','Pilot Effort',...
       'Overall Rating'});
b6 = bar(bars,[pAC,pPE,pOR]);
b6.FaceColor = 'flat';  %allow for setting color of bars
b6.CData = eye(3);  %make bars red, green, blue
grid on
ylabel('p-Value')
title('Mann-Whitney U-Test/Wilcoxon Rank Sum Test v19 vs v18')
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
axAC18 = gca; %get current axes object
[bAC18,statsAC18] = stackDivPlot(axAC18,Mac18,maxAC,medAC); %v18 ratings
title('Cessna\_152\_v18')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mac18(i,:))));
    text(axAC18.XLim(2),i,txt)
end
legend('off')   %use legend generated from second subplot
axAC18.OuterPosition = [0 0 0.35 1];   %shrink plotbox to fit legend & text
axAC18.XGrid = 'on';  %turn vertical gridlines on

subplot(1,2,2)
axAC = gca; %get current axes object
[bAC,statsAC] = stackDivPlot(axAC,Mac19,maxAC,medAC);   %v19 ratings
title('Cessna\_152\_v19')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mac19(i,:))));
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
axPE18 = gca; %get current axes object
[bPE18,statsPE18] = stackDivPlot(axPE18,Mpe18,maxPE,medPE);
title('Cessna\_152\_v18')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mpe18(i,:))));
    text(axPE18.XLim(2),i,txt)
end
legend('off')   %use legend generated from second subplot
axPE18.OuterPosition = [0 0 0.35 1];   %shrink plotbox to fit legend & text
axPE18.XGrid = 'on';  %turn vertical gridlines on

subplot(1,2,2)
axPE = gca; %get current axes object
[bPE,statsPE] = stackDivPlot(axPE,Mpe19,maxPE,medPE);
title('Cessna\_152\_v19')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mpe19(i,:))));
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
axOR18 = gca; %get current axes object
[bOR18,statsOR18] = stackDivPlot(axOR18,Mor18,maxOR,medOR);
title('Cessna\_152\_v18')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mor18(i,:))));
    text(axOR18.XLim(2),i,txt)
end
legend('off')   %use legend generated from second subplot
axOR18.OuterPosition = [0 0 0.35 1];   %shrink plotbox to fit legend & text
axOR18.XGrid = 'on';  %turn vertical gridlines on

subplot(1,2,2)
axOR = gca; %get current axes object
[bOR,statsOR] = stackDivPlot(axOR,Mor19,maxOR,medOR);
title('Cessna\_152\_v19')
xlabel('# of Responses')
ylabel('Condition #')
for i = 1:numTasks
    txt = sprintf('  N = %i',sum(~isnan(Mor19(i,:))));
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



% %Create Bar Plots
% [bAC,statsAC] = stackDivPlot(Mac19,maxAC);
% fAC = gcf;  %get current figure object
% axAC = gca; %get current axes object
% title('Aircraft Characteristics Ratings v19')
% xlabel('# of Responses')
% ylabel('Condition #')
% for i = 1:numTasks
%     txt = sprintf('  N = %i',sum(~isnan(Mac19(i,:))));
%     text(axAC.XLim(2),i,txt)
% end
% l = legend;     %get legend object
% oldpos = l.Position;    %get current legend position
% l.Position = [1-oldpos(3)-0.01,oldpos(2:4)];    %shift legend right
% axAC.OuterPosition = [0 0 0.8 1];   %shrink plotbox to fit legend & text
% axAC.XGrid = 'on';  %turn vertical gridlines on
% 
% [bPE,statsPE] = stackDivPlot(Mpe19,maxPE);
% fPE = gcf;  %get current figure object
% axPE = gca; %get current axes object
% title('Pilot Effort Ratings v19')
% xlabel('# of Responses')
% ylabel('Condition #')
% for i = 1:numTasks
%     txt = sprintf('  N = %i',sum(~isnan(Mpe19(i,:))));
%     text(axPE.XLim(2),i,txt)
% end
% l = legend;     %get legend object
% oldpos = l.Position;    %get current legend position
% l.Position = [1-oldpos(3)-0.01,oldpos(2:4)];    %shift legend right
% axPE.OuterPosition = [0 0 0.8 1];   %shrink plotbox to fit legend & text
% axPE.XGrid = 'on';  %turn vertical gridlines on
% 
% [bOR,statsOR] = stackDivPlot(Mor19,maxOR);
% fOR = gcf;  %get current figure object
% axOR = gca; %get current axes object
% title('Overall Ratings v19')
% xlabel('# of Responses')
% ylabel('Condition #')
% for i = 1:numTasks
%     txt = sprintf('  N = %i',sum(~isnan(Mor19(i,:))));
%     text(axOR.XLim(2),i,txt)
% end
% l = legend;     %get legend object
% oldpos = l.Position;    %get current legend position
% l.Position = [1-oldpos(3)-0.01,oldpos(2:4)];    %shift legend right
% axOR.OuterPosition = [0 0 0.8 1];   %shrink plotbox to fit legend & text
% axOR.XGrid = 'on';  %turn vertical gridlines on
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% v18 Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Create Bar Plots
% [bAC18,statsAC18] = stackDivPlot(Mac18,maxAC);
% fAC18 = gcf;  %get current figure object
% axAC18 = gca; %get current axes object
% title('Aircraft Characteristics Ratings v18')
% xlabel('# of Responses')
% ylabel('Condition #')
% for i = 1:numTasks
%     txt = sprintf('  N = %i',sum(~isnan(Mac18(i,:))));
%     text(axAC18.XLim(2),i,txt)
% end
% l = legend;     %get legend object
% oldpos = l.Position;    %get current legend position
% l.Position = [1-oldpos(3)-0.01,oldpos(2:4)];    %shift legend right
% axAC18.OuterPosition = [0 0 0.8 1];   %shrink plotbox to fit legend & text
% axAC18.XGrid = 'on';  %turn vertical gridlines on
% 
% [bPE18,statsPE18] = stackDivPlot(Mpe18,maxPE);
% fPE18 = gcf;  %get current figure object
% axPE18 = gca; %get current axes object
% title('Pilot Effort Ratings v18')
% xlabel('# of Responses')
% ylabel('Condition #')
% for i = 1:numTasks
%     txt = sprintf('  N = %i',sum(~isnan(Mpe18(i,:))));
%     text(axPE18.XLim(2),i,txt)
% end
% l = legend;     %get legend object
% oldpos = l.Position;    %get current legend position
% l.Position = [1-oldpos(3)-0.01,oldpos(2:4)];    %shift legend right
% axPE18.OuterPosition = [0 0 0.8 1];   %shrink plotbox to fit legend & text
% axPE18.XGrid = 'on';  %turn vertical gridlines on
% 
% [bOR18,statsOR18] = stackDivPlot(Mor18,maxOR);
% fOR18 = gcf;  %get current figure object
% axOR18 = gca; %get current axes object
% title('Overall Ratings v18')
% xlabel('# of Responses')
% ylabel('Condition #')
% for i = 1:numTasks
%     txt = sprintf('  N = %i',sum(~isnan(Mor18(i,:))));
%     text(axOR18.XLim(2),i,txt)
% end
% l = legend;     %get legend object
% oldpos = l.Position;    %get current legend position
% l.Position = [1-oldpos(3)-0.01,oldpos(2:4)];    %shift legend right
% axOR18.OuterPosition = [0 0 0.8 1];   %shrink plotbox to fit legend & text
% axOR18.XGrid = 'on';  %turn vertical gridlines on

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



