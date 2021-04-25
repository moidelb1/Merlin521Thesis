%Benjamin Moidel
%January 4, 2021
%Non-parametric statistical test visualization

close all

%% Wilcoxon Rank Sum Test/Mann-Whitney U-Test
clc;clear

%Inputs
n1 = 18;    %size of sample 1
n2 = 9;     %size of sample 2
alpha = 0.05;   %level of significance

%Calculating U-statistic
n = n1+n2;  %total number of samples
ns = min([n1 n2]);  %smaller sample size
nb = n-ns;          %bigger sample size
ranks = 1:n;  %ranks that can be assigned

poss = nchoosek(ranks,ns);  %all possible rank assignments
Rs = sum(poss,2);       %all possible rank sums for smaller sample
% Rb = n*(n+1)/2 - Rs;    %all possible rank sums for bigger sample

Us = ns*nb + ns*(ns+1)/2 - Rs;  %all possible U-statistics, smaller sample
Ub = n1*n2 - Us;                %all possible U-statistics, bigger sample

U = min([Us,Ub],[],2);  %calculate U-statistic for each case

bins = 0:(n1*n2/2);   %all possible integer U-values, # of "bins"
binlims = (0:(n1*n2/2+1))-0.5;  %min & max interval limits for each "bin"
h = histcounts(U,binlims);  %count # of occurances of each U-value
hprob = h/numel(U);     %probability for each "bin"
hcdf = cumsum(hprob);   %cumulative probability for <= each "bin"

Ucrit = floor(interp1(hcdf,bins,alpha));    %maximum U-value to reject null

%Plot PDF
figure
hold on
grid on
b1 = bar(bins,hprob);       %probability density function
title(sprintf('Probability Density, n_1=%i, n_2=%i',n1,n2))
xlabel('U-Statistic')
ylabel('Probability')
b0 = bar(NaN(2,1));         %dummy bar plot for legend entries
b0.FaceColor = [1 0 0];     %red
legend(b0,'Critical U-Statistic','Location','northwest')

%Plot CDF
figure
hold on
grid on
b2 = bar(bins,hcdf);        %cumulative distribution function
title(sprintf('Cumulative Distribution, n_1=%i, n_2=%i',n1,n2))
alphaline = plot([0,n1*n2/2],ones(1,2)*alpha,'k--');
xlabel('U-Statistic')
ylabel('Cumulative Probability')
b0 = bar(NaN(2,1));         %dummy bar plot for legend entries
b0.FaceColor = [1 0 0];     %red
legend([b0,alphaline],{'Critical U-Stat','\alpha = 0.05'},...
    'Location','northwest')

%Highlighting critical U-statistic
if ~isnan(Ucrit)
    b1.FaceColor = 'flat';
    b1.CData(Ucrit+1,:) = [1 0 0];  %make critical U-statistic red
    b2.FaceColor = 'flat';
    b2.CData(Ucrit+1,:) = [1 0 0];  %make critical U-statistic red
end

%% Wilcoxon Signed Rank Test
clc;clear

%Inputs
n = 18;     %size of sample
alpha = 0.05;   %level of significance

%Calculating W-statistic
ranks = 1:n;    %ranks that can be assigned

bincombs = dec2bin(0:2^n-1)-'0';    %permutations of 0 and 1 of length n
                                    %1 represents positive ranks, 0 = neg
poss = ranks.*bincombs;     %all possible positive rank assignments
                    %all possible negative rank assignments are not needed

Wpos = sum(poss,2);         %determine W+ (positive rank sum)
Wneg = n*(n+1)/2-Wpos;      %determine W- (negative rank sum)

W = min([Wpos,Wneg],[],2);  %calculate W-statistic for each case

bins = 0:(n*(n+1)/4-0.5);   %all possible integer W-values, # of "bins"
binlims = -0.5:(n*(n+1)/4); %min & max interval limits for each "bin"
h = histcounts(W,binlims);  %count # of occurances of each U-value
hprob = h/numel(W);     %probability for each "bin"
hcdf = cumsum(hprob);   %cumulative probability for <= each "bin"

Wcrit = floor(interp1(hcdf,bins,alpha));    %maximum W-value to reject null

%Plot PDF
figure
hold on
grid on
b1 = bar(bins,hprob);       %probability density function
title(sprintf('Probability Density, n=%i',n))
xlabel('U-Statistic')
ylabel('Probability')
b0 = bar(NaN(2,1));         %dummy bar plot for legend entries
b0.FaceColor = [1 0 0];     %red
legend(b0,'Critical W-Statistic','Location','northwest')

%Plot CDF
figure
hold on
grid on
b2 = bar(bins,hcdf);        %cumulative distribution function
title(sprintf('Cumulative Distribution, n=%i',n))
alphaline = plot([0,n*(n+1)/4-.5],ones(1,2)*alpha,'k--');
xlabel('U-Statistic')
ylabel('Cumulative Probability')
b0 = bar(NaN(2,1));         %dummy bar plot for legend entries
b0.FaceColor = [1 0 0];     %red
legend([b0,alphaline],{'Critical W-Stat','\alpha = 0.05'},...
    'Location','northwest')

%Highlighting critical W-statistic
if ~isnan(Wcrit)
    b1.FaceColor = 'flat';
    b1.CData(Wcrit+1,:) = [1 0 0];  %make critical W-statistic red
    b2.FaceColor = 'flat';
    b2.CData(Wcrit+1,:) = [1 0 0];  %make critical W-statistic red
end

