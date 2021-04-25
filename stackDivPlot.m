function [b,stats] = stackDivPlot(varargin)
%Creates a Stacked Diverging Bar Plot in a new Figure window based on
%survey responses to multiple questions that use identical rating scales
%
%   Written by:     Benjamin Moidel
%   Last Updated:   April 15, 2021
%
%Usage:
%   stackDivPlot(M)
%   stackDivPlot(M,rmax)
%   stackDivPlot(M,rmax,rmid)
%   stackDivPlot(ax,____)
%   b = stackDivPlot(____);
%   [b,stats] = stackDivPlot(____);
%Inputs:
%   ax = (optional) target axes for plotting; otherwise plots in new figure
%   M = matrix of integer survey responses with 1 as lowest score
%       rows = each question
%       columns = each response
%       *Note: if a certain question has less than the maximum number of
%        responses, pad the rest of the row with NaN*
%   rmax = (optional) maximum survey score, assuming scale is from 1 to
%          rmax; otherwise assumes highest response is maximum score
%   rmid = (optional) median score on scale, values must be greater than
%          this value to be "positive" (right of center line), otherwise if
%          <= rmid will be considered "negative" (left of center line)
%Outputs:
%   b = (optional) bar graph object
%   stats = (optional) structure containing statistics of the data:
%       rfreq = response frequency of sample mode (per question)
%       sdev = sample standard deviation (per question)
%       smean = sample mean (per question)
%       smed = sample median (per question)
%       smode = sample mode (per question)
%       sneut = neutral rating on the scale (ALL questions)
%       svar = sample variance (per question)

%Get information about function inputs
if isnumeric(varargin{1})   %plot in new figure window
    figure
    ax = gca;
    M = varargin{1};
    if nargin<2     %find max response if not given as an input
        rmax = max(M,[],'all'); %highest overall response
    else
        rmax = varargin{2};
    end
    if nargin<3     %find neutral rating if not given as an input
        rmid = floor(rmax/2);       %column index of neutral rating
        stats.sneut = (rmax+1)/2;   %neutral rating, negative if <= sneut
    else
        rmid = floor(varargin{3});  %column index of neutral rating
        stats.sneut = rmid;         %neutral rating, negative if <= sneut
    end
elseif isobject(varargin{1})    %axes was specified to plot on
    ax = varargin{1};
    M = varargin{2};
    if nargin<3     %find max response if not given as an input
        rmax = max(M,[],'all'); %highest overall response
    else
        rmax = varargin{3};
    end
    if nargin<4     %find neutral rating if not given as an input
        rmid = floor(rmax/2);       %column index of neutral rating
        stats.sneut = (rmax+1)/2;   %neutral rating, negative if <= sneut
    else
        rmid = floor(varargin{4});  %column index of neutral rating
        stats.sneut = rmid;         %neutral rating, negative if <= sneut
    end
else
    error(['Invalid data type. ',...
        'First argument must be an axes handle or numeric matrix.'])
end

%Get information about data matrix
[numQ,numR] = size(M);  %get # of questions, # of responses

%Count each distinct response to each question
F = zeros(numQ,rmax);   %preallocation
for i = 1:numQ
    for j = 1:rmax
        F(i,j) = sum(M(i,:)==j);  %frequency matrix
    end
end

%Statistics
stats.smed = median(M,2,'omitnan');     %median response per question
stats.smean = mean(M,2,'omitnan');      %mean response per question
[stats.smode,stats.rfreq] = mode(M,'all');  %mode response and # of times
    %it occurs per question; if more than 1 response is the mode, this
    %returns the smallest value!
stats.svar = var(M,0,2,'omitnan');  %sample variance per question
stats.sdev = sqrt(stats.svar);  %sample standard deviation per question

%Rearranging Frequency Matrix
newind = [rmid:-1:1,(rmid+1):rmax]; %new column index order (for plotting)
F = [-F(:,1:rmid),F(:,(rmid+1):rmax)];  %make bad responses negative
F = F(:,newind);    %reorder columns so bars stack correctly

%Colors for Plotting (blue = good, red = bad)
c = flipud(interp1(1:64,jet(64),linspace(1,64,rmax)));

%Plotting Bar Graph
b = barh(ax,F,'stacked');

%Formatting Plot
ax.XLim = [-numR,numR]; %place neutral response in center of plot
ax.YDir = 'reverse';    %put first question at top, last at bottom
lgdcell{rmax} = '';     %preallocation
for r = 1:rmax
    b(r).FaceColor = c(newind(r),:);
    lgdcell{rmax-r+1} = sprintf('%3i',r);
end

%Creating Legend
lgd = legend(b(fliplr(newind)),lgdcell,'Location','west');
title(lgd,'Responses')

%Creating Colorbar
% colormap(c)
% cb = colorbar;
% cb.Ticks = ((1:rmax)-0.5)/rmax;
% cb.TickLabels = 1:rmax;
% cb.TickDirection = 'out';
% ylabel(cb,'Responses')

if nargout<1    %clear b from memory if no outputs requested
    clearvars b
end