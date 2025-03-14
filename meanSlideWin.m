
%%% Adapted from Hardwick et al. 2019
function [ts] = meanSlideWin(x, y, winLen, timeMax, timeMin, nBoot, wantPlot, color, lineWidth, varargin)

%%% x           = vector of the variable to slide along
%%% y           = vector of the variable to be averaged
%%% timeMax     = max value of y to be covered by the window
%%% timeMin     = min value of y to be covered by the window
%%% nboot       = number of iterations for the bootstrapped confidence interval
%%% wantPlot    = if plot it will plot the result

%%% Set defaults
if nargin < 3; winLen = 100; timeMax = max(x); timeMin = 1; nboot = 1000; wantPlot = 0;
elseif nargin < 4; timeMax = max(x); timeMin = 1; nboot = 1000; wantPlot = 0; color = [0.2588 0.2863 0.2863]; linewidth = 3;
elseif nargin < 5; timeMin = 1; nboot = 1000; wantPlot = 0; color = [0.2588 0.2863 0.2863]; linewidth = 3;
elseif nargin < 6; nBoot = 1000; wantPlot = 0; color = [0.2588 0.2863 0.2863]; linewidth = 3;
elseif nargin < 7; wantPlot = 0; color = [0.2588 0.2863 0.2863]; linewidth = 3;
elseif nargin < 8; color = [0.2588 0.2863 0.2863]; linewidth = 3;
elseif nargin < 9; lineWidth = 3;
end


%%% Slide the windows
ts = [];
for i = timeMin : timeMax
    %%% Set the center of the window at 1
    lower   = i - (winLen/2);               % Lower half of the window
    upper   = i + (winLen/2);               % Upper half
    points  = find(x <= upper & x >= lower);    % points of x contained in the window of y
    sample  = y(points);                    
    
    %%% if the window contains points, average them and bootrstrap CIs
    if length(sample) > 1
        average = mean(sample);
        ci = bootci(nBoot, @mean, sample);
        ts = [ts; [ci(1), average, ci(2)]];
    else
        ts = [ts; [NaN NaN NaN]];
    end
    
end

%%% if you want plot the result
if wantPlot
    tsPlot = ts';
    time = 1:length(tsPlot);
    patch([time, fliplr(time)], [tsPlot(1, :), fliplr(tsPlot(3, :))], color, 'FaceAlpha',0.3, 'EdgeColor','none', 'HandleVisibility', 'off')
    hold on
    plot(tsPlot(2, :), 'Color', color, 'LineWidth', lineWidth)
    xlim([1, length(tsPlot)])
end

end




