function [] = make_color_bar(data, cmap)
%% Generates Color Bars for Plots
% This function takes the Data used for a plot along with the colormap used
% for the plot and then creates a scaled colorbar using the imagesc
% funciton in matlab.
data_min = min(data(:));
data_max = max(data(:));

span = data_max - data_min;
pts = 100;

colorspan = (-1*data_min):span/pts:(1*data_max);
colorspan = colorspan(end:-1:1)';

figure; imagesc(colorspan);
colormap(cmap);
ax = gca;
ax.XTick = [];
colorbar

% title('Data Colorbar')