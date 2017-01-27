function [] = fit_line_wparams(x,y)
%% Plots Y vs X with a fit line fitting parameters
scatter_marker = 'rx';
regress_color = 'b';

if size(x,1)==1
    x = x';
end
if size(y,1) == 1
   y = y';
end

figure; plot(x,y,scatter_marker)
a = lsline;
a.Color = regress_color;

[b,~,~,~,stats] = regress(y,[ones(length(x),1),x]);

dim = [.18 .5 .4 0.4];
str = {['y = ' num2str(b(2),3) '*x + ' num2str(b(1),3)];...
    ['R^2 = ' num2str(stats(1),3)];...
    ['pval = ' num2str(stats(3),3)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
