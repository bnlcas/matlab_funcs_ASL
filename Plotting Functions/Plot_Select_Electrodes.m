function [] = Plot_Select_Electrodes(plot_electrodes, ERPs, ecog, comp_inds, colorlist, tags)
%% This function takes a select list of electrodes and Plots ERPs on the basis of the indecies of a 
% n x indecies matrix of comparisons
% This function is explicitly desiged for use in the generate ASL Plots
% function.

%% Establish Paramters
comp1 = tags{1};
comp2 = tags{2};
eq_param = tags{3};

%% Find Max and Min Bounds for Relevant Traces
eq_scale = true;
if eq_scale
    for i = 1:size(comp_inds,2)
       data = ecog(plot_electrodes,:,comp_inds(:,i));
       mean_data = mean(data,3);
       sem_data = std(data,[],3)/sqrt(size(data,3));
       lower_bounds(i) = min(min(mean_data - sem_data));
       upper_bounds(i) = max(max(mean_data + sem_data));
    end
    ybounds = [min(lower_bounds) max(upper_bounds)];
end

figure;
for i = 1:length(plot_electrodes)
    % Create a subplot in the squarest grid
    subplot(ceil(sqrt(length(plot_electrodes))),ceil(sqrt(length(plot_electrodes))),i)
    data = squeeze(ecog(plot_electrodes(i),:,:));
    time_axis = ERPs.time_axis;
    for j = 1:size(comp_inds,2)
        shadedErrorBar(time_axis,mean(data(:,comp_inds(:,j)),2), nansem(data(:,comp_inds(:,j)),2),{'color', colorlist(j,:)},1);
        hold on;
    end
    axis tight;
%             set(gca,'YLim',[minbound maxbound]...
%                 ,'XTickLabel',[],'YTickLabel',[],'YTick',[minbound maxbound]);
    if eq_scale
        ylim(ybounds);
    end
    line([0 0],get(gca,'YLim'),'Color','k');
    line(get(gca,'XLim'),[0 0],'Color','k');
    xlabel('Time From Onset (ms)')
    ylabel('High Gamma Intensity')
    title({['Plot of High Gamma in Ch ', num2str(plot_electrodes(i))]; [' for ', comp1, ' and ', comp2, ' Onset ERPs with equalized ', eq_param, ' position']})
end