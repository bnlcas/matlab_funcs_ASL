function PlotECogGrid_ShadedPatces(ERPs, Means, is_sig_mat, varargin)
%% One off Function to shade in significant pvals
%% This function takes a input matricies of 256 x timepts containing
% some ststatisitic that applies over the ECoG grid channels over time,
% and plots the time evolution of this parameter.
%
% The tag means sets the function to take the mean and standard error
% of the inputs, otherwise it is assumed that the inputs are already in the
% form of some statistic.
%
% under the current settings the graph plots the one second prior to and
% the one second posterior to the start of the tagged event.
%
% The tag normalize, will normalize all plots to a maximum absolute value
% of 1. If empty, this feature is turned off.



%% Get significant channel plot:
is_sig_mat = ttest_priming_differences(varargin{1}, varargin{2}, 0.05, true);


num_plots = length(varargin);

Normalize = false;
if isempty(Normalize)
    Normalize = false;
end
Normalize = false;
if isempty(Means)
    Means = false;
end

%the first ecog_statistic is plotted in blue, while the second statistic is
%plotted in red
%The Boolean Tag Normalize is set to 1 iff the two plots should be
%normalized to the same axis


%normalize the data of statistic 1 and 2 such that the maximum displacement
%of the statistic is set to 1


%
timerange = 101:301;
%timerange = 301:801;
%timerange = 1:801;
%timerange = 101:401;
time_axis = ERPs.time_axis(timerange);
is_sig_mat = is_sig_mat(:,timerange);




if Means
%% Processing for Means
    %% Get Mean and SEM over trials of ECoG data
    trial_size = zeros(num_plots,1);
    for i = 1:num_plots
        trial_size(i) = size(varargin{i},3);
    end
    trial_size
    mean_mat = [];      % mean_mat is a 256xtimeptsx number of input ECoGs matrix
    sem_mat = [];       % standard error of mean for mean_mat
    for i = 1:num_plots
        mean_mat(:,:,i) = mean(varargin{i},3);
    end
    for i = 1:num_plots
        sem_mat(:,:,i) = nansem(varargin{i},3);
    end



    %% Clear Bad Channels:
    clear_bad_chan = true;
    if clear_bad_chan
        Bad_Channels = ERPs.BadChans{:,2};
%         Cleared_Chans = 1:256;
%         Cleared_Chans = Cleared_Chans(~SigChans);
%         Bad_Channels = [Bad_Channels, Cleared_Chans];
        for i = 1:num_plots
            mean_mat(Bad_Channels,:,:) = 0;
            sem_mat(Bad_Channels,:,:) = 0;
        end
    end
    
    %% Get axis scaling
    maxbounds = zeros(i,1);
    minbounds = zeros(i,1);
    for i = 1:num_plots
        maxbounds(i) = max(max(squeeze(mean_mat(:,timerange,i)+sem_mat(:,timerange,i))));
    end
    for i = 1:num_plots
        minbounds(i) = min(min(squeeze(mean_mat(:,timerange,i)-sem_mat(:,timerange,i))));
    end

    if Normalize
        for i = 1:num_plots
            norm = max(max(abs(mean_mat(:,timerange,i))));
            mean_mat(:,timerange,i) = mean_mat(:,timerange,i)./norm;
            sem_mat(:,timerange,i) = sem_mat(:,timerange,i)./norm;
            %reset the max and min window
            maxbounds(i) = max(max(mean_mat(:,timerange,i)+sem_mat(:,timerange,i)));
            minbounds(i) = min(min(mean_mat(:,timerange,i)-sem_mat(:,timerange,i)));
        end
    end

    % find the upper and lower bounds of the plotting windows:
    maxbound = max(maxbounds);
    minbound = min(minbounds);
    [minbound, maxbound]
    
    
    %% Plot Grid
    figure;
    hold on; 
    % colorlist = ['r' 'b' 'k' 'g' 'm' 'y' 'k']; % this list assigns colors to plots - no not plot more than 7 vars
    colorlist = [0.8 0 0; 0 0 0.8;0.8 0.8 0; 0 0.8 0.8; 0 0.8 0; 0.8 0 0.8];    % Plots of colors listed as Red, Blue, Yellow, Teal, Green, Purple
    grd = reshape(1:256,16,16)';
    for i = 1:256
        p = plotGridPosition(grd(i)); 
        subplot('position',p);
        %plots the time evolution of the ECoG statistic (say average)
        %for a given channel
        %shadedErrorBar(ERPs.time_axis(timerange),ecog_1(i,timerange), ecog_error_1(i,timerange),'b',1);
        %% Plot Shaded Patches for Channel
        patch_x = zeros(1,4);
        patch_y = [minbound maxbound maxbound minbound];
        sig_inds = find(squeeze(is_sig_mat(i,:)));
        inds_spacing = diff([0, sig_inds]);
        k = 1;
        while k < (length(inds_spacing)-1)
            if inds_spacing(k) > 1 & inds_spacing(k+1) == 1
                patch_x(1:2) = time_axis(inds_spacing(k));
                k = k +1;
                while k < (length(inds_spacing)-1) | inds_spacing(k+1) > 1
                    k = k+1;
                end
                patch_x(3:4) = time_axis(inds_spacing(k));
                
                consequtive_thresh = 4; % Demand that multiple significant events must happen
                if patch_x(3) > patch_x(1)+consequtive_thresh
                    patch(patch_x, patch_y,'black', 'FaceAlpha', 0.8, 'EdgeAlpha',0.8)
                end

            end
        end
        
        for j = 1:num_plots
            shadedErrorBar(time_axis,mean_mat(i,timerange,j), sem_mat(i,timerange,j),{'color', colorlist(j,:)},1);
        end
        %shadedErrorBar(ERPs.time_axis(timerange),ecog_2(i,timerange), ecog_error_2(i,timerange),'r',1);

        axis tight;
        set(gca,'YLim',[minbound maxbound]...
            ,'XTickLabel',[],'YTickLabel',[],'YTick',[minbound maxbound]);
        line([0 0],get(gca,'YLim'),'Color','k');
        %line([910 910],get(gca,'YLim'),'Color','b')
        line(get(gca,'XLim'),[0 0],'Color','k');
        text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim')) * 0.1)),num2str(i));

    end;

%% Plot for basic statisics
else
    %% Assemble Input Matrix
    stats_mat = [];      % mean_mat is a 256xtimeptsx number of input statsitics matrix
    for i = 1:num_plots
        stat_mat(:,:,i) = varargin{i};
    end
    
    %% Clear Bad Channels:
    clear_bad_chan = true;
    if clear_bad_chan
        Bad_Channels = ERPs.BadChans{:,2};
        for i = 1:num_plots
            stat_mat(Bad_Channels,:,:) = 0;
        end
    end
    
    %% Get axis scaling
    maxbounds = zeros(i,1);
    minbounds = zeros(i,1);
    for i = 1:num_plots
        maxbounds(i) = nanmax(nanmax(squeeze(stat_mat(:,timerange,i))));
    end
    for i = 1:num_plots
        minbounds(i) = nanmin(nanmin(squeeze(stat_mat(:,timerange,i))));
    end

   
    if Normalize
        for i = 1:num_plots
            norm = max(max(abs(stat_mat(:,timerange,i))));
            stat_mat(:,timerange,i) = stat_mat(:,timerange,i)./norm;
            %reset the max and min window
            maxbounds(i) = max(max(squeeze(stat_mat(:,timerange,i))));
            minbounds(i) = min(min(squeeze(stat_mat(:,timerange,i))));
        end
    end

    % find the upper and lower bounds of the plotting windows:
    maxbound = max(maxbounds)
    minbound = min(minbounds)
    
    
    %% plot grid:
    figure;
%    colorlist = ['r' 'b' 'c' 'g' 'm' 'y' 'k']; % this list assigns colors to plots - no not plot more than 7 vars
    colorlist = [0.8 0 0; 0 0 0.8;0.8 0.8 0; 0 0.8 0.8; 0 0.8 0];
   
    
  
    grd = reshape(1:256,16,16)';
    for i = 1:256
        p = plotGridPosition(grd(i)); 
        subplot('position',p);
        hold on; 
        %plots the time evolution of the ECoG statistic (say average)
        %for a given channel
        %shadedErrorBar(ERPs.time_axis(timerange),ecog_1(i,timerange), ecog_error_1(i,timerange),'b',1);
        %% Plot Shaded Patches for Channel
        patch_x = zeros(1,4);
        patch_y = [minbound maxbound maxbound minbound];
        sig_inds = find(squeeze(is_sig_mat(i,:)));
        inds_spacing = diff([0, sig_inds]);
        k = 1;
        while k < (length(inds_spacing)-1)
            if inds_spacing(k) > 1 & inds_spacing(k+1) == 1
                patch_x(1:2) = time_axis(sig_inds(k));
                if k < (length(inds_spacing)-1)
                    while k < (length(inds_spacing)-1) & inds_spacing(k+1) == 1
                        k = k+1;
                    end
                    if k+1 ==1
                        patch_x(3:4) = time_axis(sig_inds(k+1));
                    else
                        patch_x(3:4) = time_axis(sig_inds(k));
                    end
                patch(patch_x, patch_y,'black', 'FaceAlpha', 0.2, 'EdgeAlpha',0.2)
                end
            end
            k = k+1;
        end
        
        
        for j = 1:num_plots
            plot(time_axis,stat_mat(i,timerange,j),'Color', colorlist(j,:));
            %plot(time_axis,stat_mat(i,timerange,j),colorlist(j));
            hold on; 
        end
        axis tight;
        set(gca,'YLim',[minbound maxbound]...
            ,'XTickLabel',[],'YTickLabel',[],'YTick',[minbound maxbound]);
        line([0 0],get(gca,'YLim'),'Color','k'); 
        line(get(gca,'XLim'),[0 0],'Color','k');
        text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim')) * 0.1)),num2str(i));

    end;
end

a = 1;