function [] = create_movie_time_bar(timeval)
%% This function is made for the brain plot movies, however it could be
% adpoted for use in other cinematic displays.
axes('Position', [0.1, 0.05, 0.8 0.05])

%% Color and Timing Settings:

%% Real-Psuedo or Trans Lex
trans_lex = false;
real_pseudo = true;

%% TransLex
if trans_lex
    color_1 =  [0 0 1];
    color_2 = [1 0 0];

    duration_1 = 262; % milliseconds) duration of transitional translex
    duration_2 = 277; % duration of lexical
    %% Real Values shifted over 5 for display purposes true are 267 and 272...
 
    stimulus_time_1 = -10e+03; % reaction time of real signs
    stimulus_time_2 = -100.5453e+03; % reaction time of pseudo signs

end


%% Real Pseudo
if real_pseudo
    color_1 =  [0     1     1];
    color_2 = [0.8192    0.3851    0.5000];

    duration_1 = 409; % duration ofmilliseconds of duplicate real signs
    duration_2 = 372; % duration of duplicate pseudo signs

    stimulus_time_1 = -1.4281e+03; % reaction time of real signs
    stimulus_time_2 = -1.6396e+03; % reaction time of pseudo signs
end

xrange = [-2000, 2000]; yrange = [-1 1];

%% Plot Current time point

plot([timeval, timeval], yrange, 'r','LineWidth',2.5)
hold on;

%% Plot Duration of Events
plot([0 0], yrange,'k','LineWidth',2) % Plot StartPt
plot([duration_1 duration_1], yrange,'Color', color_1,'LineWidth',1.5)
plot([duration_2 duration_2], yrange,'Color', color_2,'LineWidth',1.5)

%% Plot Time To Stimulus
plot([stimulus_time_1, stimulus_time_1], yrange, 'Color', color_1,'LineWidth',1.5)
plot([stimulus_time_2, stimulus_time_2], yrange, 'Color', color_2,'LineWidth',1.5)
%plot([stimulus_time stimulus_time], yrange, 'b')


ax = gca;
ax.YTick = [];
[relevant_pts,order] = sort([stimulus_time_2, duration_2, 0,]); % stimulus_time]); %, timeval]);
tags = {'        Mean Stimulus Onset',   '               Mean Duration', 'Onset'}; %'Mean Stimulus Onset'}; %, {'';'time'}};
ax.XTick = relevant_pts;
ax.XTickLabel = tags(order);
ax.TickLength = [0 0]; % Clear Tick Boxes


xlim(xrange)
ylim(yrange)
