function [] = ClusteringTimeSeries(ERPs, sig_chans)
%% Run through Events to Calculate The Ward Linkagestructre
time_axis = ERPs.time_axis;

timepts = length(time_axis);
advance = 5;
window_length = 20;
frames = floor((timepts-window_length)/advance);

%% Output Waves from Linkage
maxes = ones(frames, 1);
mins = ones(frames,1);
means = ones(frames, 1);
time_val = ones(frames,1);

%% Loop through and calcualte linkage at different times
for k = 1:frames
    timewin = advance*k:(advance*k+window_length);
    center = floor(mean(timewin));
    
    link_Order = SampleSize_of_Centroids(ERPs, sig_chans, 50, timewin);
    
    time_val(k) = time_axis(center);
    maxes(k) = max(link_Order);
    mins(k) = min(link_Order);
    means(k) = mean(link_Order);
end

figure;
plot(time_val, maxes, time_val, mins, time_val, means)
legend('Maximum', 'Minimum', 'Mean')
ylabel('Euclidean Distance')
xlabel('Time from Location Tag Onset (ms)')
hold

plot([0 0], get(gca, 'ylim'),'k', 'LineWidth',2)
a = 1;

end
