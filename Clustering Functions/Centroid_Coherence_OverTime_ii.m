function [] = Centroid_Coherence_OverTime_ii(ERPs, sig_chans)
%% Plots the Mean distance of a cluster from its centroid and plot this over time



class_thresh = 10;

isgood = is_good_trial(ERPs);
%isgood(:) = true;
labels = ERPs.annot.handshape(isgood);
categories = unique(labels);
categories = categories(~strcmpi(categories,''));
cat_count = get_category_size(labels);
categories(cat_count < class_thresh) = [];

% categories(strcmp(categories,'neutral')) = [];
% categories(strcmp(categories,'changing')) = [];

centroids = [];

%% Loop through and cal norm distance between two centroids overtime:
time_axis = ERPs.time_axis;

timepts = length(time_axis);
advance = 5;
window_length = 100;
frames = floor((timepts-window_length)/advance);

%% Output Waves from Linkage
mean_dist_cent = [];
time_val = ones(frames,1);

%% Loop through and calcualte linkage at different times
for k = 1:frames
    twin = advance*k:(advance*k+window_length);
    centroids = [];
    for i = 1:length(categories)
        centroids(i,:) = mean(mean(ERPs.ecog(sig_chans,twin,strcmpi(labels, categories(i))),2),3);
    end
    for i = 1:length(categories)
        rep_cent = repmat(centroids(i,:),sum(strcmpi(labels, categories(i))),1);

        cluster_dists = pdist2(squeeze(mean(ERPs.ecog(sig_chans,twin,strcmpi(labels, categories(i))),2))',rep_cent );
%        mean_dist_cent(k,i) = trace(cluster_dists)./size(cluster_dists,1);
        diag_els = max((eye(size(cluster_dists)).*cluster_dists)); % only rows are the distance of point the centroid
        %Alt: diag_els = cluster_dists(:,1);  Above is more general to fire
        %distances between two sets
        mean_dist_cent(k,i) = mean(diag_els);
        std_dist_cent(k,i)=std(diag_els);
        stderr_dist_cent(k,i) = std_dist_cent(k,i)/mean_dist_cent(k,i);
    end
    time_val(k) = time_axis(floor(mean(twin)));
end
time_val = repmat(time_val,1,length(categories));

figure;
plot(time_val, mean_dist_cent)
ylabel('Euclidean Distance')
xlabel('Time from Location Tag Onset (ms)')
hold
plot([0 0], get(gca, 'ylim'),'k', 'LineWidth',2)
title('Plot of Distance Between Centroid of Two Closest Sign Locations Over Time')
a = 1;


end