function [] = Centroid_Coherence_OverTime(ERPs, sig_chans)
%% Finds the Two Closest Categories at the Peak and Traces their TimeEvo
Max_Sample_Size = 50;
iterations = 1000;

%% Find the closest two Categories at peak
twin = 191:211;
class_thresh = 50;

isgood = is_good_trial(ERPs);
%isgood(:) = true;
labels = ERPs.annot.loc(isgood);
categories = unique(labels);
categories = categories(~strcmpi(categories,''));
cat_count = get_category_size(labels);
categories(cat_count < class_thresh) = [];

% categories(strcmp(categories,'neutral')) = [];
% categories(strcmp(categories,'changing')) = [];

centroids = [];
for i = 1:length(categories)
    centroids(i,:) = mean(mean(ERPs.ecog(sig_chans,twin,strcmpi(labels, categories(i))),2),3);
end

distances = pdist(centroids);
% [~,mincoord] = min(distances);
% distances(:) = 0;
% distances(mincoord) =1;
% distances = squareform(distances);
% [rowmin,colmin] = find(distances,1);
% 
% categories = categories([rowmin, colmin]);


%% Loop through and cal norm distance between two centroids overtime:
time_axis = ERPs.time_axis;

timepts = length(time_axis);
advance = 5;
window_length = 5;
frames = floor((timepts-window_length)/advance);

%% Output Waves from Linkage
cent_dist = ones(frames, length(distances));
time_val = ones(frames,1);

%% Loop through and calcualte linkage at different times
for k = 1:frames
    twin = advance*k:(advance*k+window_length);
    centroids = [];
    for i = 1:length(categories)
        centroids(i,:) = mean(mean(ERPs.ecog(sig_chans,twin,strcmpi(labels, categories(i))),2),3);
    end
    cent_dist(k,:) = pdist(centroids);
    time_val(k) = time_axis(floor(mean(twin)));
end
time_val = repmat(time_val,1,size(cent_dist,2));

figure;
plot(time_val, cent_dist)
ylabel('Euclidean Distance')
xlabel('Time from Location Tag Onset (ms)')
hold
plot([0 0], get(gca, 'ylim'),'k', 'LineWidth',2)
title('Plot of Distance Between Location Centroids Over Time (No Bad Trials')a = 1;

end