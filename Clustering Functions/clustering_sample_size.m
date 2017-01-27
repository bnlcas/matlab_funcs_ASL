function [] = clustering_sample_size(ERPs, sig_chan)
%% Plot the effect of reducing sample size on Clustering ordering
meandist = [];
class_size = 50;
center = ClusterOneClass(ERPs, class_size, sig_chan, class_size, false);
for i = 1:2:class_size 
    i
    % order = ClusterOneClass(ERPs, class_size, sig_chan, i, true);
    samplings = 500;
    order = zeros(samplings, length(center));
    for j = 1:samplings;
        order(j,:) = ClusterOneClass(ERPs, class_size, sig_chan, i, true);
    end
    % center = mode(order,1);
    stdev = 0;
    for j = 1:samplings
        stdev = stdev + pdist([order(j,:); center], 'hamming');
    end
    stdev = stdev/samplings;
    meandist = [meandist, stdev];
end
figure; plot(1:2:class_size, meandist, 'rx')
a = 1;


