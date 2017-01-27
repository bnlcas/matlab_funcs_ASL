function [] = DownSampling_DistanceVar(ERPs, sig_chans)
%% Calculate the average spread of the distance matrix (as a vector) vs Sample Size

Max_Sample_Size = 50;
Sampling_Points = 2:2:Max_Sample_Size;
iterations = 1000;

Distance_Variance = zeros(length(Sampling_Points),1);
test = SampleSize_of_Centroids(ERPs, sig_chans, 50, 1);
%% Loop through and calcualte linkage at different times
for k = Sampling_Points
    distances = zeros(iterations, length(test));
    for i = 1:iterations
        distances(i,:) = SampleSize_of_Centroids(ERPs, sig_chans, 50, k);
    end
%     Distance_Variance(k) = sum(sum(gsubtract(distances, mean(distances)).^2));
%    Distance_Variance(k) = sum(var(distances))./size(distances,2);
    
    mean_dist = mean(distances);    
    d_squared = 0;
    for i = 1:iterations;
        d_squared = d_squared + (distances(i,:) - mean_dist)*(distances(i,:) - mean_dist)';
    end
%    Distance_Variance(k) = sqrt(d_squared)/iterations;
    Distance_Variance(k) = (sqrt(d_squared)/iterations)./sqrt(mean_dist*mean_dist');
end

figure;
plot(Sampling_Points, Distance_Variance(Sampling_Points))
% legend('Maximum', 'Minimum', 'Mean')
ylabel('Standard Deviation of Distance between Centroids')
xlabel('Samples Size')

title('Plot of the Variance of Euclidean Distance Among Location Classifier Centroids vs. Sample Size')
a = 1;

end
