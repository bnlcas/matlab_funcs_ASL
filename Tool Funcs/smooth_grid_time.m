function [out] = smooth_grid_time(input, smoothing);
%% Runs Smoothing on the time axis of the grid
out = input;
time_dimension = 2;
for i = 1:size(input,1)
    out(i,:) = smooth(input(i,:),smoothing);
end

end