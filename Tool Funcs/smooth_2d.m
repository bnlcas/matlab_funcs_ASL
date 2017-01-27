function [out] = smooth_2d(input, smoothing, smooth_dimension);
%% Takes an 2d matrix and smooths the data in one dimension
% inputs:
% input: n x m matrix
% 
% smoothing: degree of smoothing
% 
% smooth_dimension: the axis to smooth the data on (1, 2)
out = zeros(size(input));

if smooth_dimension == 1
    for i = 1:size(input,2)
        out(:,i) = smooth(input(:,i), smoothing);
    end
elseif smooth_dimension == 2
    for i = 1:size(input,1)
        out(i,:) = smooth(input(i,:), smoothing);
    end
end



end