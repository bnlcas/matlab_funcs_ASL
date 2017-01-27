function [out] = smooth_3d(input, smoothing, smooth_dimension);
%% Takes an 3d matrix and smooths the data in one dimension
% inputs:
% input: n x m x p matrix
% 
% smoothing: degree of smoothing
% 
% smooth_dimension: the axis to smooth the data on (1, 2, 3)
out = zeros(size(input));

if smooth_dimension == 1
    for i = 1:size(input,2)
        for j = 1:size(input,3)
            out(:,i,j) = smooth(squeeze(input(:,i,j)), smoothing);
        end
    end
elseif smooth_dimension == 2
    for i = 1:size(input,1)
        for j = 1:size(input,3)
            out(i,:,j) = smooth(squeeze(input(i,:,j)), smoothing);
        end
    end
elseif smooth_dimension == 3
    for i = 1:size(input,1)
        for j = 1:size(input,2)
            out(i,j,:) = smooth(squeeze(input(i,j,:)), smoothing);
        end
    end
end




end
