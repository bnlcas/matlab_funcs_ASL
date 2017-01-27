function out_array = interleave(varargin)
%% This function takes m 1xm arrays and interleaves them into a single array
% Ex:
% dat1 = [1 2 3 4 5];
% dat2 = [10 20 30 40 50];
% out = interleaves(dat1, dat2);
% out will be [1 10 2 20 3 30 4 40 5 50];
%
arrays_combined = [];
for k = 1:length(varargin)
    tmp = varargin{k};
    arrays_combined = [arrays_combined tmp(:)'];
end

out_array = arrays_combined(:);
end