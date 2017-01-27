function [eq_inds] = Equalize_Tag_Sizes(varargin)
%% This function takes m (n x 1) boolean vectors denoting some information (neutral location vs face location)
% and returns a (n x m) matrix, where thes columns are subsets of the
% inputs, and each column has the same number of On Taggings as the others
% the down sampling to equalize the subset sizes is randomized.

num_inputs = length(varargin);
size_inputs = length(varargin{1});   % assumes input is non-empty...
eq_inds = false(size_inputs, num_inputs);

%% Set the size of subsets
in_set_sizes = zeros(num_inputs,1);
for i = 1:num_inputs
    data = varargin{i};
    in_set_sizes(i) = sum(data);
end
min_ds = 0; % extra down sampling to avoid biasing on the smallest class
out_set_size = min(in_set_sizes) - min_ds;


ind1 = [];
for i = 1:num_inputs
    data = varargin{i};
    set_inds = find(data);   
    shuffle = randperm(length(set_inds));   % shuffle on the indecies to be used
    set_inds = set_inds(shuffle);
    eq_inds(set_inds(1:out_set_size),i) = true;
    a=1;
end



