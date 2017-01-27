function [is_good] = is_good_chan(ERPs)
%% This function returns a 256 lenght boolean array that is true
% when the channel is not a bad channel

trials = size(ERPs.BadChans,1); % number of trials in data
bad_trials = [];
for i = 1:trials
    bad_trials = [bad_trials ERPs.BadChans{i,2}];
end
bad_trials = unique(bad_trials);

is_good = true(size(ERPs.ecog,1),1);
is_good(bad_trials) = false;

end