function [is_good] = is_good_trial_block(ERPs)
%% This function returns a Trails length boolean array that is true
% for all sub-erps in a stimulus trial if all are good
% and false for all sub-erps in stimulus section is any suberps is bad

%% Get a List of Good Trials
bad_trials = cell2mat(ERPs.BadTrials(:,2));
is_good = true(size(ERPs.ecog,3),1);
%is_good(bad_trials) = false;

trial_nums = ERPs.annot.trl;
for i = 1:length(bad_trials)
    bad_block = strcmpi(trial_nums, trial_nums(bad_trials(i)));
    is_good(bad_block) = false;
end

end