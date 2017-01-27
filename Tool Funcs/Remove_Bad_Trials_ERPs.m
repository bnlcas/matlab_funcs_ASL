function [ERPs_Trim] = Remove_Bad_Trials_ERPs(ERPs)
%%This function creates a structure ERPs_Trim which copies the Data from
%%ERPs with rows for bad channels removed

ERPs_Trim = ERPs;
num_trails = length(ERPs.Blocks);
bad_trials = cell2mat(ERPs.BadTrials(:,2));


num_removed = 0;
for i = 1:length(bad_trials)
    delete_index = bad_trials(i) - num_removed;
    EPRs_Trim.ecog(:,:,delete_index) = [];
    EPRs_Trim.stims(delete_index) = [];
    ERPs_Trim.Blocks(delete_index) = [];
    ERPs_Trim.cond_num(delete_index) = [];
    ERPs_Trim.freq(delete_index) = [];
    ERPs_Trim.aoa(delete_index) = [];
    
    removed = removed + 1;
end

end