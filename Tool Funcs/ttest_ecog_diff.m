function [pvals] = ttest_ecog_diff(ecog1, ecog2)
%% Function takes two 256x401xtrials ecog matricies and along with a pvalue threshold
% and returns a 256x401 boolean matrix that is true when a ch & time point
% has a statisticall significant difference in two way ttest.

chans = size(ecog1,1);
timepts = size(ecog1,2);
pvals = zeros(chans, timepts);
for i = 1:chans
    for j = 1:timepts
        dat1 = squeeze(ecog1(i,j,:));
        dat2 = squeeze(ecog2(i,j,:));
        
        [~,pvals(i,j)] = ttest2(dat1, dat2);
    end
end
pmat = pvals;
           

