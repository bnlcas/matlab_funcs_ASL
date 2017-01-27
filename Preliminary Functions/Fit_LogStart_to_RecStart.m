function [] = Fit_LogStart_to_RecStart(LogFile_StartTimes, Rec_StartTimes)
%%This function loops through the 4 sectors of stimuli and plots linear
%%regression models for each region to convert the time recorded on the
%%logfile (divied by 10^4) into the recorded stimulus times.
pts_per_sample = 79;     % magic constant representing the 79 trials per block
ch = [1 2 3 4];          % magic constant representing the 4 blocks
for i = 1:length(ch)
    start_bound = 1+pts_per_sample*(ch(i)-1);
    end_bound = pts_per_sample*ch(i);
    Y = Rec_StartTimes(start_bound:end_bound)';
    X = [ones(1, pts_per_sample); LogFile_StartTimes(start_bound:end_bound)./(10^4)]';
    [fit_co,~,~,~,stats] = regress(Y, X);
    fit_co
    r2=stats(1)
end
end