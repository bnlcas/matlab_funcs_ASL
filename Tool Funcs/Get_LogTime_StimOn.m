function [LogTime_StimOn] = Get_LogTime_StimOn(MasterTable)
%%This function returns a boolean array which is ON iff the entry for that
%%index in the MasterFile is a Stimulus 'S' event


%NOTE THIS FUNCTION IS NO LONGER CALLED: DELETE

Stimulus = table2array(MasterTable(:,11));
LogTime_StimOn = (strcmp(Stimulus, 'S'));

end