function [StimOnTimes] = Get_LogStimTimes(MasterTable)
%%This function takes the MasterTable as an input and determines the
%%LogFile times at which stimuli are given

Stimulus = table2array(MasterTable(:,11));      % array listing of whether an action is Stimulus by the marker 'S' 
LogTimes = table2array(MasterTable(:,4));       % Array listing when each event occurs

StimOn = (strcmp(Stimulus, 'S'));               % boolean array: is ON iff the entry of its index in the MasterFile is a Stimulus 'S' event

StimOnTimes = LogTimes(StimOn) 

end
