function RecTime = LogFileTimeConvert(MasterTable, LogTime, RowNum)
%%This function takes the MasterTable of annotated ECoG data is converts
%%the time value of a RowNum from and a point in LogTime into the scaling of the 
%%RecTime in seconds.
%%This function assumes four blocks of testing, labled in the MasterTable
%%as Log File : 5.2, 5.35, 5.4, 5.5

Col_BlockLabel = 3;     % Column in Masterfile that contains info on what block an event occured in.

Sec_num = table2array(MasterTable(:,Col_BlockLabel));     % array of section number of each event

if Sec_num(RowNum) == 5.2
    RecTime = 1.0002*(LogTime./10000) + 2.5067;
elseif Sec_num(RowNum) == 5.35
    RecTime = 1.0*(LogTime./10000) + 15.9492;
elseif Sec_num(RowNum) == 5.4
    RecTime = 1.0*(LogTime./10000) + 14.9606;
elseif Sec_num(RowNum) == 5.5
    RecTime = 1.0001*(LogTime./10000) + 9.5850;
end

end