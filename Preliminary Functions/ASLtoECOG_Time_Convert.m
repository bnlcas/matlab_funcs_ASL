function [out_time] = ASLtoECOG_Time_Convert(annot, in_time)


Block_Label = unique(annot.log);            % Array of the tags for each block of recording in the master file
Slopes = [1.0002 1 1 1.0001];               % Slope of scaling 
Offsets = [2.5067 15.9492 14.9606 9.5850];  % Offset of timing

out_time = zeros(size(in_time));
for i = 1:length(Block_Label) 
    Is_Section = strcmpi(annot.log,Block_Label{i}); %(annot.log == Block_Label(i));
    out_time = out_time + Slopes(i).*((in_time.*double(Is_Section))./10000)+Offsets(i).*double(Is_Section);
end
