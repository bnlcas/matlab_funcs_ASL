function [pvals] = ttest_over_time(ERPs, sig_chans, ecog1, ecog2)
% This function runs a two way T-Test on the significant Channels of a the
% ECoG on the Trails from Label1 and Label2 over the timecourse of the ERPs
% 
%Data_Tag = is_good_trial(ERPs) & strcmpi(ERPs.annot.filledLexTrans,'lexical');

%ecog1 = ERPs.ecog(sig_chans, :, Data_Tag & label1);
%ecog2 = ERPs.ecog(sig_chans, :,  Data_Tag & label2);

normalize = true;
if normalize
    ecog1 = ecog1./(max(abs(ecog1(:))));
    ecog2 = ecog2./(max(abs(ecog2(:))));
end

pvals = zeros(length(ERPs.time_axis),size(ecog1,1));
for i = 1:length(ERPs.time_axis)
    data1 = squeeze(ecog1(:,i,:))';
    data2 = squeeze(ecog2(:,i,:))';
    [~,p] = ttest2(data1, data2,'Vartype','unequal');
    pvals(i,:) = p;
end
a =1;
end