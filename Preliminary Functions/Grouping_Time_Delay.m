function [] = Grouping_Time_Delay(ERPs)
%% Find the time delay between Grouped Lexical Event Onsets and the Onset of a subtag Event

tag = ERPs.annot.handshape; % Variable 
tag = strrep(tag, ' ', '');
tag = strrep(tag, 'changing', '');
is_tag = ~strcmpi(tag, '');
group_data_ind = find(is_tag & strcmpi(ERPs.alt_annot.group_tag, 'G'));


tag_delay = zeros(length(group_data_ind),1);
for i = 1:length(group_data_ind);
    tag_onset = ERPs.annot.start_ms(group_data_ind(i));
    
    found_sub = false;
    j = group_data_ind(i);
    while ~found_sub
        if is_tag(j);
            found_sub = true;
            subtag_onset = ERPs.annot.start_ms(j);
        end
        j = j+1;
    end
    tag_delay(i) = subtag_onset - tag_onset;
end

figure;
histogram(tag_delay)

a = 1;

end
    