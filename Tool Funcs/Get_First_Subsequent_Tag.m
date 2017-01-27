function [First_Secondary_Tag, Filled_Secondary_Tag] = Get_First_Subsequent_Tag(Primary_Tag, Secondary_Tag)
%% Takes a Boolean Vector for the ERP.annot structures and returns a related
% Boolean vector which is on for the first index that has Secondary Param
% after each instance of Primary Tag
%% Example:
% y = Get_First_Subsequent_Tag(is_stim_onset, is_lexical_onset);
% returns the subset of is_lexical onset such that each is the first
% lexical onset tag to occur after a stimulus - may have few entries than
% primary tag in the event of no lexical resposne

% The output Filled_Secondary_Tag is on for all primary tags that had at
% least one subsequent secondary tag.

First_Secondary_Tag = false(size(Primary_Tag));
Filled_Secondary_Tag = Primary_Tag;

primary_inds = find(Primary_Tag);
secondary_inds = find(Secondary_Tag);

%% Find the indecies of secondary tag between each clump of primary tags
j = 1;
for i = 1:(length(primary_inds)-1)
    possible_inds = find(secondary_inds >= primary_inds(i) & secondary_inds <= primary_inds(i+1),1);
    if ~isempty(possible_inds)
        first_secondary_inds(j) = secondary_inds(possible_inds);
        j = j+1;
    else
        Filled_Secondary_Tag(primary_inds(i)) = false;
    end
end
% End Case:
possible_inds = find(secondary_inds >= primary_inds(end),1);
if ~isempty(possible_inds)
    first_secondary_inds(j) = secondary_inds(possible_inds);
else
    Filled_Secondary_Tag(primary_inds(end)) = false;
end

First_Secondary_Tag(first_secondary_inds) = true;

end