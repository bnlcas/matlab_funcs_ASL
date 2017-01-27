function [is_grouped] = is_grouped_row(ERPs)
%% This function takes a grouped annotation and returns a boolean value 
% depending on whether a given row contains some set of annotations
% then reduces this list to only include one from each set of consequtive
% redundant lists
% Non-Lingusitic ERPs will be excluded from groupings


annot = ERPs.filled_change_annot;

%% Which annotations will be used
use_hs = true;
use_loc = true;
use_intmov = false

relevant_tags = [];
if use_hs
    relevant_tags = [relevant_tags, annot.handshape];
end
if use_loc
    relevant_tags = [relevant_tags, annot.loc];
end
if use_intmov
    relevant_tags = [relevant_tags, annot.intMov];
end

%% Loop through and find rows with only non-reject tags
is_grouped = true(size(relevant_tags,1),1);
ungroup_lables = {'changing', ''};      % annotations that will disqualify a row
row_size = length(relevant_tags(1,:));  % Determine if
for i = 1:size(relevant_tags,1)
    row_data = relevant_tags(i,:);
    for k = 1:length(ungroup_lables)
        if sum(~strcmpi(row_data, ungroup_lables(k))) < row_size;
            is_grouped(i) = false;
        end
    end

end

%% Remove non-linguistic Gestures
is_grouped = is_grouped & strcmpi(ERPs.annot.filledLexTrans,'lexical');

%% Remove Redundancy
group_inds = find(is_grouped);

gap = [0; group_inds(2:end) - group_inds(1:(end-1))];
is_duplicate = (gap == 1);
dup_inds = group_inds(is_duplicate);
for i = 1:length(dup_inds)
    predecesor = mat2str(cell2mat(relevant_tags(dup_inds(i)-1,:)));
    element = mat2str(cell2mat(relevant_tags(dup_inds(i),:)));
    if strcmpi(element, predecesor)
        is_grouped(dup_inds(i)) = false;
    end
end
    
a = 1;