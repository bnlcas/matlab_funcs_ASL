function [loc_eq_HS_compare_Table] = EqualizeHSCompTable(ERPs)
%% Takes class of Data (ERPs.grouped_annot.loc) and creates a table
% Where each column is a type within that Colunm
hs_data = ERPs.grouped_annot.hs;
Data_Tag = is_good_trial(ERPs) & strcmpi(ERPs.annot.filledLexTrans,'lexical');

%% GetList of Locations for column rows:
handshapes = unique(hs_data(Data_Tag));
handshapes(strcmpi(handshapes, '')) = [];
handshapes(strcmpi(handshapes,'changing')) = [];
handshape(strcmpi(loc_names, 'lax')) = [];

hs_size = zeros(length(loc_names),1);
% Find biggest_hs
for i = 1:length(loc_names)
    hs_size(i) = sum(strcmpi(ERPs.grouped_annot.handshape(Data_Tag), handshape(i)));
end
[~, order] = sort(hs_size,'descend');
handshape = handshape(order(1:10)); % takes the 10 largest categories



%% Assemble List of Handshapes as rows
loc_names = unique(ERPs.grouped_annot.loc(Data_Tag));
loc_names(strcmpi(loc_names, '')) = [];
loc_names(strcmpi(loc_names, 'lax')) = [];
loc_names(strcmpi(loc_names, 'changing')) = [];

count = 1;
loc_combs = [];
for i = 1:length(loc_names)
    for j = (i+1):length(loc_names)
        loc_combs{count} = [loc_names{i}, ' - ', loc_names{j}];
        count = count+1;
    end
end
data_table.handshapes = loc_combs';


%% Find the Number of comparisons between each combination of HS for each location
for i = 1:length(handshapes);
    hs_name = handshapes{i};
    is_hs = Data_Tag & strcmpi(ERPs.grouped_annot.handshape, handshapes(i));
    count = 1;
    for k = 1:length(loc_names)
        for j = (k+1):length(loc_names)
            % Create is_Ongoing and is_onset booleans for hs(k) and hs(j)
            loc_1 = loc_names(k);
            is_loc1_onset = strcmpi(ERPs.grouped_annot.loc, loc_1);
            is_loc1_ongoing = ~is_loc1_onset & strcmpi(ERPs.filled_change_annot.loc, loc_1);
            
            loc_2 = loc_names(j);
            is_loc2_onset = strcmpi(ERPs.grouped_annot.loc, loc_2);
            is_loc2_ongoing = ~is_loc2_onset & strcmpi(ERPs.filled_change_annot.loc, loc_2);
                        
            Loc12_Comp_on_Loc = Equalize_Tag_Sizes(is_hs & is_loc1_onset, is_hs & is_loc2_onset) |  Equalize_Tag_Sizes(is_hs & is_loc1_ongoing, is_hs & is_loc2_ongoing);
            loc_comb_count(count) = sum(Loc12_Comp_on_Loc(:,1));
                        count = count+1;
        end
    end
    if size(loc_comb_count,1) < size(loc_comb_count,2)
        loc_comb_count = loc_comb_count';
    end
    eval(['data_table.' hs_name,' = loc_comb_count']);

end

hs_eq_loc_compare_Table = struct2table(data_table);
