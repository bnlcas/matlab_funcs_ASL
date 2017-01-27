function [loc_eq_HS_compare_Table] = make_hs_comp_loc_eq_table(ERPs)
%% Takes class of Data (ERPs.grouped_annot.loc) and creates a table
% Where each column is a type within that Colunm
%% Indexes the ERPs of Handshape Onset with Location Equalized
loc_data = ERPs.grouped_annot.loc;
Data_Tag = is_good_trial(ERPs) & strcmpi(ERPs.annot.filledLexTrans,'lexical');

%% GetList of Locations for column rows:
locations = unique(loc_data(Data_Tag));
locations(strcmpi(locations, '')) = [];
locations(strcmpi(locations,'changing')) = [];

%% Assemble List of Handshapes as rows
hs_names = unique(ERPs.grouped_annot.handshape(Data_Tag));
hs_names(strcmpi(hs_names, '')) = [];
hs_names(strcmpi(hs_names, 'lax')) = [];
hs_names(strcmpi(hs_names, 'changing')) = [];
hs_size = zeros(length(hs_names),1);
% Find biggest_hs
for i = 1:length(hs_names)
    hs_size(i) = sum(strcmpi(ERPs.grouped_annot.handshape(Data_Tag), hs_names(i)));
end
[~, order] = sort(hs_size,'descend');
hs_names = hs_names(order(1:10)); % takes the 10 largest categories

count = 1;
hand_shape_combs = [];
for i = 1:length(hs_names)
    for j = (i+1):length(hs_names)
        hand_shape_combs{count} = [hs_names{i}, ' - ', hs_names{j}];
        count = count+1;
    end
end
data_table.handshapes = hand_shape_combs';


%% Find the Number of comparisons between each combination of HS for each location
for i = 1:length(locations);
    loc_name = locations{i};
    is_loc_onset = Data_Tag & strcmpi(ERPs.grouped_annot.loc, locations(i));
    is_loc_ongoing = Data_Tag & ~ is_loc_onset & strcmpi(ERPs.filled_change_annot.loc, locations(i));
    count = 1;
    for k = 1:length(hs_names)
        for j = (k+1):length(hs_names)
            % Create is_Ongoing and is_onset booleans for hs(k) and hs(j)
            hs_1 = hs_names(k);
            is_hs1_onset = strcmpi(ERPs.grouped_annot.handshape, hs_1);
            %is_hs1_ongoing = ~is_hs1_onset & strcmpi(ERPs.filled_change_annot.handshape, hs_1);
            
            hs_2 = hs_names(j);
            is_hs2_onset = strcmpi(ERPs.grouped_annot.handshape, hs_2);
            %is_hs2_ongoing = ~is_hs2_onset & strcmpi(ERPs.filled_change_annot.handshape, hs_2);
                        
            HS12_Comp_on_Loc = Equalize_Tag_Sizes(is_loc_onset & is_hs1_onset, is_loc_onset & is_hs2_onset) |  Equalize_Tag_Sizes(is_loc_ongoing & is_hs1_onset, is_loc_ongoing & is_hs2_onset);
            hs_comb_count(count) = sum(HS12_Comp_on_Loc(:,1));
                        count = count+1;
                        if strcmpi(loc_name,'hand') & strcmpi(hs_1,'5') & strcmpi(hs_2,'B')
                             PlotECogGrid_Gen(ERPs, true, ERPs.ecog(:,:,HS12_Comp_on_Loc(:,2)), ERPs.ecog(:,:,HS12_Comp_on_Loc(:,1)))
                        end

%                         if strcmpi(loc_name ,'neutral') & strcmpi(hs_1,'O') & strcmpi(hs_2,'5')
%                              PlotECogGrid_Gen(ERPs, true, ERPs.ecog(:,:,HS12_Comp_on_Loc(:,1)), ERPs.ecog(:,:,HS12_Comp_on_Loc(:,2)))
%                         end
%                         if strcmpi(loc_name ,'neutral') & strcmpi(hs_1,'O') & strcmpi(hs_2,'S')
%                              PlotECogGrid_Gen(ERPs, true, ERPs.ecog(:,:,HS12_Comp_on_Loc(:,1)), ERPs.ecog(:,:,HS12_Comp_on_Loc(:,2)))
%                         end
        end
    end
    if size(hs_comb_count,1) < size(hs_comb_count,2)
        hs_comb_count = hs_comb_count';
    end
    eval(['data_table.' loc_name,' = hs_comb_count;']);

end

loc_eq_HS_compare_Table = struct2table(data_table);
