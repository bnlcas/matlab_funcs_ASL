function [grouped_annot] = group_annot_alt(annot)
%% This function takes the annot class and returns a re grouped annot - a new data structure
% Grouped Annot is Constructed with the same subcategories, however the
% values have been re odered in the folloing ways;

% Changing Handshapes and locations have been eliminated

% Events with Empty Fillings that contain
% times at the index of the event in the grouping with the longest duration

% The value of each colum in the group is assigned to the member event's
% value with highest duration and a non-empty value in that colum

% The Resulting Structure also has the additional sub_class Grouping,
% denoted as with a 'G' at the rows with filled in data.

%% Convert annotations into a Table
data_table = struct2table(annot);
data_table_grouped = data_table;


%% Group Non-Fingerspelling Lexical Responses:
LexTag = strrep(table2array(data_table(:,11)),' ',''); % Lexical Row with Spaces removed
Lexical_Index = strcmpi(LexTag, 'lexical'); % Returns true on rows with nonfilled Lex tag

FSrespTag = strrep(table2array(data_table(:,15)), ' ','');
Non_FS_Index = ~strcmpi(FSrespTag, 'fs');

Is_Grouping = Lexical_Index & Non_FS_Index;
Grouping_Index = find(Is_Grouping);

%% Create a Data Colum Group tag that is 'G' on Grouping Indecies
group_tag = repmat({''},size(data_table,1),1);
group_tag(Grouping_Index) = {'G'};

%% Assemble Associated Parameters
IsFilled_lex = strcmpi(strrep(table2array(data_table(:,12)),' ',''), 'lexical');
Duration = table2array(data_table(:,7));

%% Loop through Columns of DataMatrix and Cluster things for each grouping
for i = 1:size(data_table,2)
    col = table2array(data_table(:,i));
    if ~isnumeric(col(1))  % only deal with string columns
        % Preprocess
        col = strrep(col,' ','');            %clear spaces
        col(strcmpi(col,'changing')) = {''};   % clear changing tags
        
        for k = 1:sum(Is_Grouping)
            % if the column value the grouping index is empty fill it in
            if strcmp(col(Grouping_Index(k)),'')
                j = Grouping_Index(k);
                longest_entry = 0;
                while IsFilled_lex(j)
                    if ~strcmp(col(j),'')
                       if Duration(j)>longest_entry
                            col(Grouping_Index(k)) = col(j);
                       end
                       %data_table_grouped(Grouping_Index(k),i) = {col(j)};
                        %col(Grouping_Index(k)) = col(j); % Replaces with last non zero entry
                    end
                    j = j+1;
                end
            end
        end  
        data_table_grouped(:,i) = col;
    end
end
 
grouped_annot = table2struct(data_table_grouped, 'ToScalar', true);
grouped_annot.group_tag = group_tag;
         





%% Partition With Time Overlaps ~ 324 Time Clumps

%% Find Start and EndTime Regions
% start_times = table2array(data_table(:,4));                 % Column 3 is the row of starttimes: Magic Constant
% end_times = table2array(data_table(:,6));
% 
% group_start_times = [];
% group_end_times = [];
% k = 1;
% new_group = true;
% for i = 1:length(start_times)
%     if i~=1
%         if start_times(i)<start_times(i-1)
%             new_group = true;
%         end
%     end
%     if new_group
%         group_start_times(k) = start_times(i);
%         group_end_time(k) = end_times(i);
%         new_group = false;
%     else
%         if(start_times(i) > group_end_time(k))
%             new_group = true;
%             k = k+1;
%         else
%             if(end_times(i)> group_end_time(k))
%                 group_end_time(k) = end_times(i);
%             end
%         end
%     end
%     
% end
        


%%

% 
% %% Get Stat and End pts of each event
% start_times = table2array(data_table(:,4));                 % Column 3 is the row of starttimes: Magic Constant
% end_times = table2array(data_table(:,6));
% 
% %% Find size of Table and create Group Tag
% table_length = length(start_times);
% group = repmat({''},table_length,1);
% 
% start_times_cat = unique(start_times);    % Since the Data_Table is ordered by start time, first instance is a helpful index.
% 
% for i = 1:length(start_times_cat)
%     relevant_rows = find((start_times == start_times_cat(i)));
%     [~,longest_row] = max(end_times(relevant_rows));
%     write_index = relevant_rows(longest_row); % write to relevant_rows(longest_row) in Data_Table
%     group(write_index) = {'G'};
%     
%     %% Loop through each column and find longest nonempty tag to fill into grouping
%     for k = 1:size(data_table,2)
%         col = table2array(data_table(relevant_rows,k)); %
%         durations = end_times(relevant_rows);
%         if isnumeric(col(1))
%             is_written = (col ~= 0); % select non-empty values
%             col = col(is_written);
%             if ~isempty(col)
%                 durations = durations(is_written);
%                 [~,data] = max(durations);
%             
%                 data_table_grouped(write_index,k) = {col(data)};
%             end
%         else
%             is_written = ~(strcmpi(col,'') | strcmpi(col,' '));
%             col = col(is_written);
%             if ~isempty(col)
%                 durations = durations(is_written);
%                 [~,data] = max(durations);
%                 data_table_grouped(write_index,k) = {col(data)};
%             end
%         end    
%     end
% end
% grouped_annot = table2struct(data_table_grouped, 'ToScalar', true);
% grouped_annot.group_index = group;
% end