function [grouped_annot] = annot_group_temp_prox(annot)
%% This function takes the annot class and returns a grouped annot structure
% The new data structure, Grouped Annot, is Constructed with the same subcategories, however the
% values for location have been added into filled handshape entries provided that their onset falls
% within a timespan set by a threshold (here it is set at 50ms)

% The Resulting Structure also has the additional sub_class Grouping,
% denoted as with a 'G' at the rows with filled in data for the locations.
% data_table = struct2table(annot);
% data_table_grouped = data_table;

%% Get Stat and End pts of each event
start_times = annot.start_ms*1000;
%start_times = table2array(data_table(:,4))*1000;    % Column 4 is the row of starttimes in seconds: Magic Constant
time_prox_thresh = 50;  % 50 ms start window
search_window = 8;     % number of entries to search for filled locations

%% Find size of Table and create Group Tag
table_length = length(start_times);
group = repmat({''},table_length,1);

%% For all filled handshapes test if there is a location to group
handshapes = strrep(annot.handshape,' ','');
locations = strrep(annot.loc,' ','');

filled_hs = find(~strcmpi(handshapes,''));
for i = 1:length(filled_hs)
    hs_start = start_times(filled_hs(i));
    range = intersect((filled_hs(i)-floor(search_window/2)):(filled_hs(i)+floor(search_window/2)),1:length(locations));
    for j = range;
        if (abs(start_times(j) - start_times(i)) < time_prox_thresh) & ~strcmpi(locations(j),'')
            locations(i) = locations(j);
        end
    end
    
    % Tag as grouped
    if ~strcmpi(locations(filled_hs(i)),'')
        group(filled_hs(i)) = {'G'};
    end
    
end
grouped_annot = annot;
grouped_annot.loc = locations;
grouped_annot.grouped_tag = group;

end