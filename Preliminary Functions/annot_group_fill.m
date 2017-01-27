function [grouped_filled_annot] = annot_group_fill(annot)
%% This function takes the annot class and returns a grouped annot structure
% The new data structure, Grouped Annot, is Constructed with the same subcategories, however the
% values for location and handshape have been duplicated for erps with the
% same onset time.
% Following this, Taggings are filled in until a non-empty erp results.
% ERPs with a unique tagging, or the same time as another ERP with unique
% tagging are labeled Onset, while those whose tag is carried over are
% labeled sustained.

% The Resulting Structure also has the additional sub_class Grouping,
% denoted as with a 'G' at the rows with filled in data for the locations.
% data_table = struct2table(annot);
% data_table_grouped = data_table;

%% Get Stat and End pts of each event
start_times = annot.start_ms*1000;
%start_times = table2array(data_table(:,4))*1000;    % Column 4 is the row of starttimes in seconds: Magic Constant

%% Find size of Table and create Group Tag
table_length = length(start_times);
hs_group = repmat({''},table_length,1);
loc_group = repmat({''},table_length,1);

%% For all filled handshapes test if there is a location to group
handshapes = strrep(annot.handshape,' ','');
locations = strrep(annot.loc,' ','');
intmovs = strrep(annot.intMov,' ','');
contact = strrep(annot.contact,' ','');

%% Label already filled in Tags as onsets
filled_hs = find(~strcmpi(handshapes,''));
hs_group(filled_hs) = {'onset'};
filled_loc = find(~strcmpi(locations,''));
loc_group(filled_loc) = {'onset'};

%% Combine tags if start time is the same
%% Label the first of a set of common ERPs as a unique ERP
event_start_times = unique(start_times);
is_unique_erp = false(length(start_times),1);
for i = 1:length(event_start_times)
    event_inds = find(start_times == event_start_times(i));
    
    is_unique_erp(event_inds(1)) = true;
    if length(event_inds > 1)
        if ~isempty(find(~strcmpi(handshapes(event_inds),''),1))    % Boolean only fills in when there is at least 1 annotation on a set of redundant erps
            handshapes(event_inds) = handshapes(event_inds(find(~strcmpi(handshapes(event_inds),''),1)));
        end  
        if ~isempty(find(~strcmpi(locations(event_inds),''),1)) 
            locations(event_inds) = locations(event_inds(find(~strcmpi(locations(event_inds),''),1)));
        end
        if ~isempty(find(strcmpi(hs_group(event_inds),'onset')))
            hs_group(event_inds) = {'onset'};
        end
        if ~isempty(find(strcmpi(loc_group(event_inds),'onset')))
            loc_group(event_inds) = {'onset'};
        end
        
        if ~isempty(find(~strcmpi(intmovs(event_inds),''),1))
            intmovs(event_inds) = intmovs(event_inds(find(~strcmpi(intmovs(event_inds),''),1)));
        end
        if ~isempty(find(~strcmpi(contact(event_inds),''),1))
            contact(event_inds) = contact(event_inds(find(~strcmpi(contact(event_inds),''),1)));
        end
    end
end
%% Make Tag for Unique Events
% is_unique_erp = false(length(start_times),1);
% for i = 1:length(event_start_times)
%     event_ind = find(start_times == event_start_times(i),1);
%     if ~isempty(event_ind)
%         is_unique_erp(event_ind) = true;
%     end
% end

%% Fill in Previous Classes



for i = 2:length(handshapes)
    if strcmp(handshapes(i),'') & ~strcmpi(handshapes(i-1),'')
        handshapes(i) = handshapes(i-1);
        hs_group(i) = {'sustained'};
    end
end

for i = 2:length(locations)
    if strcmp(locations(i),'') & ~strcmpi(locations(i-1),'')
        locations(i) = locations(i-1);
        loc_group(i) = {'sustained'};
    end
 
end
grouped_filled_annot = annot;
grouped_filled_annot.handshape = handshapes;
grouped_filled_annot.loc = locations;
grouped_filled_annot.grouping_hs = hs_group;
grouped_filled_annot.grouping_loc = loc_group;
grouped_filled_annot.is_unique_erp = is_unique_erp;

grouped_filled_annot.intMov = intmovs;
grouped_filled_annot.contact = contact;


end