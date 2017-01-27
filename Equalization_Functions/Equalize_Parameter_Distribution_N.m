function [eq_inds] = Equalize_Parameter_Distribution_N(ERPs, Data_Tag, equalize_by, varargin)
%% This function takes the ECoG table, a restricting boolean Data_Tag
% and N booleans ranging over the size of the ERP annotation table, and 
% returns a matrix of (length(annotations)) rows by N comparators, whose entries are such
% That each column contains the indecies of ERPs such that each row has an
% equal number of each value of the equalization_by parameter

num_labels = length(varargin);
annot = ERPs.grouped_fill_annot;

%% set equalizer parameter (handshapes or locations or...)
if strcmpi(equalize_by,'handshape')
    eq_param = annot.handshape; % Arbitary; could be loc or intmov or whatever
end
if strcmpi(equalize_by, 'loc')
    eq_param = annot.loc;
end
eq_param = strrep(eq_param, ' ','');

%% restrict parameter and labels to Data_Tag erps by setting other values to false
for i = 1:num_labels
    label = varargin{i};
    label(~Data_Tag) = false;
    varargin{i} = label;
end

% Get Listing of Handshapes in n Classes
categories = unique(eq_param); % full list of possible categories;
exclude_tags = strcmpi(categories,'') | strcmpi(categories,'lax') | strcmpi(categories,'changing');
categories(exclude_tags) = []; % Clear Null Class

for i = 1:num_labels
    label= varargin{i};
    eq_param_i = eq_param(label);
    categories_i = unique(eq_param_i);
    categories = intersect(categories, categories_i);

end
%% Get largest number of each handshape that exists in both classes
cat_count_sub = zeros(num_labels,1); % The size of A category for each label
cat_count = zeros(size(categories)); % The minimum size of EACH category across ALL Labels
for i = 1:length(categories);
    for j = 1:num_labels
        label= varargin{j};
        eq_param_i = eq_param(label);  % Equalization parameter on the jth data tag
        cat_count_sub(j) = sum(strcmpi(eq_param_i,categories(i)));
    end
    cat_count(i) = min(cat_count_sub);
end

%% Define Output Matrix with Equal Indecies
eq_inds = zeros(sum(cat_count), num_labels); % 

for j = 1:num_labels
    label= varargin{j};
    inds_sub = [];
    for i = 1:length(categories)
        ind_cat = find(label & strcmpi(eq_param, categories(i)));
        %% Shuffle and Take cat_count(i)
        shuffle = randperm(length(ind_cat));
        ind_cat = ind_cat(shuffle(1:cat_count(i)));
        inds_sub = [inds_sub; ind_cat];   
    end
    eq_inds(:,j) = inds_sub;
end



