function [eq_inds] = Equalize_Parameter_Distribution(eq_param, Data_Tag, labels1, labels2)
%% This function takes the ECoG table, as well as two data Booleans and Returns 
% Returns of 2 x n vector whose colums represent indecies for 
% parameter equalized ERPs of the label reduced ECoG matrix

%% set equalizer parameter (handshapes or locations or...)
%eq_param = ERPs.annot.handshape; % Arbitary; could be loc or intmov or whatever
eq_param = strrep(eq_param,' ', '');

%% restrict parameter and labels to Data_Tag erps by setting other values to false
%eq_param = eq_param(Data_Tag);
labels1(~Data_Tag) = false;
labels2(~Data_Tag) = false;

% Get Listing of Handshapes in both Classes
eq_param1 = eq_param(labels1);
eq_param2 = eq_param(labels2);

categories1 = unique(eq_param1);
categories2 = unique(eq_param2);

categories = intersect(categories1, categories2);
exclude_tags = strcmpi(categories,'') | strcmpi(categories,'lax') | strcmpi(categories,'changing');
categories(exclude_tags) = []; % Clear Null Class

% Get largest number of each handshape that exists in both classes
cat_count = zeros(size(categories));
for i = 1:length(categories);
    cat_count(i) = min(sum(strcmpi(eq_param1,categories(i))), sum(strcmpi(eq_param2,categories(i))));
end

%% Take the first
ind1 = [];
for i = 1:length(categories)
    ind_cat = find(labels1 & strcmpi(eq_param, categories(i)));
    %% Shuffle and Take cat_count(i)
    shuffle = randperm(length(ind_cat));
    ind_cat = ind_cat(shuffle(1:cat_count(i)));
    ind1 = [ind1; ind_cat];   
end

ind2 = [];
for i = 1:length(categories)
    ind_cat = find(labels2 & strcmpi(eq_param, categories(i)));
    %% Shuffle and Take cat_count(i)
    shuffle = randperm(length(ind_cat));
    ind_cat = ind_cat(shuffle(1:cat_count(i)));
    ind2 = [ind2; ind_cat];   
end

eq_inds= [ind1, ind2];


