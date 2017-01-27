function [categories, cat_count] = Get_Categories(ERPs, Data_Tag, cat_thresh, loc_hs)
%% This function Takes a data Data and a parameter ('loc' or 'handshape' in this incarnation)
% and generates a list of distinct categories in that parameter that are
% marked by Data_Tag more than cat_thresh number of times


%% Get Annotation of Handshapes involved in duplicates
use_loc = false;    use_handshape = false;    use_intmov = false;
if substrcmp(loc_hs,'loc')
    use_loc = true;
end
if substrcmp(loc_hs,'handshape')
    use_handshape = true;
end
if substrcmp(loc_hs,'intMov')
    use_intmov = true;
end

if use_loc
    sign_data = ERPs.grouped_fill_annot.loc(Data_Tag); % Segments on the Basis of
end
if use_handshape
    sign_data = ERPs.grouped_fill_annot.handshape(Data_Tag); % segment with location
end
if use_intmov
    sign_data = ERPs.grouped_fill_annot.intMov(Data_Tag);
end


sign_data = strrep(sign_data,' ','');               % Removes blank spaces
categories = unique(sign_data);
categories = categories(~strcmp(categories,''));    % Clear blanks from categories
sparse_category = false(size(categories));          % Will remove categories with too few instances
category_size_thresh = cat_thresh;                           % Min of Three instances per category
cat_count = zeros(1,length(categories));
for i = 1:length(categories)
    cat_count(i) = sum(strcmpi(sign_data, categories(i)));
end
sparse_category = (cat_count<category_size_thresh);
cat_count = cat_count(~sparse_category);
categories = categories(~sparse_category);

%%Eliminate Changing and Lax Categoreis
cat_count = cat_count(~(strcmpi(categories, 'changing') | strcmpi(categories, 'lax')));

categories = categories(~strcmpi(categories, 'changing'));
categories = categories(~strcmpi(categories, 'lax'));

cat_count = cat_count(~strcmpi(categories,'fs'));
categories = categories(~strcmpi(categories,'fs'));
%cat_count(find(strcmpi(categories, 'neutral'))) = [];
%categories(find(strcmpi(categories, 'neutral'))) = [];


end