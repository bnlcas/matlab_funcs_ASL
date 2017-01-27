function cat_count = get_category_size(labels)
%% Takes a list of labels data and returns the number of times that each label appears
categories = unique(labels);
categories = categories(~strcmpi(categories,''));
num_cats = length(categories);
cat_count = zeros(num_cats,1);
for i = 1:num_cats
    cat_count(i) = sum(strcmpi(labels, categories(i)));
end

end
