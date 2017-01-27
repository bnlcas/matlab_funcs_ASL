function [] = ClusterGestures_Comb_grouped_variable_comp(ERPs, Data_Tag, event_thresh, sig_chans, color_thresh, group1, group2)
%% Clusters Gestures
% This is equivalent to ClusterGestures_Comb_grouped, function, taking ERPs
% for some Data_Taging, however it takes an input group1, group2,
% group1, group2, must be given as 'loc', 'intMov', or 'handshapes'
% and these must be different due to lazyness on the part of Ben Lucas


% event thresh specifies minimum number of occurances required to classify a data tag
timewin = 191:211; % Specifies relevant Time window as being 200 ms about the event onset
ecog = ERPs.ecog(sig_chans,timewin,Data_Tag);
ecog = squeeze(mean(ecog,2));

%% USE filled list of annotations or not
annot = ERPs.grouped_fill_annot;
%annot = ERPs.alt_annot;
% Data_Tag = strcmpi(annot.group_tag, 'G');

%% Assemble List of Categories
if substrcmp(group1, 'handshape') | substrcmp(group2, 'handshape')
    handshapes = annot.handshape(Data_Tag);
else
    handshapes = annot.intMov(Data_Tag); % if neither group is handshape then one must be intMov, and it takes the place of handshape
end
handshapes = strrep(handshapes, ' ',''); % List of Handshapes with (spaces remvoed)
%included_hs = ~(strcmpi(handshapes, '')) & strcmpi(ERPs.annot.filledLexTrans,'lexical'); % List of all meaningful handshapes
exclude_hs = strcmpi(handshapes, 'lax') | strcmpi(handshapes, 'changing') | strcmpi(handshapes, '') | strcmpi(handshapes,'fs');
%locations(exclude_loc) = {''};

if substrcmp(group1, 'loc') | substrcmp(group2, 'loc')
    locations = annot.loc(Data_Tag);
else
    locations = annot.intMov(Data_Tag);
end
locations = strrep(locations, ' ','');
exclude_loc = strcmpi(locations, '') | strcmpi(locations,'changing') | strcmpi(locations,'fs'); %strcmpi(locations,'fingerspelling');
%locations(exclude_loc) = {''};

%% List of every HandShape - Location Combination
% dashes = repmat('-',length(locations),1);
% class_list_raw = strcat(locations, dashes, handshapes);


%% Assemble list of Unique Categories with blanks excluded
included_erps = ~exclude_hs & ~exclude_loc;
locations = locations(included_erps);
handshapes = handshapes(included_erps);
ecog = ecog(:,included_erps);

dashes = repmat('-',length(locations),1);
class_names = strcat(locations, dashes, handshapes); % add a dash in the middle for legibility
categories = unique(class_names);



%% Threshold Relevant Categories above EVENT_THRESH instances
cat_count = zeros(size(categories));
for i = 1:length(cat_count)
   cat_count(i) = sum(strcmpi(class_names, categories(i))); 
end
sig_cat = (cat_count > event_thresh);
categories = categories(sig_cat);


centroids = zeros(length(categories),size(ecog,1));

for i = 1:length(categories)
    is_cat = strcmpi(class_names, categories(i));
    centroids(i,:) = mean(ecog(:,is_cat),2); % mean ecog in each channel for ERPs in categorie - given category centroid
%    centroids(i,:) = median(ecog(:,is_cat),2);
end

%% Reduce Dimensionality of Data
% Run PCA on all included ERPs
run_pca = false;
if run_pca
    include_inpca = strcmpi(ERPs.annot.filledLexTrans,'lexical');
    %include_inpca = included_erps;
    [reduced,~,~,~,explained] = pca(ecog(:,include_inpca)');
    exp_thresh = 90;
    num_pcs = find(cumsum(explained) > exp_thresh,1)
    %num_pcs = 25;
    sum(explained(1:num_pcs))
    centroids = centroids*reduced(:,1:num_pcs);
end

%From Confusion Matrix:

%% Cluster centroids Data
distances = pdist(centroids, 'euclidean'); % option to use some other distance metric i.g. mahalanobis

%distances = pdist(centroids, 'mahalanobis', cov(ecog(:,included_erps)')); % option to use some other distance metric i.g. mahalanobis

link = linkage(distances, 'ward');
%figure; dendrogram(link,length(link))
leafOrder = optimalleaforder(link, distances);
figure; [~,~,order] = dendrogram(link,length(link), 'Reorder',leafOrder, 'Orientation','right', 'ColorThreshold', color_thresh); % 0.5*max(link(:,3)));
ax = gca;
%order = str2num(ax.XTickLabel);
ax.YTickLabel = categories(order)';
%xticklabel_rotate()
%xlabel({'Sign Category'});
xlabel('Euclidean Distance')

title(['Clustering of Lexical Sign Classifications Using Euclidean Distance On Significant Channels']) % on ' num2str(num_pcs) ' PCs'])
%% 

a=1;

%% Addendum - plot groupings of data:
ecog_full = ERPs.ecog(:,:,Data_Tag);
ecog_full = ecog_full(:,:,included_erps);


is_fs_open = strcmpi(class_names, 'fingerspelling-W')|strcmpi(class_names, 'fingerspelling-C')|strcmpi(class_names, 'fingerspelling-B');
is_fs_closed = strcmpi(class_names, 'fingerspelling-S')|strcmpi(class_names, 'fingerspelling-A')|strcmpi(class_names, 'fingerspelling-T')|strcmpi(class_names, 'fingerspelling-M');
%PlotECogGrid_Gen(ERPs,true, ecog_full(:,:,is_fs_closed),ecog_full(:,:,is_fs_open))

is_face_open = strcmpi(class_names,'face-5')| strcmpi(class_names,'face-Y')|strcmpi(class_names,'face-B');
is_neut_open = strcmpi(class_names, 'neutral-5')|strcmpi(class_names, 'neutral-C')|strcmpi(class_names, 'neutral-B');
%PlotECogGrid_Gen(ERPs,true, ecog_full(:,:,is_neut_open),ecog_full(:,:,is_face_open))

open_hs = {'W', 'C','B', 'F', 'O'};
closed_hs = {'A', 'M','n','t','S','E'};

fsn = 'fingerspelling';

is_fs_open = false(size(ecog_full,3),1);
for i = 1:length(open_hs)
    is_fs_open = is_fs_open | strcmpi(class_names, [fsn, '-', open_hs{i}]);
end
is_fs_closed = false(size(ecog_full,3),1);
for i = 1:length(open_hs)
    is_fs_closed = is_fs_closed | strcmpi(class_names, [fsn, '-', closed_hs{i}]);
end

a=1;
PlotECogGrid_Gen(ERPs,true, ecog_full(:,:,is_fs_closed),ecog_full(:,:,is_fs_open))



%%
PlotECogGrid_Gen(ERPs,true, ecog_full(:,:,is_fs_closed),ecog_full(:,:,is_fs_open), ecog_palm)
% EDIT in Debug
tmp = mean_mat(:,:,1) - mean_mat(:,:,2);
mean_mat(:,:,1) = tmp;
mean_mat(:,:,2) = mean_mat(:,:,3);
mean_mat(:,:,3) = [];
tmp = sqrt(sem_mat(:,:,1).^2 + sem_mat(:,:,2).^2);
sem_mat(:,:,1) = tmp;
sem_mat(:,:,2) = sem_mat(:,:,3);
sem_mat(:,:,3) = [];
num_plots = 2;


