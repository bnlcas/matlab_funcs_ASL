function [] = ClusterGestures_Comb_grouped(ERPs, Data_Tag, event_thresh, sig_chans, color_thresh)
%% Clusters Gestures
% event thresh specifies minimum number of occurances required to classify a data tag
timewin = 191:211; % Specifies relevant Time window as being 200 ms about the event onset
ecog = ERPs.ecog(sig_chans,timewin,Data_Tag);
ecog = squeeze(mean(ecog,2));

%% USE filled list of annotations or not
annot = ERPs.grouped_fill_annot;
%annot = ERPs.alt_annot;
% Data_Tag = strcmpi(annot.group_tag, 'G');

%% Assemble List of Categories

handshapes = annot.handshape(Data_Tag);
%handshapes = annot.intMov(Data_Tag);
handshapes = strrep(handshapes, ' ',''); % List of Handshapes with (spaces remvoed)
%included_hs = ~(strcmpi(handshapes, '')) & strcmpi(ERPs.annot.filledLexTrans,'lexical'); % List of all meaningful handshapes
exclude_hs = strcmpi(handshapes, 'lax') | strcmpi(handshapes, 'changing') | strcmpi(handshapes, '') | strcmpi(handshapes,'fs');
%locations(exclude_loc) = {''};


locations = annot.loc(Data_Tag);
locations = annot.intMov(Data_Tag);
locations = strrep(locations, ' ','');
exclude_loc = strcmpi(locations, '') | strcmpi(locations,'changing') | strcmpi(locations,'fingerspelling') | strcmpi(locations,'fs');
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
