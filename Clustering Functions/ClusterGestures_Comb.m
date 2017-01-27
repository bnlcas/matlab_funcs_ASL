function [] = ClusterGestures_Comb(ERPs, event_thresh, sig_chans)
%% Clusters Gestures
% event thresh specifies minimum number of occurances required to classify a data tag
timewin = 191:211; % Specifies relevant Time window as being 200 ms about the event onset
ecog = ERPs.ecog(sig_chans,timewin,:);
ecog = squeeze(mean(ecog,2));

%% USE filled list of annotations or not
use_filled_annot = true;
if use_filled_annot
    annot = ERPs.filled_annot;
else
    annot = ERPs.annot;
end

%% Assemble List of Categories

handshapes = annot.handshape;
handshapes = strrep(handshapes, ' ',''); % List of Handshapes with (spaces remvoed)
included_hs = ~(strcmpi(handshapes, '')) & strcmpi(ERPs.annot.filledLexTrans,'lexical'); % List of all meaningful handshapes

locations = annot.loc;
locations = strrep(locations, ' ','');
included_loc = ~(strcmpi(locations, '')) & strcmpi(ERPs.annot.filledLexTrans,'lexical');
%%Only segment Fingerspelling
%included_loc = included_loc & ~strcmpi(locations, 'fingerspelling');
included_loc = strcmpi(locations, 'fingerspelling');

%% List of every HandShape - Location Combination
dashes = repmat('-',length(locations),1);
class_list_raw = strcat(locations, dashes, handshapes);


%% Assemble list of Unique Categories with blanks excluded
included_erps = included_hs & included_loc;
locations = locations(included_erps);
handshapes = handshapes(included_erps);

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
    is_cat = strcmpi(class_list_raw, categories(i));
    centroids(i,:) = mean(ecog(:,is_cat),2); % mean ecog in each channel for ERPs in categorie - given category centroid
end

%% Reduce Dimensionality of Data
% Run PCA on all included ERPs
run_pca = false;
if run_pca
    [reduced,~,~,~,explained] = pca(ecog(:,included_erps)');
    exp_thresh = 90;
    num_pcs = find(cumsum(explained) > exp_thresh,1);
    %num_pcs = 25;
    sum(explained(1:num_pcs))
    centroids = centroids*reduced(:,1:num_pcs);
end

%From Confusion Matrix:

%% Cluster centroids Data
distances = pdist(centroids, 'euclidean'); % option to use some other distance metric i.g. mahalanobis

%distances = pdist(centroids, 'mahalanobis', cov(ecog(:,included_erps)')); % option to use some other distance metric i.g. mahalanobis

link = linkage(distances);
%figure; dendrogram(link,length(link))
leafOrder = optimalleaforder(link, distances);
figure; [~,~,order] = dendrogram(link,length(link), 'Reorder',leafOrder, 'Orientation', 'right'); % 'ColorThreshold', 3.5);
ax = gca;
%order = str2num(ax.XTickLabel);
ax.YTickLabel = categories(order)';
%xticklabel_rotate()
%xlabel({'Sign Category'});
xlabel('Euclidean Distance')

title(['Clustering of Lexical Sign Classifications Using Euclidean Distance On Significant Channels (50+Occurances)']) % on ' num2str(num_pcs) ' PCs'])
%% 

