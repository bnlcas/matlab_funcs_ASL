function [] = ClusterGestures_Comb_clustered(ERPs, event_thresh, sig_chans)
%% Clusters Gestures
% event thresh specifies minimum number of occurances required to classify a data tag
timewin = 191:211; % Specifies relevant Time window as being 200 ms about the event onset
ecog = ERPs.ecog(sig_chans,timewin,:);
ecog = squeeze(mean(ecog,2));

%% USE grouped list of annotations or not
annot = ERPs.grouped_annot;
Data_Tag = strcmpi(annot.group_tag, 'G') & is_good_trial(ERPs);

%% Assemble List of Categories

handshapes = annot.handshape;
handshapes = strrep(handshapes, ' ',''); % List of Handshapes with (spaces remvoed)
included_hs = ~(strcmpi(handshapes, '')) & strcmpi(ERPs.annot.filledLexTrans,'lexical'); % List of all meaningful handshapes

locations = annot.loc;
locations = strrep(locations, ' ','');
included_loc = ~(strcmpi(locations, '')) & strcmpi(ERPs.annot.filledLexTrans,'lexical');
% 
% movDir = annot.movDir;
% movDir = strrep(movDir, ' ','');
% 
% movPath = annot.movPath;
% movPath = strrep(movPath, ' ','');

intMov = annot.intMov;
intMov = strrep(intMov, ' ','');

%% Include/Exclude Locations
%included_loc = included_loc & ~strcmpi(locations, 'fingerspelling');
%included_loc = strcmpi(locations, 'fingerspelling');

%% List of every HandShape - Location Combination
dashes = repmat('-',length(locations),1);
class_list_raw = strcat(locations, dashes, handshapes, dashes, intMov); %movPath);% , dashes, movDir); % add a dash in the middle for legibility



%% Assemble list of Unique Categories with blanks excluded
included_erps = Data_Tag;
included_erps = ~strcmp(intMov, '') & Data_Tag;

locations = locations(included_erps);
handshapes = handshapes(included_erps);
% movDir = movDir(included_erps);
% movPath = movPath(included_erps);
intMov = intMov(included_erps);



dashes = repmat('-',length(locations),1);
class_names = strcat(locations, dashes, handshapes, dashes, intMov); %movPath); %, dashes, movDir); % add a dash in the middle for legibility
categories = unique(class_names);

%% Threshold Relevant Categories above EVENT_THRESH instances
cat_count = get_category_size(class_names);
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
    include_inpca = strcmpi(ERPs.annot.filledLexTrans,'lexical');
    %include_inpca = included_erps;
    [reduced,~,~,~,explained] = pca(ecog(:,include_inpca)');
    exp_thresh = 90;
    num_pcs = find(cumsum(explained) > exp_thresh,1)
    %num_pcs = 25;
    sum(explained(1:num_pcs))
    centroids = centroids*reduced(:,1:num_pcs);
end

%% Cluster centroids Data
distances = pdist(centroids, 'euclidean'); % option to use some other distance metric i.g. mahalanobis

%distances = pdist(centroids, 'mahalanobis', cov(ecog(:,included_erps)')); % option to use some other distance metric i.g. mahalanobis

link = linkage(distances, 'ward');
%figure; dendrogram(link,length(link))
leafOrder = optimalleaforder(link, distances);
figure; [~,~,order] = dendrogram(link,length(link), 'Reorder',leafOrder, 'Orientation','right', 'ColorThreshold', 0.5*max(link(:,3)));
ax = gca;
%order = str2num(ax.XTickLabel);
ax.YTickLabel = categories(order)';
%xticklabel_rotate()
%xlabel({'Sign Category'});
xlabel('Euclidean Distance')

title(['Clustering of Lexical Sign Classifications Using Euclidean Distance On Significant Channels'])
%% 

