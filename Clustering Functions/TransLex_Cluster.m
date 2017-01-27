function [] = Cluster_Trans_Lex(ERPs, Data_Tag, sig_chans_noLoc)
%% This function Takes the largest handshapes that are common to
% transitional and lexical movements and them clusters them on the basis
% of their centroids in the significant_electectrodes that do not encode
% location

if isempty(Data_Tag)
    Data_Tag = is_good_trial(ERPs);
end
if isempty(sig_chans_noLoc)
    sig_chans_noLoc = true(256,1);
end

handshapes = strrep(ERPs.annot.handshape,' ','');

Threshold_HS_Size = 10;

is_trans = strcmpi(ERPs.annot.filledLexTrans,'transitional');
is_ling = strcmpi(ERPs.annot.filledLexTrans,'lexical');
is_fs = strcmpi(ERPs.grouped_fill_annot.loc,'fingerspelling');
is_lex = is_ling & ~ is_fs;

%% Get indecies containg equal HS distribution in Trans/Lex
comp_inds = Equalize_Parameter_Distribution(ERPs, Data_Tag, is_trans, is_lex);

inds_trans = comp_inds(:,1);
inds_lex = comp_inds(:,2);

hs_trans = handshapes(comp_inds(:,1));
hs_lex = handshapes(comp_inds(:,2));

cats = unique(hs_trans);
cat_sizes = zeros(length(cats),1);
for i = 1:length(cats)
    cat_sizes(i) = sum(strcmpi(hs_trans, cats(i)));
end

%% List of categories to remove from comp_inds
exclude_cat = cats(cat_sizes < Threshold_HS_Size);
cats(cat_sizes < Threshold_HS_Size) = [];

clear_inds = false(size(inds_trans));
for i = 1:length(exclude_cat)
    clear_inds = clear_inds | strcmpi(hs_trans, exclude_cat(i));
end
inds_trans(clear_inds) = [];
clear_inds = false(size(inds_lex));
for i = 1:length(exclude_cat)
    clear_inds = clear_inds | strcmpi(hs_lex, exclude_cat(i));
end
inds_lex(clear_inds) = [];

%% Get Centroids for clustering
timewin = 191:211;
Cluster_Data = ERPs.ecog(sig_chans_noLoc, timewin,:);
Cluster_Data = squeeze(mean(Cluster_Data,2));



Centroids_Lex = zeros(length(cats),size(Cluster_Data,1));
hs_lex = handshapes(inds_lex);
for i = 1:length(cats)
    inds_cat = inds_lex(strcmpi(hs_lex, cats(i)));
    Centroids_Lex(i,:) = mean(Cluster_Data(:,inds_cat),2); % mean ecog in each channel for ERPs in categorie - given category centroid
end

Centroids_Trans = zeros(length(cats),size(Cluster_Data,1));
hs_trans = handshapes(inds_trans);
for i = 1:length(cats)
    inds_cat = inds_trans(strcmpi(hs_trans, cats(i)));
    Centroids_Trans(i,:) = mean(Cluster_Data(:,inds_cat),2); % mean ecog in each channel for ERPs in categorie - given category centroid
end


%Centroids_Lex = mean(Cluster_Data(:,inds_lex),2);
%Centroids_Trans = mean(Cluster_Data(:,inds_trans),2);



%% Cluster Lex
distances = pdist(Centroids_Lex, 'euclidean'); % option to use some other distance metric i.g. mahalanobis

link = linkage(distances, 'ward');
leafOrder = optimalleaforder(link, distances);
figure; [~,~,order] = dendrogram(link,length(link), 'Reorder',leafOrder, 'Orientation','right', 'ColorThreshold', 0.5*max(link(:,3)));
ax = gca;
ax.YTickLabel = cats(order)';
xlabel('Euclidean Distance')
title('Hierarchical Clustering of Lexical Handshapes')



distances = pdist(Centroids_Trans, 'euclidean'); % option to use some other distance metric i.g. mahalanobis

link = linkage(distances, 'ward');
leafOrder = optimalleaforder(link, distances);
figure; [~,~,order] = dendrogram(link,length(link), 'Reorder',leafOrder, 'Orientation','right', 'ColorThreshold', 0.5*max(link(:,3)));
ax = gca;
ax.YTickLabel = cats(order)';
xlabel('Euclidean Distance')
title('Hierarchical Clustering of Transitional Handshapes')
