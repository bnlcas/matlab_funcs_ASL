function [order] = ClusterOneClass(ERPs, event_thresh, sig_chans, sample_size, down_sample)
%% Clusters Gestures
% event thresh specifies minimum number of occurances required to classify a data tag
timewin = 191:211; % Specifies relevant Time window as being 200 ms about the event onset
ecog = ERPs.ecog(sig_chans,timewin,:);
ecog = squeeze(mean(ecog,2));

%% USE filled list of annotations or not
%annot = ERPs.grouped_annot;
annot = ERPs.annot;

%% Assemble List of Categories

Labels = annot.intMov;
Labels = strrep(Labels, ' ','');

restrict_to_ling = true;
if restrict_to_ling
    Ling = strcmpi(annot.filledLexTrans, 'lexical');
else
    Ling = true(size(annot.filledLexTrans));
end

Data_Tag = ~strcmpi(Labels, '') & Ling & is_good_trial(ERPs); % Restrict to non-empty lablings of good trials of lingusitic_data
%Data_Tag = Data_Tag & strcmpi(ERPs.annot.respType, 'fs');
%% Reduce Relevant Data to Tagged Data
Labels = Labels(Data_Tag);
ecog = ecog(:,Data_Tag);

%% Get Categories to Cluster
categories = unique(Labels);


%% Threshold Relevant Categories above EVENT_THRESH instances
cat_count = get_category_size(Labels);
sig_cat = (cat_count > event_thresh);
categories = categories(sig_cat);
categories = categories(~strcmpi(categories, 'changing')); % remove Changing


%% Calclate Centroids of the Categories
centroids = zeros(length(categories),size(ecog,1));

%% DownSample
if down_sample
    RandStream.setGlobalStream(RandStream('mcg16807', 'seed', sum(clock)));
    for i = 1:length(categories)

        is_cat = strcmpi(Labels, categories(i));
        relevant_ecog = ecog(:,is_cat);
        shuffle = randperm(sum(is_cat)); % Shuffle
        centroids(i,:) = mean(relevant_ecog(:,shuffle(1:sample_size)),2); % mean ecog in each channel for ERPs in categorie - given category centroid
    end
else
    for i = 1:length(categories)
        is_cat = strcmpi(Labels, categories(i));
        centroids(i,:) = mean(ecog(:,is_cat),2); % mean ecog in each channel for ERPs in categorie - given category centroid

    end
end

%% Cluster centroids Data
distances = pdist(centroids, 'euclidean'); % option to use some other distance metric i.g. mahalanobis

link = linkage(distances, 'ward');
leafOrder = optimalleaforder(link, distances);
figure;
[~,~,order] = dendrogram(link,length(link), 'Reorder',leafOrder, 'Orientation','right', 'ColorThreshold', 0.5*max(link(:,3)));
ax = gca;
ax.YTickLabel = categories(order)';


xlabel('Euclidean Distance')

title(['Clustering of      Using Euclidean Distance On Significant Channels']) % on ' num2str(num_pcs) ' PCs'])
title('Handshape Clustering (Bad Trails Included)')
%% 





% 
% 
% 
% timewin = 191:211; % Specifies relevant Time window as being 200 ms about the event onset
% ecog = ERPs.ecog(sig_chans,timewin,:);
% ecog = squeeze(mean(ecog,2));
% 
% %% USE filled list of annotations or not
% %annot = ERPs.grouped_annot;
% annot = ERPs.annot;
% 
% %% Assemble List of Categories
% 
% Labels = annot.handshape;
% Labels = strrep(Labels, ' ','');
% 
% restrict_to_ling = true;
% if restrict_to_ling
%     Ling = strcmpi(annot.filledLexTrans, 'lexical');
% else
%     Ling = true(size(annot.filledLexTrans));
% end
% 
% Data_Tag = ~strcmpi(Labels, '') & Ling & is_good_trial(ERPs); % Restrict to non-empty labings
% Data_Tag = Data_Tag & strcmpi(ERPs.annot.respType, 'fs');
% %% Reduce Relevant Data to Tagged Data
% Labels = Labels(Data_Tag);
% ecog = ecog(:,Data_Tag);
% 
% %% Get Categories to Cluster
% categories = unique(Labels);
% 
% 
% %% Threshold Relevant Categories above EVENT_THRESH instances
% cat_count = get_category_size(Labels);
% sig_cat = (cat_count > event_thresh);
% categories = categories(sig_cat);
% categories = categories(~strcmpi(categories, 'changing'));
% 
% 
% %% Calclate Centroids of the Categories
% centroids = zeros(length(categories),size(ecog,1));
% 
% %% DownSample
% if down_sample
%     RandStream.setGlobalStream(RandStream('mcg16807', 'seed', sum(clock)));
%     for i = 1:length(categories)
% 
%         is_cat = strcmpi(Labels, categories(i));
%         relevant_ecog = ecog(:,is_cat);
%         shuffle = randperm(sum(is_cat)); % Shuffle
%         centroids(i,:) = mean(relevant_ecog(:,shuffle(1:sample_size)),2); % mean ecog in each channel for ERPs in categorie - given category centroid
%     end
% else
%     for i = 1:length(categories)
%         is_cat = strcmpi(Labels, categories(i));
%         centroids(i,:) = mean(ecog(:,is_cat),2); % mean ecog in each channel for ERPs in categorie - given category centroid
% 
%     end
% end
% 
% %% Cluster centroids Data
% distances = pdist(centroids, 'euclidean'); % option to use some other distance metric i.g. mahalanobis
% 
% link = linkage(distances, 'ward');
% leafOrder = optimalleaforder(link, distances);
% subplot(1,2,2)
% [~,~,order] = dendrogram(link,length(link), 'Reorder',leafOrder, 'Orientation','right', 'ColorThreshold', median(link(:,3)));
% ax = gca;
% ax.YTickLabel = categories(order)';
% xticklabel_rotate()
% 
% xlabel('Euclidean Distance')
% 
% title('Fingerspelling Clustering (Good Trails Only)')