function [distances] = SampleSize_of_Centroids(ERPs, sig_chans, event_thresh, sample_size)



timewin = 191:211; % Specifies relevant Time window as being 200 ms about the event onset
ecog = ERPs.ecog(sig_chans,timewin,:);
ecog = squeeze(mean(ecog,2));

%% USE filled list of annotations or not
%annot = ERPs.grouped_annot;
annot = ERPs.annot;
%annot = ERPs.alt_annot;

%% Assemble List of Categories
% Data_Tag = strcmpi(annot.group_tag, 'G'); % Restrict to Group Signs
Data_Tag = strcmpi(annot.filledLexTrans, 'lexical'); % Restrict to Lexical Signs

Labels = annot.loc;
Labels = strrep(Labels, ' ','');

Data_Tag = Data_Tag & ~ strcmpi(Labels, ''); % Restrict to non-empty labings

%% Reduce Relevant Data to Tagged Data
Labels = Labels(Data_Tag);
ecog = ecog(:,Data_Tag);

%% Get Categories to Cluster
categories = unique(Labels);

%% Threshold Relevant Categories above EVENT_THRESH instances
cat_count = get_category_size(Labels);
sig_cat = (cat_count > event_thresh);
categories = categories(sig_cat);


%% Calclate Centroids of the Categories
centroids = zeros(length(categories),size(ecog,1));

%% DownSample
down_sample = true;
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



%[~,~,order] = dendrogram(link,length(link), 'Reorder',leafOrder, 'Orientation','right'); % 'ColorThreshold', 8);
%ax = gca;
%ax.YTickLabel = categories(order)';
%xticklabel_rotate()
%xlabel('Euclidean Distance')
%title(['Clustering of      Using Euclidean Distance On Significant Channels']) % on ' num2str(num_pcs) ' PCs'])
%% 

