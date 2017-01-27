function [] = ClusterGestures(ERPs, event_thresh)
%% Clusters Gestures
% event thresh specifies minimum number of occurances required to classify a data tag
Data_Tag = strcmpi(ERPs.annot.filledLexTrans, 'lexical');

timewin = 191:211; % Specifies relevant Time window as being 200 ms about the event onset
ecog = ERPs.ecog(:,timewin,Data_Tag);
ecog = squeeze(mean(ecog,2));

%% USE filled list of annotations or not
use_filled_annot = false;
if use_filled_annot
    annot = ERPs.filled_annot;
else
    annot = ERPs.annot;
end


%% Assemble List of Categories

use_handshapes = true;
use_locations = false;
use_movements = false;


include_ERP = false(size(annot.handshape(Data_Tag)));
categories = {};
centroids = [];

if use_handshapes
    handshapes = annot.handshape(Data_Tag);
    handshapes = strrep(handshapes, ' ','');
    include_hs = ~(strcmpi(handshapes, ''));
    categories_hs = unique(handshapes(include_hs));
    cat_count = zeros(size(categories_hs));
    for i = 1:length(categories_hs)
        cat_count(i) = sum(strcmpi(handshapes, categories_hs(i)));
    end
    above_thresh = cat_count > event_thresh;
    
    categories_hs = categories_hs(above_thresh);
    cat_count = cat_count(above_thresh);
    
    % Take mean of relevant Data Points in Class
    cat_centroids = zeros(length(categories_hs),256);
    for i = 1:length(categories_hs)
        is_cat = strcmpi(handshapes, categories_hs(i));
        cat_centroids(i,:) = mean(ecog(:,is_cat),2);
    end
    
    
    categories =  [categories; categories_hs];
    include_ERP = include_ERP & include_hs;
    centroids = [centroids; cat_centroids];
end

if use_locations
    locations = annot.loc(Data_Tag);
    locations = strrep(locations, ' ','');
    include_loc = ~(strcmpi(locations, ''));
    categories_loc = unique(locations(include_loc));
    cat_count = zeros(size(categories_loc));
    for i = 1:length(categories_loc)
        cat_count(i) = sum(strcmpi(locations, categories_loc(i)));
    end
    above_thresh = cat_count > event_thresh;
    
    categories_loc = categories_loc(above_thresh);
    cat_count = cat_count(above_thresh);
    
    %categories_loc(strcmpi(categories_loc,'fingerspelling')) = [];
    % Take mean of relevant Data Points in Class
    cat_centroids = zeros(length(categories_loc),256);
    for i = 1:length(categories_loc)
        is_cat = strcmpi(locations, categories_loc(i));
        cat_centroids(i,:) = mean(ecog(:,is_cat),2);
    end
        
    for i = 1:length(categories_loc);
        if sum(strcmpi(categories, categories_loc(i))) > 0
            categories_loc(i) = strcat(categories_loc(i), '_{loc}');
        end
    end
    
    categories =  [categories; categories_loc];
    include_ERP = include_ERP & include_loc;
    centroids = [centroids; cat_centroids];
end
%% Add categories for  Movements
distances = pdist(centroids, 'euclidean');
link = linkage(distances);

leafOrder = optimalleaforder(link, distances);
figure; [~,~,order] = dendrogram(link,length(link), 'Reorder',leafOrder, 'Orientation','right');
ax = gca;
ax.YTickLabel = categories(order)';
%ylabel({'';'Sign Category'});
xlabel('Euclidean Distance')

title('Clustering of Lexical Handshapes with >75 Occurances Using Euclidean Distance')
%% 

