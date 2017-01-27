function [] = KmeansOnAll(ERPs, sig_chans)
Data_Tag = strcmpi(ERPs.annot.filledLexTrans,'lexical');
Data_Tag = Data_Tag & is_good_trial(ERPs);
is_hs = ~strcmpi(ERPs.annot.handshape,'') & Data_Tag;
Significant_Electrodes = sig_chans;
%Significant_Electrodes = [201 202 217 218];

%% Get List of Handshapes for Analysis
Handshapes = ERPs.annot.handshape(is_hs);
categories = unique(Handshapes);            % List of HandShapes
cat_count = get_category_size(Handshapes);  % Num Instances of Each HShape
cat_thresh = 10; % Limit on Num of Instances
categories(cat_count < cat_thresh) = [];
cat_count(cat_count < cat_thresh) = [];
reject_cats = strcmpi(categories,'changing') | strcmpi(categories,'lax');
categories = categories(~reject_cats);
cat_count = cat_count(~reject_cats);

%% Restrict Is_HS
Handshapes = ERPs.annot.handshape;
included = false(size(Handshapes));
for i = 1:length(categories);
    included = included | strcmpi(Handshapes, categories(i));
end
a = 1;
is_hs = is_hs & included;   % The above does not care is something is a proper hs or not
Handshapes = Handshapes(is_hs);

Relevant_ecog = squeeze(mean(ERPs.ecog(Significant_Electrodes,191:211,is_hs),2));

data = Relevant_ecog';

cluster = kmeans(data,2);

mds = cmdscale(pdist(data));
figure; plot(mds((cluster==1),1),mds((cluster==1),2),'r.')
hold
plot(mds((cluster==2),1),mds((cluster==2),2),'b.')
xlabel('Dimension 1')
ylabel('Dimension 2')
legend('Cluster 1', 'Cluster 2')
title({'Plot of K-means Clustering of Handshape ERPs on First Two Dimensions of'; 'Multi-Dimensional Scaling of Electrodes'})

a = 1;




hs_c1 = zeros(length(categories),1);        % The number of Instances of Each HShape in Cluster 1
hs_c2 = zeros(length(categories),1);

for i = 1:length(categories)
    hs_inds = strcmpi(Handshapes, categories(i));
    hs_c1(i) = sum(cluster(hs_inds)==1);
    hs_c2(i) = sum(cluster(hs_inds)==2);

end



figure; barh(hs_c1./cat_count)
xlabel('Fraction of Handshape ERPs in Cluster 1')
ax = gca;
ax.YTick = 1:length(categories)';
ax.YTickLabel = categories;

figure; barh(hs_c2./cat_count)
xlabel('Fraction of Handshape ERPs in Cluster 2')
ax = gca;
ax.YTick = 1:length(categories)';
ax.YTickLabel = categories;

a =1;


