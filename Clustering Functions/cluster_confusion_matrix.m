function [] = cluster_confusion_matrix(confusion_matrix, categories)
distances = pdist(confusion_matrix, 'cityblock');

link = linkage(distances);
%figure; dendrogram(link,length(link))
leafOrder = optimalleaforder(link, distances);
figure; [~,~,order] = dendrogram(link,length(link), 'Reorder',leafOrder, 'ColorThreshold', 3.5, 'Orientation','right');
ax = gca;
%order = str2num(ax.XTickLabel);
ax.YTickLabel = categories(order)';
%xlabel({'Sign Category'});
xlabel('Distance')

title(['Clustering of Handshape Sign Classifications Using CityBlock Distances on Confusion Matrix']) % on ' num2str(num_pcs) ' PCs'])


end
