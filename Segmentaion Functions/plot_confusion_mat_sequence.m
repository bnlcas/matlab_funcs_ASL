function [] = plot_confusion_mat_sequence(svm_time_axis, categories, confusion_mat) % confusion_mat_1) %, confusion_mat_2)
%% Function Decomposes A categories x categories x timeframes series of confusion matricies
% and generates a movie and timing information.
% The graph can take a variable number of confusion matricies and plot
% their time series information along side eachother.


cmax = max(confusion_mat(:));
num_mats = size(confusion_mat,3);
cmap = cbrewer('seq','Reds',101);


%% Get List of Accuracy for all Times:
class_accuracy = [];
for j = 1:num_mats
    confusion_mat_local = squeeze(confusion_mat(:,:,j));
    normalized_accuracy = false;
    if normalized_accuracy
        confusion_mat_local = gdivide(confusion_mat_local, sum(confusion_mat_local,2));
    end
    class_accuracy(j) = trace(confusion_mat_local)/(sum(confusion_mat_local(:)));
end
    
 
 %% Draw Beginning Middle End Confusion Matricies:
 [~,ind]=max(class_accuracy);   % Set Max index to be all peak accuracy
 
%  use_median = false;
%  if use_median
%     peak = (ind-3)+ find((median(class_accuracy((ind-2):(ind+2))) == class_accuracy((ind-2):(ind+2))),1);
%     [~,ind] = min(abs(timept+1000));
%     start = (ind-3) + find((median(class_accuracy((ind-2):(ind+2))) == class_accuracy((ind-2):(ind+2))),1);    
%     [~,ind] = min(abs(timept-1000));
%     final = (ind-3) + find((median(class_accuracy((ind-2):(ind+2))) == class_accuracy((ind-2):(ind+2))),1);
%  end
     [~,peak]= max(class_accuracy);
     [~,start] = min(abs(svm_time_axis+1000));
     [~,final] = min(abs(svm_time_axis-1000));

 indecies = [start peak final];

 mean_confusion_mats = true;
 if mean_confusion_mats
     meaning_window = 4; % points on each side of index to mean over
     for i = 1:length(indecies)
         window = (indecies(i)-meaning_window):(indecies(i)+meaning_window);
         relevant_mats(:,:,i) = mean(confusion_mat(:,:,window),3);
     end
 end
     
 
 %Normalize:
 Normalize_confusion = false;
 if Normalize_confusion
     for i = 1:length(indecies)
        confusion_mat_frame = squeeze(relevant_mats(:,:,i));
        confusion_mat_frame = gdivide(confusion_mat_frame,sum(confusion_mat_frame,2));
        normalized_mats(:,:,i) = confusion_mat_frame;
     end
     relevant_mats = normalized_mats;
 else
     relevant_mats = confusion_mat(:,:,indecies);
 end
 
 maxval = max(relevant_mats(:));
 clims = [0, maxval];
 %% Plot Confusion Matrix
 figure;
 for i = 1:length(indecies)
    subplot(1,length(indecies),i);
    confusion_mat_frame = squeeze(relevant_mats(:,:,i));
    accuracy = trace(confusion_mat_frame)/sum(confusion_mat_frame(:));
    imagesc(confusion_mat_frame)
    % colormap(gray);
    colormap(cmap);
    
    ax = gca;
    ax.XTick = 1:length(categories);
    ax.YTick = 1:length(categories);
    ax.XTickLabel = categories';
    ax.YTickLabel = categories';
    caxis([0 maxval]);
    %caxis([0 0.4])
    title({['Normalized Confusion Matrix around t = ' num2str(svm_time_axis(indecies(i))) ' (ms)'];['Accuracy = ' num2str(100*accuracy,3) '%']});
 end

plot_color_bar = true;
if plot_color_bar
    crange = 0:(maxval/100): maxval;
    figure; imagesc(flipud(crange'))
    colormap(cmap);
    colorbar
end