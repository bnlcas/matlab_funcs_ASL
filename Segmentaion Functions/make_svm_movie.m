function [Film] = make_svm_movie(ERPs, categories, varargin) % confusion_mat_1) %, confusion_mat_2)
%% Function Decomposes A categories x categories x timeframes series of confusion matricies
% and generates a movie and timing information.
% The graph can take a variable number of confusion matricies and plot
% their time series information along side eachother.

Analysis_frame = 10;        % Average over 200ms window
time_advance = 5;           % Advance the frame by 50ms each time
total_time_pts = 401;
% frames = floor((total_time_pts - Analysis_frame)/time_advance) - 1; % clippout the last frame
num_mats = length(varargin);
confusion_mat = varargin{1};
frames = size(confusion_mat,3);
%confusion_mat = zeros(length(category_data),length(category_data),frames);
cmax = max(confusion_mat(:));

fig = figure;
clear F
confusion_mat_frame_all = [];
class_accuracy = [];

for k = 1:frames
    
    for j = 1:num_mats
        confusion_mat = varargin{j};
        confusion_mat = confusion_mat;
        confusion_mat_frame_all(:,:,j) = squeeze(confusion_mat(:,:,k));
    end

%     confusion_mat_frame = squeeze(confusion_mat(:,:,k));
%     confusion_mat_frame_1 = squeeze(confusion_mat_1(:,:,k));
%     %confusion_mat_frame_2 = squeeze(confusion_mat_2(:,:,k));
    
    for j = 1:num_mats
        confusion_mat_local = squeeze(confusion_mat_frame_all(:,:,j));
        normalized_accuracy = true;
        if normalized_accuracy
            confusion_mat_local = gdivide(confusion_mat_local, sum(confusion_mat_local,2));
        end

        %confusion_mat_local = diag(cat_weight)*confusion_mat_local;
        class_accuracy(k,j) = trace(confusion_mat_local)/(sum(confusion_mat_local(:)));
    end
       % class_accuracy(k) = trace(confusion_mat_frame)/(sum(confusion_mat_frame(:)));
%     class_accuracy_1(k) = trace(confusion_mat_frame_1)/(sum(confusion_mat_frame_1(:)));
%     %class_accuracy_2(k) = trace(confusion_mat_frame_2)/(sum(confusion_mat_frame_2(:)));

    
    
    timept(k) = ERPs.time_axis(floor((k-1)*time_advance + Analysis_frame/2));
     confusion_mat_frame = squeeze(confusion_mat(:,:,k));
    confusion_mat_frame = gdivide(confusion_mat_frame, sum(confusion_mat_frame,2));
%      imagesc(confusion_mat_frame);
%      caxis([0, 0.4]);
%      title({'Normalized Confusion Matrix'; ['Time = ' num2str(timept(k)/1000, '%.2f') 's']; ['Accuracy =  ' num2str(class_accuracy(k)*100, '%.2f') '%']}, 'FontSize', 8)   
%      colorbar;
%           
%       ax = gca;
%       ax.XTick = 1:length(categories);
%       ax.YTick = 1:length(categories);
%       ax.XTickLabel = categories';
%       ax.YTickLabel = categories';
% %     
% %     title({'Confusion Matrix of Handshape Classifications for LDA with a Linear Kernel';'Using Average ECoG over 200ms on 29 Channels with Significant PValues';['Time = ' num2str(timept(k))]})   
%      drawnow;
%     F(k) = getframe(fig);
% %    
% %    
end
% figure;
% axis off
% movie(F,1,2)
%Film = F;
% 
%  

plot_accruacy = true;
if plot_accruacy
    figure;
    plot(timept, class_accuracy)%, timept, class_accuracy_1)%, timept, class_accuracy_2)
    title({'Time Evolution of Handshape Classification Accuracy';'On Channels'})
    xlabel('Time From Handshape Onset (ms)')
    ylabel('Classification Accuracy')
    ax = gca;
    ax.YLim(1)=0;
    hold
    plot(ax.XLim, [1/length(categories), 1/length(categories)],'k--')
    
end
% 
 [~,ind]=max(class_accuracy(:,1));
% timept(ind)
class_accuracy(ind);
 
 
 %% Draw Beginning Middle End Confusion Matricies:
 %Find typical points at start middle and end and plot their confusion
 %matricies
 [~,ind]=max(class_accuracy);
 
 use_median = false;
 if use_median
    peak = (ind-3)+ find((median(class_accuracy((ind-2):(ind+2))) == class_accuracy((ind-2):(ind+2))),1);
    [~,ind] = min(abs(timept+1000));
    start = (ind-3) + find((median(class_accuracy((ind-2):(ind+2))) == class_accuracy((ind-2):(ind+2))),1);    
    [~,ind] = min(abs(timept-1000));
    final = (ind-3) + find((median(class_accuracy((ind-2):(ind+2))) == class_accuracy((ind-2):(ind+2))),1);
 else
     [~,peak]= max(class_accuracy(:,1));
     [~,start] = min(abs(timept+1000));
     [~,final] = min(abs(timept-1000));
 end
 indecies = [start peak final];

 mean_confusion_mats = true;
 if mean_confusion_mats
     meaning_window = 1; % points on each side of index to mean over
     for i = 1:length(indecies)
         window = (indecies(i)-meaning_window):(indecies(i)+meaning_window);
         relevant_mats(:,:,i) = mean(confusion_mat(:,:,window),3);
     end
 end
     
 
 %Normalize:
 Normalize_confusion = true;
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
 %% Rearrange ConfusionMatrix
%  distances = pdist2(test_mat, test_mat, 'cityblock');
%  for 
%  [~,min] = min(distances(:,1:4),[],2)
%  
 figure;
 for i = 1:length(indecies)
    subplot(1,3,i);
    confusion_mat_frame = squeeze(relevant_mats(:,:,i));
    accuracy = trace(confusion_mat_frame)/sum(confusion_mat_frame(:));
    imagesc(confusion_mat_frame)
    % colormap(gray);
    
    ax = gca;
    ax.XTick = 1:length(categories);
    ax.YTick = 1:length(categories);
    ax.XTickLabel = categories';
    ax.YTickLabel = categories';
    caxis([0 maxval]);
    %caxis([0 0.4])
    title({['Normalized Confusion Matrix around t = ' num2str(timept(indecies(i))+floor(Analysis_frame*10/2)) ' (ms)'];['Accuracy = ' num2str(100*accuracy,3) '%']});
 end
%  subplot(1,3,3);
%  colorbar('eastoutside')
make_color_bar(clims, 'parula');

%  figure;
%  
%  imagesc(squeeze(confusion_mat(:,:,ind)));
%  title(['T = ' num2str(timept(k)/1000, 3) 's'], 'FontSize', 8)   
%      %colorbar;
%          
%      ax = gca;
%      ax.XTick = 1:length(categories);
%      ax.YTick = 1:length(categories);
%      ax.XTickLabel = categories';
%      ax.YTickLabel = categories';
%     
%     title({'Confusion Matrix of Handshape Classifications at Maximum Accuracy';'Using Average ECoG over 200ms on 29 Channels with Significant PValues';['Time = ' num2str(timept(ind)) ' (ms post onset)   '  '   Classification Accuracy = ' num2str(class_accuracy(ind))]})   
%     drawnow;
% 
