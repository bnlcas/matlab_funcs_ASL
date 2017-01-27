function [output] = Segment_Handshapes_Slide(ERPs, sig_chans, occurances, components)
%confusion_mat
output_confusion_mat = true;    % output is confusion_mat when true, and categories otherwise
use_CV = false;
output = [];
%% This function attempts to classify the handshapes of various handshapes of respones.
% the sig chans is currently a boolean 256x1 array that is TRUE for
% channels that are relevant to handshape segmentation.

% Occurances represents the minimum number of instances that a classifier occurs in the data
% for it to be analyzed 

% time win is the range of time points where the classifier does it's thing

% components specifies the number of PCA components to use. If empty this
% will simply use the full ecog data


% restrict signs to responses
Data_Tag = strcmpi(ERPs.annot.respType,'dup');
ecog = ERPs.ecog(:,:,Data_Tag);


% Get Annotation of signs
sign_data = ERPs.annot.handshape(Data_Tag);
sign_data = strrep(sign_data,' ','');               % Removes blank spaces
categories = unique(sign_data);
categories = categories(~strcmp(categories,''));    % Clear blanks from categories
sparse_category = false(size(categories));          % Will remove categories with too few instances
category_size_thresh = occurances;                           % Min of Three instances per category
for i = 1:length(categories)
    cat_count = sum(strcmpi(sign_data, categories(i)));
    if cat_count < category_size_thresh
        sparse_category(i) = true;
    end
end
categories = categories(~sparse_category);
%categories = {'A'; 'S'};
if ~output_confusion_mat
    output = categories;
else


%% Assemble Y data for Classifiers
% Create a Sign_datapts x Num_Categories Matrix whose colums are data
% vectors containing info on whether a data point is in a category or not
category_data = false(length(sign_data), length(categories));
for i = 1:length(categories)
    category_data(:,i) = strcmp(sign_data, categories(i));   % Boolean Vector
end
%% Randomly Permute Category Data to test predictie model
Randomize_Data = false;
if Randomize_Data
    new_order = randperm(length(sign_data));
    category_data = category_data(new_order,:);
end

%% Assemble X data to to PCA for Classifier
% relevant channels is function call,
% taken here as sig_chans = (sum((stat_grid_dups_p(:,191:211,52)<0.01),2) > 12);
% or Chans with >2/3 significant linear regression of ECoG to to handshape
relevant_chans = sig_chans;
%TRIM BAD CHANNELS:
is_good = is_good_chan(ERPs);
relevant_chans = sig_chans & is_good;

ecog_test = ecog(relevant_chans,:,:);
dims = size(ecog_test);

if isempty(components) | ((components == length(categories))|(components == 0)) 
    ecog_data_pca = ecog_test;
else    
    ecog_shift = reshape(ecog_test,dims(1),(dims(2).*dims(3)));    % crunch timepts and trails together for PCA
    ecog_shift = ecog_shift';    % flips data s.t. channels are cols
    pca_mat = pca(ecog_shift);


    %components = 20;        % select number of components. Must be >1 & <chans
    pca_data = ecog_shift*pca_mat(:,1:components);
    pca_data = pca_data';
    ecog_data_pca = reshape(pca_data, components, dims(2), dims(3)); % switch back into traditional form
end

%% Loop through and make a series of Confusion Matricis:

Analysis_frame = 20;        % Average over 200ms window
time_advance = 5;           % Advance the frame by 50ms each time
total_time_pts = size(ecog,2);
frames = floor((total_time_pts - Analysis_frame)/time_advance); % clippout the last frame
%confusion_mat = zeros(length(category_data),length(category_data),frames);



%% Remove Data From Excluded HandShapes:
groupings = zeros(size(sign_data));
for i = 1:length(categories)
    groupings = groupings + i.*category_data(:,i); % the ith classifier is listed as an i in a vector
end
exclude = (groupings == 0);
groupings = groupings(~exclude);
ecog_data_pca = ecog_data_pca(:,:,~exclude);
category_data = category_data(~exclude,:);

%% Assemble X data for each handshape:

% select relevant fingerspellings to partition and take the 
% mean High Gamma:

% relevant_chans = 1:256;
% relevant_timepts = 101:301;


for k = 1:frames
time_start = (1+time_advance*(k-1));
time_end = time_start + Analysis_frame;
timewin = time_start:time_end;
k

%relevant_timepts = 191:211;
 relevant_timepts= timewin;

ecog_mean = squeeze(mean(ecog_data_pca(:,relevant_timepts,:),2));
% ecog_mean = squeeze(mean(ecog(relevant_chans,relevant_timepts,:),2));

X_data = ecog_mean';

% Convert Data into Priciple Component Channels:
% pca_mat = pca(X_data);
% components = 15;
% X_data = X_data*pca_mat(:,1:(1+components));





%% Create SVM model for is/isnot category for each of the categories
SVMModels = [];     % Initializing this would be prudent...
for i = 1:length(categories)
    Y_data = category_data(:,i);
    % Fit Support Vector Machine:
    
    %SVMModel = fitcdiscr(X_data,Y_data, 'Prior', 'empirical', 'Gamma', 1);
    
    SVMModel = fitcsvm(X_data,Y_data,'KernelFunction','linear', 'ClassNames',logical([0,1]), 'BoxConstraint', 0.1, 'Leaveout', 'on');
    %SVMModel = fitcsvm(X_data,Y_data,'KernelFunction','rbf', 'ClassNames',logical([0,1]), 'BoxConstraint',0.001);
    

   SVMModels{i} = SVMModel;
end


%% GET CONFUSION MATRIX

% assemble groupings
% groupings = zeros(size(category_data,1));
% for i = 1:length(categories)
%     groupings = groupings + i.*category_data(:,i); % the ith classifier is listed as an i in a vector
% end

% assemble predictions:
predictions = zeros(size(sign_data));

%ecog_mean = squeeze(mean(ecog(relevant_chans,timewin,:),2));
%X_data = ecog_mean';
%X_data = PCA_data;

%for i = 1:length(sign_data)
%for i = 1:size(X_data,1);


% groupings_test = groupings(test_indecies);
    %data_point = X_data(i,:);
    scores = zeros(size(X_data,1), length(categories));
    for j = 1:length(categories)


           %[~,score] = predict(SVMModels{j}.Trained{i},data_point);
           %[~,score] = predict(SVMModels{j},X_data_test);
           [~,score] = kfoldPredict(SVMModels{j});
           scores(:,j) = score(:,2);    % score of TRUE for each data point for category j
    end
    [~,ind]=max(scores,[],2);
    predictions = ind;


% Clear uncategorized inputs (handshapes with a singular occurace, etc
%  excluded = (groupings == 0);
%  groupings = groupings(~excluded);
%  predictions = predictions(~excluded);

confusion_mat_temp = confusionmat(groupings, predictions);
confusion_mat_temp = gdivide(confusion_mat_temp,sum(confusion_mat_temp,2)); % normalize Confusion Matrix
confusion_mat(:,:,k) = confusion_mat_temp;

trace(confusion_mat_temp)/sum(confusion_mat_temp(:));
% 
% 
% %% Plot Confusion Matrix:
% % figure;
% 
% %BORROWED FROM STACK EXCHANGE:
% imagesc(confusion_mat);            %# Create a colored plot of the matrix values
% %colormap(gray);  %# Change the colormap to gray (so higher values are
% colorbar;                       
% %   textStrings = num2str(confusion_mat(:),'%0.1f');  %# Create strings from the matrix values to 2 decimals
% %   textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
% %   [x,y] = meshgrid(1:length(categories));   %# Create x and y coordinates for the strings
% %  hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');   
% %  midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
% %  textColors = repmat(confusion_mat(:) < midValue,1,3);  %# Choose white or black for the
%                                               %#   text color of the strings so
%                                               %#   they can be easily seen over
%                                               %#   the background color
%  %set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
%  
%  %***End StackExchange Borrow***
%
%  ax = gca;
%  ax.XTick = 1:length(categories);
%  ax.YTick = 1:length(categories);
%  ax.XTickLabel = categories';
%  ax.YTickLabel = categories';
%  title({'Confusion Matrix of SVM Classifier with Gaussian Kernel for Handshapes wit 5+ Occurances';'With Leaveout Cross Validation Using Mean ECoG HG on [-0.1s,0.1s] on Channels with p<0.01 on Linear Regression'})
% %  
% %  drawnow;
% %  F(k) = getframe;
% %  
% %  class_accuracy(k) = trace(confusion_mat)/(sum(confusion_mat(:)));
% % end
% %  figure; 
% %  plot(ERPs.time_axis(floor(midpt)), class_accuracy)
% %  title({'Classification Accuracy of Gaussian Over Time';'(No Sliding)'})
% %  xlabel('ERPs Time Axis')
% %  
% %  figure;
% %  movie(F)
%  
%  end
%  
% 
% 

end 

output = confusion_mat;
end

end
    
    