function [] = compare_equal_distribution(ERPs)
%% This function uses the ERP structure to to select ecog from trails containing properties
% Like Lexical vs Transitional, and downsamples to select an equal
% distribution of handshapes for the sake of comparison

%% Get Parameter List
data_type1 = is_good_trial(ERPs) & (strcmpi(ERPs.annot.filledLexTrans, 'lexical') & ~strcmpi(ERPs.annot.respType,'fs')); %NonFingerspelling Lexical Signs
data_type2 = is_good_trial(ERPs) & (strcmpi(ERPs.annot.filledLexTrans, 'lexical') & strcmpi(ERPs.annot.respType,'fs'));  %Lexical Fingerspellings

%Past Comparison:
% data_type1 = strcmpi(ERPs.annot.filledLexTrans, 'transitional');
% data_type2 = strcmpi(ERPs.annot.filledLexTrans, 'lexical');



%% Get handshapes involved in each
handshapes1 = ERPs.annot.handshape(data_type1);
%handshapes1 = ERPs.annot.handshape(data_type1);
handshapes1 = strrep(handshapes1, ' ', '');  % remove emptyspaces
categories1 = unique(handshapes1);
categories1(strcmpi(categories1,'')) = []; % clear empty cells


handshapes2 = ERPs.annot.handshape(data_type2);
handshapes2 = strrep(handshapes2, ' ', ''); %remove emptyspaces from handshape names
categories2 = unique(handshapes2);
categories1(strcmpi(categories1,'')) = [];

categories = intersect(categories1, categories2);
categories(strcmpi(categories,'neutral'))=[];

%% Assemble a count of each handshape used taken to be the minimun number
% of occurances of this category for each 
for i = 1:length(categories);
    cat_count(i) = min(sum(strcmpi(handshapes1,categories(i))), sum(strcmpi(handshapes2,categories(i))));
end

%% Shuffle order of events and take the first cat_count occurances of each
%%handshape into the average

ecog = ERPs.ecog;


%data_type1 = data_type1(neworder); % shuffle items in data_type1

ecog1 = ecog(:,:,data_type1);

shuffle = randperm(length(handshapes1));
shuffle = randperm(length(handshapes1));
ecog1 = ecog1(:,:,shuffle);
handshapes1 = handshapes1(shuffle);

indecies = [];
for i = 1:length(categories)
    indecies = [indecies; find(strcmpi(handshapes1, categories(i)), cat_count(i))];   
end
%indecies = neworder(indecies);

ecog1 = ecog1(:,:,indecies);


%% shuffle for data_type2

% data_type2 = data_type2(neworder); % shuffle items in data_type1


ecog2 = ecog(:,:,data_type2);

shuffle = randperm(length(handshapes2));
shuffle = randperm(length(handshapes2));
ecog2 = ecog2(:,:,shuffle);
handshapes2 = handshapes2(shuffle);

indecies = [];
for i = 1:length(categories)
    indecies = [indecies; find(strcmpi(handshapes2, categories(i)), cat_count(i))];
end
%indecies = indecies(neworder);
ecog2 = ecog2(:,:,indecies);

%% Plot ERPs
PlotECogGrid_Gen(ERPs, false, true, ecog1, ecog2)


mean1_raw = mean(ecog(:,:,data_type1),3);
mean2_raw = mean(ecog(:,:,data_type2),3);
diff_squared = (mean1_raw - mean2_raw).^2;
raw_sum_diff_squared = sum(diff_squared(:))

mean1_eq = mean(ecog1,3);
mean2_eq = mean(ecog2,3);
diff_squared = (mean1_eq - mean2_eq).^2;
eq_sum_diff_squared = sum(diff_squared(:))
% 
% %% RUN SVM ON DATA:
% 
% 
% %% Get Annotation of Handshapes involved in duplicates
% sign_data = ERPs.annot.handshape(Data_Tag);
% sign_data = strrep(sign_data,' ','');               % Removes blank spaces
% categories = unique(sign_data);
% categories = categories(~strcmp(categories,''));    % Clear blanks from categories
% sparse_category = false(size(categories));          % Will remove categories with too few instances
% category_size_thresh = occurances;                           % Min of Three instances per category
% cat_count = zeros(1,length(categories));
% for i = 1:length(categories)
%     cat_count(i) = sum(strcmpi(sign_data, categories(i)));
%     if cat_count(i) < category_size_thresh
%         sparse_category(i) = true;
%     end
% end
% cat_count = cat_count(~sparse_category);
% categories = categories(~sparse_category);
% % categories = {'S'; 'A'}; %- useful to manually select categories
% if ~output_confusion_mat
%     output = categories;    % Return list of categories
%     out_cat_size = false; % optional output of magnetude of each category (may be
%     if out_cat_size
%         output = cat_count;   
%     end
% else
% 
% 
% %% Assemble Y data for Classifiers
% % Create a num_handshapes x Num_Categories Matrix whose colums are data
% % vectors containing info on whether a data point is in a category or not
% category_data = false(length(sign_data), length(categories));
% for i = 1:length(categories)
%     category_data(:,i) = strcmp(sign_data, categories(i));   % Boolean Vector
% end
% 
% groupings = double(category_data)*transpose(1:length(categories)); % numbering for each category (ith category is i...)
% Is_categorized = (sum(category_data,2) ~= 0);   % boolean is false if a Data point is not categorized
% 
% 
% %% Randomly Permute Category Data to test model
% Randomize_Data = false;
% if Randomize_Data
%     new_order = randperm(length(sign_data));
%     category_data = category_data(new_order,:);
%     groupings = groupings(new_order);
%     Is_categorized = Is_categorized(new_order);
% end
% 
% 
% %% Assemble X data to to PCA for Classifier
% % relevant channels is function call,
% % taken here as sig_chans = (sum((stat_grid_dups_p(:,191:211,52)<0.01),2) > 12);
% % or Chans with >2/3 significant linear regression of ECoG to to handshape
% relevant_chans = sig_chans;
% %TRIM BAD CHANNELS:
% is_good = is_good_chan(ERPs); % Finds the not-bad channels of Recording
% relevant_chans = sig_chans & is_good;
% 
% ecog_test = ecog(relevant_chans,:,:);
% dims = size(ecog_test);
% 
% if isempty(components) | ((components == dims(1))|(components == 0)) 
%     ecog_data = ecog_test; % Does ignores PCA depending on user imput
% else    
%     ecog_shift = reshape(ecog_test,dims(1),(dims(2).*dims(3)));    % crunch timepts and trails together for PCA
%     ecog_shift = ecog_shift';    % flips data s.t. channels are cols
%     pca_mat = pca(ecog_shift);
%     %components = 20;        % select number of components. Must be >1 & <chans
%     pca_data = ecog_shift*pca_mat(:,1:components);
%     pca_data = pca_data';
%     ecog_data = reshape(pca_data, components, dims(2), dims(3)); % switch back into traditional form
% end
% 
% %% Remove uncategorized rows from Data:
% category_data = category_data(Is_categorized,:);
% ecog_data = ecog_data(:,:,Is_categorized);
% groupings = groupings(Is_categorized);
% 
% %% Loop through and make a series of Confusion Matricies:
% %Analysis_frame = 10;        % Analysis dones on average HG over 200ms window
% time_advance = 1;           % Advance the frame by 50ms each iteration
% total_time_pts = size(ecog,2);
% frames = floor((total_time_pts - Analysis_frame)/time_advance); % clippout the last frame
% 
% 
% confusion_mat = zeros(length(categories),length(categories),frames);
% SVMModels = [];     % Initializing this would be prudent...
% 
% for k = 1:frames
%     frame_num = ['Frame Number: ', num2str(k)];
%     save('status.mat','frame_num') % Print frame number - helpful to track runtime progress
%     time_start = (1+time_advance*(k-1)); % Starts at index = 1
%     time_end = time_start + Analysis_frame;
%     timewin = time_start:time_end;
%     relevant_timepts= timewin;
%     
% %     reps = length(timewin);
% %     X_data = [];
% %     for i = timewin
% %        datapt = squeeze(ecog_data(:,i,:));
% %         datapt = datapt';
% %         X_data = [X_data; datapt];
% %     end
%     
%     ecog_mean = squeeze(mean(ecog_data(:,relevant_timepts,:),2)); % Average HG in TimeWindow
%     X_data = ecog_mean';        % Important For formatting
% 
% 
%     %% Create SVM model for is/isnot category for each of the categories
% if k == 1
%     for i = 1:length(categories)
%         %Loops through for Each Category and Creates a kFoldPredict Struct
%         %for that class using a Linear SVM
%         %BoxConstraint can be lowered to precent over fitting;   Kernel Might Be Adjusted to 'rbf'
%         
%         Y_data = category_data(:,i);
%         %Y_data = repmat(Y_data, reps,1);
%         
%         SVMModel = fitcsvm(X_data,Y_data,'KernelFunction','linear', 'ClassNames',logical([0,1]), 'BoxConstraint', 0.01, 'Leaveout', 'on');
%         SVMModels{i} = SVMModel;
%     end
% else
%     for i = 1:length(categories)
%         Y_data = category_data(:,i);
%         %Y_data = repmat(Y_data, reps,1);
%         SVMModel = fitcsvm(X_data,Y_data,'KernelFunction','linear', 'ClassNames',logical([0,1]), 'BoxConstraint', 0.01, 'Leaveout', 'on');
%         SVMModels{i} = SVMModel;
%     end
% end
% 
%     %% GET Predictions from Classifier
%     predictions = zeros(size(groupings));
%     scores = zeros(size(X_data,1), length(categories)); % Input data_pts x categories matrix of Classifier Scores
%     for j = 1:length(categories)
%         [label,score] = kfoldPredict(SVMModels{j});
%         %[~,score] = predict(SVMModels{j},X_data_test);   
%         scores(:,j) = score(:,2);    % Classifier score of TRUE for each data point for category j    
%         labels(:,j) = label;
%     end
%     [~,ind]=max(scores,[],2);   % Classify to the classifier with the highest score
%     predictions = ind;
%     
%     %% Make ConfusionMatrix for TimeFrame
%      %groups = repmat(groupings, reps,1 );
%     groups = groupings;
%     confusion_mat_temp = confusionmat(groups, predictions);
%     normalize_confusion = false; % Normalizes the Rows of the ConfusionMatrix to 1
%     if normalize_confusion    
%         confusion_mat_temp = gdivide(confusion_mat_temp,sum(confusion_mat_temp,2)); % normalize Confusion Matrix
%     end
%     confusion_mat(:,:,k) = confusion_mat_temp;
%     
%     %Optional Print Positive ID accuracy (Precision)
%     accuracy(k) = trace(confusion_mat_temp)/sum(confusion_mat_temp(:));
% end 
% output = confusion_mat;
% end
% 
% 
% end
%     
%   
% 
% end
%     