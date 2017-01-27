function [output] = Segment_Handshapes_Slide_Cleaned(ERPs, Data_Tag, sig_chans, occurances, rand, components, loc_hs)
%% This function attempts to classify the handshapes of various handshapes of respones from the (ERPs) structure
% and create a sequence of confusion matricies charting the time evolution
% of this classification scheme.
% Only handshapes which occur (Ocurrances) number of times are included in
% the Analysis. The anaysis uses a one-vs-all classification on the relevant ERPs 
%
% The feature vector for each handshape is taken from (sig_chans) and is
% given as the Average ECoG HG response in a specified timewindow.
%
% It is possible to reduce the total number of feature vectors by performing PCA
% on the Data_Matrix. This can modulated by the input (components) which specifies the number of PCA components to use.
% If components is empty or set to 0 or set to the full number of components, then the program
% will simply use the full ecog data for the analysis
%% Preliminary Processing:
output_confusion_mat = true;    % output is confusion_mat when true, and categories otherwise - useful for obtaining a list of categories used in classificaiton
output = [];

%rng('shuffle');
%% restrict signs to duplicate responses - This could be changed
%Data_Tag = strcmpi(ERPs.annot.respType,'dup');
%Data_Tag = ~strcmpi(ERPs.annot.handshape,'');
%Data_Tag = (~strcmpi(ERPs.annot.handshape, '') & ~strcmpi(ERPs.annot.respType, 'fs'));
%Data_Tag = is_good_trial(ERPs) & (strcmpi(ERPs.annot.filledLexTrans, 'lexical')) & (ERPs.annot.dur_samp > 0);
%Data_Tag = (strcmpi(ERPs.annot.respType, 'fs')) & (ERPs.annot.dur_samp > 0);
ecog = ERPs.ecog(:,:,Data_Tag);


%% Get Annotation of Handshapes involved in duplicates
use_loc = false;    use_handshape = false;    use_intmov = false;
if substrcmp(loc_hs,'loc')
    use_loc = true;
end
if substrcmp(loc_hs,'handshape')
    use_handshape = true;
end
if strcmpi(loc_hs, 'intmov')
    use_intmov = true;
end


if use_loc
    sign_data = ERPs.grouped_fill_annot.loc(Data_Tag); % Segments on the Basis of
end
if use_handshape
    sign_data = ERPs.grouped_fill_annot.handshape(Data_Tag); % segment with location
end
if use_intmov
    sign_data = ERPs.grouped_fill_annot.intMov(Data_Tag);
end


sign_data = strrep(sign_data,' ','');               % Removes blank spaces
categories = unique(sign_data);
categories = categories(~strcmp(categories,''));    % Clear blanks from categories
sparse_category = false(size(categories));          % Will remove categories with too few instances
category_size_thresh = occurances;                           % Min of Three instances per category
cat_count = zeros(1,length(categories));
for i = 1:length(categories)
    cat_count(i) = sum(strcmpi(sign_data, categories(i)));
end
sparse_category = (cat_count<category_size_thresh);
cat_count = cat_count(~sparse_category);
categories = categories(~sparse_category);

%%Eliminate Changing and Lax Categoreis
cat_count = cat_count(~(strcmpi(categories, 'changing') | strcmpi(categories, 'lax')));

categories = categories(~strcmpi(categories, 'changing'));
categories = categories(~strcmpi(categories, 'lax'));

cat_count = cat_count(~strcmpi(categories,'fs'));
categories = categories(~strcmpi(categories,'fs'));
%cat_count(find(strcmpi(categories, 'neutral'))) = [];
%categories(find(strcmpi(categories, 'neutral'))) = [];

% categories = {'S'; 'A'}; %- useful to manually select categories
if ~output_confusion_mat
    output = categories;    % Return list of categories
    out_cat_size = false; % optional output of magnetude of each category (may be
    if out_cat_size
        output = cat_count;   
    end
else

    


%% Assemble Y data for Classifiers
% Create a num_handshapes x Num_Categories Matrix whose colums are data
% vectors containing info on whether a data point is in a category or not
category_data = false(length(sign_data), length(categories));
for i = 1:length(categories)
    category_data(:,i) = strcmp(sign_data, categories(i));   % Boolean Vector
end

groupings = double(category_data)*transpose(1:length(categories)); % numbering for each category (ith category is i...)
Is_categorized = (sum(category_data,2) ~= 0);   % boolean is false if a Data point is not categorized


%% Randomly Permute Category Data to test model
Randomize_Data = rand;
if Randomize_Data
    new_order = randperm(length(sign_data));
    category_data = category_data(new_order,:);
    groupings = groupings(new_order);
    Is_categorized = Is_categorized(new_order);
end


%% Assemble X data to to PCA for Classifier
% relevant channels is function call,
% taken here as sig_chans = (sum((stat_grid_dups_p(:,191:211,52)<0.01),2) > 12);
% or Chans with >2/3 significant linear regression of ECoG to to handshape
relevant_chans = sig_chans;
%TRIM BAD CHANNELS:
is_good = is_good_chan(ERPs); % Finds the not-bad channels of Recording
relevant_chans = sig_chans & is_good;

ecog_test = ecog(relevant_chans,:,:);
dims = size(ecog_test);

if isempty(components) | ((components == dims(1))|(components == 0)) 
    ecog_data = ecog_test; % Does ignores PCA depending on user imput
else    
    ecog_shift = reshape(ecog_test,dims(1),(dims(2).*dims(3)));    % crunch timepts and trails together for PCA
    ecog_shift = ecog_shift';    % flips data s.t. channels are cols
    ecog_shift = gsubtract(ecog_shift, mean(ecog_shift,2)); % Center the Data
    pca_mat = pca(ecog_shift);
    %components = 20;        % select number of components. Must be >1 & <chans
    pca_data = ecog_shift*pca_mat;
    pca_data = pca_data';
    ecog_data = reshape(pca_data, dims(1), dims(2), dims(3)); % switch back into traditional form
    ecog_data = ecog_data(1:components,:,:);
end

%% Remove uncategorized rows from Data:
category_data = category_data(Is_categorized,:);
ecog_data = ecog_data(:,:,Is_categorized);
groupings = groupings(Is_categorized);

%% Randomly DownSample to even_sample_size
even_sample_size = min(cat_count)-5; % the -10 term helps mitigate the sampling issues on the smallest category
use_even_sample_size = false;
if use_even_sample_size
    Shuffle = randperm(length(groupings));
    category_data = category_data(Shuffle,:);
    ecog_data = ecog_data(:,:,Shuffle);
    groupings = groupings(Shuffle);
    for i = 1:length(categories)
        % total of 
        total_samples = sum(category_data(:,i));
        remove = find(category_data(:,i), (total_samples-even_sample_size));

        category_data(remove,:) = [];
        ecog_data(:,:,remove) = [];
        groupings(remove) = [];
    end
end

%% Loop through and make a series of Confusion Matricies:
Analysis_frame = 10;        % Analysis dones on average HG over 200ms window
time_advance = 2;           % Advance the frame by 50ms each iteration
total_time_pts = size(ecog,2);
frames = floor((total_time_pts - Analysis_frame)/time_advance); % clippout the last frame


confusion_mat = zeros(length(categories),length(categories),frames);
SVMModels = [];     % Initializing this would be prudent...
SVMModels_b = [];

for k = 1:frames
%for k = 50:50:150
    %k% Print frame number - helpful to track runtime progress
    
    time_start = (1+time_advance*(k-1)); % Starts at index = 1
    time_end = time_start + Analysis_frame;
    timewin = time_start:time_end;
    relevant_timepts= timewin;
    
    
%     reps = length(timewin);
%     X_data = [];
%     for i = timewin
%        datapt = squeeze(ecog_data(:,i,:));
%         datapt = datapt';
%         X_data = [X_data; datapt];
%     end
    
    ecog_mean = squeeze(mean(ecog_data(:,relevant_timepts,:),2)); % Average HG in TimeWindow
    X_data = ecog_mean';        % Important For formatting
    
    

    %% Create SVM model for is/isnot category for each of the categories
if k == 1
    for i = 1:length(categories)
        %Loops through for Each Category and Creates a kFoldPredict Struct
        %for that class using a Linear SVM
        %BoxConstraint can be lowered to precent over fitting;   Kernel Might Be Adjusted to 'rbf'
        
        Y_data = category_data(:,i);
        
        
        %SVMModel = fitcsvm(X_data,Y_data,'KernelFunction','linear', 'ClassNames',logical([0,1]), 'BoxConstraint', 0.01, 'KernelScale', 1, 'Leaveout', 'on');
        %SVMModel = fitcsvm(X_data,Y_data,'KernelFunction','linear', 'BoxConstraint', 0.01, 'KernelScale', 1, 'KFold', 10);

        SVMModel = fitcsvm(X_data,double(Y_data),'KernelFunction','linear', 'BoxConstraint', 0.01, 'KernelScale', 1, 'KFold', 10,'prior','uniform');

        
        SVMModels{i} = SVMModel;
    end
else
    parfor i = 1:length(categories)
        Y_data = category_data(:,i);
 %       SVMModel = fitcsvm(X_data,Y_data,'KernelFunction','linear', 'ClassNames',logical([0,1]), 'BoxConstraint', 0.01, 'KernelScale', 1, 'Leaveout', 'on');
        %SVMModel = fitcsvm(X_data,Y_data,'KernelFunction','linear', 'BoxConstraint', 0.01, 'KernelScale', 1, 'KFold', 10);

        SVMModel = fitcsvm(X_data,double(Y_data),'KernelFunction','linear', 'BoxConstraint', 0.01, 'KernelScale', 1, 'KFold', 10,'prior','uniform');

        SVMModels{i} = SVMModel;

    end
end



    %% GET Predictions from Classifier
    predictions = zeros(size(groupings));
    scores = zeros(size(X_data,1), length(categories)); % Input data_pts x categories matrix of Classifier Scores
    parfor j = 1:length(categories)
        [label,score] = kfoldPredict(SVMModels{j});
        %[~,score] = predict(SVMModels{j},X_data_test);   
        scores(:,j) = score(:,2);    % Classifier score of TRUE for each data point for category j    
        labels(:,j) = label;
    end
    [~,ind]=max(scores,[],2);   % Classify to the classifier with the highest score
    predictions = ind;
    
    %% Make ConfusionMatrix for TimeFrame
     %groups = repmat(groupings, reps,1 );
    groups = groupings;
    confusion_mat_temp = confusionmat(groups, predictions);
    normalize_confusion = false; % Normalizes the Rows of the ConfusionMatrix to 1
    if normalize_confusion    
        confusion_mat_temp = gdivide(confusion_mat_temp,sum(confusion_mat_temp,2)); % normalize Confusion Matrix
    end
    confusion_mat(:,:,k) = confusion_mat_temp;
    
    %Optional Print Positive ID accuracy (Precision)
    accuracy(k) = trace(confusion_mat_temp)/sum(confusion_mat_temp(:));
end 

output = confusion_mat;
end


end
    
    