function [output] = SVM_Betas(ERPs, sig_chans, occurances, components, Analysis_frame)
%% This function attempts Runs SVM classification on the ERPs in sig_chans on classes of data which ocucur >occuraces times
% and then returns a Num(Sig_Chans) x Num_Categories x Timepts listing of the vector of the
% the hyperplane involved in the partitions among the classes involved
% the hyperplanes have been chosen without using crossvalidation, since
% this would prohibit the use of a single hyperplane.
%% Preliminary Processing:
output_confusion_mat = true;    % output is confusion_mat when true, and categories otherwise - useful for obtaining a list of categories used in classificaiton
output = [];

% restrict signs to lexical responses - This could be changed
Data_Tag = is_good_trial(ERPs) & (strcmpi(ERPs.annot.filledLexTrans, 'lexical')) & (ERPs.annot.dur_samp > 0);
%Data_Tag = (strcmpi(ERPs.annot.respType, 'fs')) & (ERPs.annot.dur_samp > 0);
ecog = ERPs.ecog(:,:,Data_Tag);


%% Get Annotation of Handshapes involved in duplicates
%sign_data = ERPs.annot.handshape(Data_Tag); % Segments on the Basis of
%sign_data = ERPs.annot.handshape(Data_Tag); % segment with location
sign_data = ERPs.annot.intMov(Data_Tag);
%sign_data = ERPs.annot.respType(Data_Tag);

sign_data = strrep(sign_data,' ','');               % Removes blank spaces
categories = unique(sign_data);
categories = categories(~strcmp(categories,''));    % Clear blanks from categories
sparse_category = false(size(categories));          % Will remove categories with too few instances
category_size_thresh = occurances;                           % Min of Three instances per category
cat_count = zeros(1,length(categories));
for i = 1:length(categories)
    cat_count(i) = sum(strcmpi(sign_data, categories(i)));
    if cat_count(i) < category_size_thresh
        sparse_category(i) = true;
    end
end
cat_count = cat_count(~sparse_category);
categories = categories(~sparse_category);

%%Eliminate Changing and Lax Categoreis
cat_count = cat_count(~(strcmpi(categories, 'changing') | strcmpi(categories, 'lax')));
categories = categories(~strcmpi(categories, 'changing'));
categories = categories(~strcmpi(categories, 'lax'));
%cat_count(find(strcmpi(categories, 'neutral'))) = [];
%categories(find(strcmpi(categories, 'neutral'))) = [];


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
Randomize_Data = false;
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
    pca_mat = pca(ecog_shift);
    %components = 20;        % select number of components. Must be >1 & <chans
    pca_data = ecog_shift*pca_mat(:,1:components);
    pca_data = pca_data';
    ecog_data = reshape(pca_data, components, dims(2), dims(3)); % switch back into traditional form
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
%Analysis_frame = 10;        % Analysis dones on average HG over 100ms
time_advance = 5;           % Advance the frame by 50ms each iteration
total_time_pts = size(ecog,2);
frames = floor((total_time_pts - Analysis_frame)/time_advance); % clippout the last frame


%SVM_betas = zeros(length(sig_chans), length(categories),frames);
SVM_betas = zeros(length(sig_chans), length(categories));

%for k = 1:frames

    k = 40% Print frame number - helpful to track runtime progress
    
    time_start = (1+time_advance*(k-1)); % Starts at index = 1
    time_end = time_start + Analysis_frame;
    relevant_timepts= time_start:time_end;
        
    ecog_mean = squeeze(mean(ecog_data(:,relevant_timepts,:),2)); % Average HG in TimeWindow
    X_data = ecog_mean';        % Important For formatting
    

    %% Create SVM model for is/isnot category for each of the categories

%    scores = zeros(size(category_data));
    for i = 1:length(categories)
        %Loops through for Each Category and Creates a kFoldPredict Struct
        %for that class using a Linear SVM
        %BoxConstraint can be lowered to precent over fitting;   Kernel Might Be Adjusted to 'rbf'
        
        Y_data = category_data(:,i);       
        
%         SVMModel = fitcsvm(X_data,Y_data,'KernelFunction','linear', 'ClassNames',logical([0,1]), 'BoxConstraint', 0.01, 'KernelScale',1);      
%         weights_single = SVMModel.Beta;
%         SVM_betas(sig_chans,i,k) = weights_single.^2;  % take the square of the beta weights into the sig_chans
        
        %% ALT CODE - AVERAGE OVER LEAVEOUT
        data_mask = true(size(Y_data));
        weights = zeros(length(Y_data), sum(sig_chans));
        
        if k>=1
            for j = 1:length(Y_data)
                data_mask(j) = false;
                SVMModel = fitcsvm(X_data(data_mask,:),Y_data(data_mask),'KernelFunction','linear', 'ClassNames',logical([0,1]), 'BoxConstraint', 0.01, 'KernelScale',1);
                weights(j,:) = SVMModel.Beta;
                data_mask(j) = true;
                
%                 [~,score] = predict(SVMModel, X_data(j,:));
%                 scores(j,i) = score(2);
            end
        else
            parfor j = 1:length(Y_data)
                tmp_mask = data_mask;
                tmp_mask(j) = false;
                SVMModel = fitcsvm(X_data(tmp_mask,:),Y_data(tmp_mask),'KernelFunction','linear', 'ClassNames',logical([0,1]), 'BoxConstraint', 0.01, 'KernelScale',1);
                weights(j,:) = SVMModel.Beta;
            end
        end
        SVM_betas(sig_chans,i) = mean(weights,1).^2;
        
        
    end
%     [~,ind] = max(scores, [], 2);
%     is_cor = false(size(ind));
%     for i = 1:length(ind)
%         is_cor(i) = (category_data(i, ind(i)) ~=0);
%     end
%     accuracy = correct/length(ind)
end



output = SVM_betas;

end

% end