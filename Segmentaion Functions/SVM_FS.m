function [output] = SVM_FS(ERPs, sig_chans)
%confusion_mat
output_confusion_mat = true;    % output is confusion_mat when true, and categories otherwise
output = [];
%% This function attempts to classify the handshapes of various handshapes of respones.
% the sig chans is currently a boolean 256x1 array that is TRUE for
% channels that are relevant to handshape segmentation.

% restrict signs to responses

Data_Tag = strcmpi(ERPs.annot.filledLexTrans, 'lexical');
ecog = ERPs.ecog(:,:,Data_Tag);

categories = {'fingerspelling'; 'lexical'};
class1 = strcmpi(ERPs.annot.respType(Data_Tag), 'fs');
class2 = ~strcmpi(ERPs.annot.respType(Data_Tag), 'fs');

%% Equalize for Handshape
eq_inds = Equalize_Parameter_Distribution(ERPs, Data_Tag, class1, class1);

class1 = class1(eq_inds);
class2 = class2(eq_inds);
ecog = ecog(:,:,eq_inds);

category_data = [class1, class2];





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

%% Loop through and make a series of Confusion Matricis:

Analysis_frame = 20;        % Average over 200ms window
time_advance = 5;           % Advance the frame by 50ms each time
total_time_pts = size(ecog,2);
frames = floor((total_time_pts - Analysis_frame)/time_advance); % clippout the last frame
%confusion_mat = zeros(length(category_data),length(category_data),frames);



%% Remove Data From Excluded HandShapes:
groupings = zeros(size(category_data,1),1);
for i = 1:length(categories)
    groupings = groupings + i.*category_data(:,i); % the ith classifier is listed as an i in a vector
end
exclude = (groupings == 0);
groupings = groupings(~exclude);
ecog_data = ecog_test(:,:,~exclude);
category_data = category_data(~exclude,:);



for k = 1:frames
    time_start = (1+time_advance*(k-1));
    time_end = time_start + Analysis_frame;
    timewin = time_start:time_end;
    k
     relevant_timepts = timewin;

    ecog_mean = squeeze(mean(ecog_data(:,relevant_timepts,:),2));

    X_data = ecog_mean';


    %% Create SVM model for is/isnot category for each of the categories
    SVMModels = [];     % Initializing this would be prudent...
    for i = 1:length(categories)
        Y_data = category_data(:,i);

        SVMModel = fitcsvm(X_data,Y_data,'KernelFunction','linear', 'ClassNames',logical([0,1]), 'BoxConstraint', 0.01, 'KernelScale', 1, 'KFold', 10);
        SVMModels{i} = SVMModel;

    end

    %% GET Predictions from Classifier
    predictions = zeros(size(groupings));
    scores = zeros(size(X_data,1), length(categories)); % Input data_pts x categories matrix of Classifier Scores
    for j = 1:length(categories)
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


    