function [stat_grid] = Regress_ECoG(ERPs, Data_Tag, varargin)
%% This function takes the ECoG data from the ERP and regresses
%% each_channel's time evolution to the values of the indicator parameters
%% listed in the MasterFile
%% Use the input data tag to select Events from the Master ERP of interest

% EX Data_Tag = strcmpi(ERPs.annot.stimOnset, 'S');
Data_Tag = Data_Tag & is_good_trial(ERPs);

if isempty(varargin)
    zero_null_betas = false;
else
    zero_null_betas = varargin{1};
end


if isempty(Data_Tag)
    Data_Tag = true(size(ERPs.annot.stimOnset)) & is_good(ERPs); % Regress Everything (all ERP trials) to Features
end



%% Get relevant ECOG DATA
ecog = ERPs.ecog(:,:,Data_Tag); % Y DATA
dims = size(ecog);

%% Perform Regression on dECoG/dt
regress_derivative = false;
if regress_derivative
    %loop through all trials and differentiate ECoG Grid w/resp2time
    ecog_diff = zeros(dims(1),dims(2),dims(3));
    for i = 1:dims(3)
        ecog_frame = squeeze(ecog(:,:,i));
        smoothing = 26;  % smoothing window length 
        fir = ones(1,smoothing)./smoothing;
        ecog_smth = filter(fir,1,ecog_frame, [],2);  %smooth in the second dimension (time)
        ecog_diff(:,1:(dims(2)-1),i)  = diff(ecog_smth,1,2);
    end
    ecog = ecog_diff;
    dims = size(ecog);
end

%%Assemble Predictor Matrix:

regress_mat = Assemble_predictor_mat_gen(ERPs, Data_Tag, true,'data'); % TRUE tag Normalizes the predictors
num_predictors = size(regress_mat,2);



%% Clear rows of regress that are all zeros (save the offset)
remove_null_predictors = true;
if remove_null_predictors
    row_total = sum(regress_mat,2);
    live_data = (row_total ~= 1); % includes offset - will giltch in case of superpositions
    
    ecog = ecog(:,:,live_data);
    regress_mat = regress_mat(live_data,:);
end

%% Collapse the Data to only make predictions for events occuring in
%% larger blocks (e.g. stimulus blocks or trans/lexical blocks)
collapse_predictors = false;
if collapse_predictors
 %    collapser = ERPs.annot.lexTrans;                % Collapse RegressMat to start of trans/Lex Events 
 %    block_on_trigger = {'transitional' 'lexical'};  % Trigger blocks on these identifies
    collapser = ERPs.annot.stimOnset;
    block_on_trigger = {'S'};
    regress_mat =  Collapse_Regression_Matrix(collapser, block_on_trigger, regress_mat);
    ecog = Collapse_ECOG_Data(collapser, block_on_trigger, ecog);
    size(regress_mat)
end
% tic
% [coeffs,~,~,~,stats] = regress(squeeze(ecog(100,250,:)), regress_mat);
% stat_grid = coeffs;
% toc
%% Loop through the Grid to regress each Channel and TimePt to predictors
% *** VERY Poor solution...
dims = size(ecog);
stat_grid = zeros(dims(1),dims(2), num_predictors+2);


for i = 1:dims(1) % Channels %  
    for j = 1:dims(2) % % 
        ecog_trial = squeeze(ecog(i,j,:));
        if zero_null_betas
            [coeffs,bint,~,~,stats]=regress(ecog_trial, regress_mat);
            is_null = (bint(:,1) < 0) & (bint(:,2) >0);
            coeffs = coeffs.*double(~is_null);
        else
            [coeffs,~,~,~,stats]=regress(ecog_trial, regress_mat);
        end
        stat_grid(i,j,1:num_predictors) = coeffs;
        stat_grid(i,j, 1+ num_predictors) = stats(1); % R^2
        stat_grid(i,j, 2+ num_predictors) = stats(3); % pval
    end
    %i % output to display progress
end


%% Clear Bad Channels
Bad_Channels = ERPs.BadChans{:,2};
stat_grid(Bad_Channels,:,1:(1+num_predictors)) = 0;
stat_grid(Bad_Channels,:,(2+num_predictors)) = 1; % set Pvals to max for badchans

end
