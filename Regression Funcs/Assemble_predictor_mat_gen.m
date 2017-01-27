function [output] = Assemble_predictor_mat_gen(ERPs, Data_Tag, normalize, out_type)
%% Assembles Predictor Matrix
% Takes the ERP and creates regress_mat, a trials x predictor matrix which uses
% dummy variables for different parts of speech taggings and annotations.

% The tag 'normalize' scales the range of every predictor variable to [0,1]

% The tag 'out_type' specifies whether to output a regression data set, or
% a list of names. The default is regression matrix. To output names, set
% out_type = 'names'

%%Used Filled in Data Annotations:




%     annot = ERPs.filled_annot;
%     annot = ERPs.grouped_annot;
%     annot = ERPs.annot;
%annot = ERPs.grouped_annot;
annot = ERPs.grouped_fill_annot;
%annot = ERPs.annot;

%% Turn on/off Classes of predictors

% Event Related/Response
use_Timings = false;
use_Reaction_Times = false;
use_stimOnset = false;
use_movOnset = false;
use_TransLex = false;
use_FilledTransLex = false;

% Lexical/Word Frequency Stats:
use_stimType = false;
use_lex_freq = false;  % was true in first hs regree
use_lex_aoa = false;
use_english_freq = false;

% Duplicate, finger spelling etc
use_respType = false;
use_contact = false; % was TRUE is first hs regress

% Phonetic Parameters:
use_HandShape = true;
use_signType = false; % was TRUE in first hs reg

% Physical Action of gestures
use_loc = true;
use_movPath = false; % 1 ld col with mDir and intMov
use_movDir = false;  % movDir and inMov had an LD col
use_intMov = true; % movDir and int Mov are LI

apply_sparse_thresh = true;


num_trials = size(ERPs.ecog,3);
regress_mat = ones(num_trials,1);       % initialize regression matrix and provide an offset
predictor_var_names = [{'offset'}]; % list of variable names

%% Run Parameters (very few possible dimensions):

%Timing Related Predictors
if use_Timings
    durations = annot.dur_samp./1000;
    regress_mat = [regress_mat durations]; %Sample Duration
    predictor_var_names = [predictor_var_names; {'SampDur'}];
    
%     last_dur = [0; durations(2:end)];
%     regress_mat = [regress_mat last_dur]; %Prior Duration
%     predictor_var_names = [predictor_var_names; {'PriorSampDur'}];
end

if use_Reaction_Times
    Stim_On = strcmpi(annot.stimOnset,'S');
    react_mov = zeros(size(Stim_On));
    react_lex = zeros(size(Stim_On));
    start_times = annot.start_samp/10; % start times in (ms)
 
    Mov_On = strcmpi(annot.movOnset,'M');
    [Movement_Tag, Filled_Movement_Tag] = Get_First_Subsequent_Tag(Stim_On, Mov_On);
    Stim_start = start_times(Filled_Movement_Tag);
    Mov_start = start_times(Movement_Tag);
    react_mov(Filled_Movement_Tag) = Mov_start - Stim_start;
    
    Lex_On = strcmpi(annot.lexTrans,'lexical');
    [Lex_Tag, Filled_Lex_Tag] = Get_First_Subsequent_Tag(Stim_On, Lex_On);
    Stim_start = start_times(Filled_Lex_Tag);
    Lex_start = start_times(Lex_Tag);
    react_lex(Filled_Lex_Tag) = Lex_start - Stim_start;
    
    % Restrict Data Tag to Reaction Events
    Data_Tag = Data_Tag & Filled_Lex_Tag & Filled_Movement_Tag;
    
    
    regress_mat = [regress_mat react_mov react_lex];
    predictor_var_names = [predictor_var_names; {'React_Mov'}; {'React_Lex'}];
    
end

%Stimulus Onset
if use_stimOnset
    occurance1 = double(strcmpi(annot.stimOnset,'S'));
    regress_mat = [regress_mat occurance1];
    predictor_var_names = [predictor_var_names; {'StimOn'}];
end

%Movement Onset
if use_movOnset
    occurance1 = double(strcmpi(annot.movOnset,'M'));
    regress_mat = [regress_mat occurance1];
    predictor_var_names = [predictor_var_names; {'movOn'}];
end

%Trans/Lex
if use_TransLex
    occurance1 = double(strcmpi(annot.lexTrans, 'lexical'));
    occurance2 = double(strcmpi(annot.lexTrans, 'transitional'));
    regress_mat = [regress_mat occurance1 occurance2];
    predictor_var_names = [predictor_var_names; {'lex'}; {'trans'}];
end

%Filled in Trans/Lex
if use_FilledTransLex
    occurance1 = double(strcmpi(annot.filledLexTrans, 'lexical'));
    occurance2 = double(strcmpi(annot.filledLexTrans, 'transitional'));
    regress_mat = [regress_mat occurance1 occurance2];
    predictor_var_names = [predictor_var_names; {'lex_fill'}; {'trans_fill'}];
end

% Stim Type: filler; early acquired; late acquired; low frequency; high frequency
if use_stimType
    occurance1 = double(strcmpi(annot.stimType, 'f'));
    occurance2 = double(strcmpi(annot.stimType, 'ea'));
    occurance3 = double(strcmpi(annot.stimType, 'la'));
    occurance4 = double(strcmpi(annot.stimType, 'lf'));
    occurance5 = double(strcmpi(annot.stimType, 'hf'));
    regress_mat = [regress_mat occurance1 occurance2 occurance3 occurance4 occurance5];
    predictor_var_names = [predictor_var_names; {'filler_Stim'}; {'early acquire'}; {'late acquire'}; {'lowfreq'}; {'highfreq'}];
end
%Alternatively to Stim Type, use numeric data on word frequency.
if use_lex_freq
    occurance2 = annot.frequency;
    regress_mat = [regress_mat occurance2];
    predictor_var_names = [predictor_var_names; {'lex_freq'}];
end
if use_lex_aoa
    occurance3 = annot.aoa;
    regress_mat = [regress_mat occurance3];
    predictor_var_names = [predictor_var_names; {'aoa'}];

end

if use_english_freq
    occurance1 = annot.frequency_eng;
    %occurance1 = annot.frequency;

    regress_mat = [regress_mat occurance1];
    predictor_var_names = [predictor_var_names; {'eng_freq'}];
end

% Response Type: duplicated stimuli; fingerspelled; repair; Yes; No;
% comment
if use_respType
    occurance1 = double(strcmpi(annot.respType, 'dup'));
    occurance2 = double(strcmpi(annot.respType, 'fs'));
    occurance3 = double(strcmpi(annot.respType, 'repair'));
    occurance4 = double(strcmpi(annot.respType, 'YES'));
    occurance5 = double(strcmpi(annot.respType, 'NO'));
    occurance6 = double(strcmpi(annot.respType, 'comment'));
    regress_mat = [regress_mat occurance1 occurance2 occurance3 occurance4 occurance5 occurance6];
    predictor_var_names = [predictor_var_names; {'dup_stim'}; {'finger_sp'}; {'repair'}; {'yes_resp'}; {'no_resp'}; {'comment'}];
end


    
%% Arbitrarly large number of possible predictors:

%% 'Phonetic' Classifiers:

% Hand Shape (of response?) lax(relaxed); changing; (Also includes letter
% spelling representing signs.
if use_HandShape
    sign_data = annot.handshape;
    categories = unique(strrep(sign_data,' ',''));  % strrep used to remove spaces
    categories = categories(~strcmp(categories,''));           % Clear '' from categories
    categories = categories(~strcmpi(categories,'?'));
    categories = categories(~strcmpi(categories,'changing'));
    categories = categories(~strcmpi(categories,'lax'));
    
    sparse_category = false(size(categories));          % Will remove categories with too few instances
    category_size_thresh = 5;                           % Min of Three instances per category
    for i = 1:length(categories)
        % cat_count = sum(strcmpi(annot.respType, 'dup')  & strcmpi(sign_data, categories(i)));
        cat_count = sum(strcmpi(sign_data, categories(i)));
        if cat_count < category_size_thresh
            sparse_category(i) = true;
        end
    end
    categories = categories(~sparse_category);
    categories = categories(~strcmpi(categories,'changing')); % Clear Changing Handshape
    
    for i = 1:length(categories)
        occurance = double(strcmp(sign_data, categories(i)));
        regress_mat = [regress_mat occurance];
        predictor_var_names = [predictor_var_names; categories(i)];
    end
end

% Sign Type: Type of sign being gestured:
if use_signType
    sign_data = annot.signType;
    
    categories = unique(strrep(sign_data,' ',''));    % strrep used to remove spaces
    categories = categories(~strcmp(categories,''));  % Clear '' from categories
    
    for i = 1:length(categories)
        occurance = double(strcmpi(sign_data, categories(i)));
        regress_mat = [regress_mat occurance];
        predictor_var_names = [predictor_var_names; categories(i)];
    end
end

%% Physical Representation of Gestures (less hand based)

% Location of signing: neutral; chest; face; arm; hand; changing (also includes
% fingerspelling, but this seems redundant
if use_loc
    sign_data = annot.loc;
    categories = unique(strrep(sign_data,' ',''));  % strrep used to remove spaces
    categories = categories(~strcmp(categories,''));          % Clear '' from categories
    category_size_thresh = 5;                           % Min of Three instances per category
    for i = 1:length(categories)
        cat_count(i) = sum(strcmpi(sign_data, categories(i)));
    end
    
    categories = categories(cat_count>category_size_thresh);
    categories = categories(~strcmpi(categories,'changing'));
 
    
    
    for i = 1:length(categories)
        occurance = double(strcmpi(sign_data, categories(i)));
        regress_mat = [regress_mat occurance];
        predictor_var_names = [predictor_var_names; categories(i)];
    end
end

% MovementPath:arc = arched path; straight = straight path; circle = path traces a circle
if use_movPath
    sign_data = annot.movPath;
    categories = unique(strrep(sign_data,' ',''));  % strrep used to remove spaces
    categories = categories(~strcmp(categories,''));          % Clear '' from categories
    for i = 1:length(categories)
        occurance = double(strcmpi(sign_data, categories(i)));
        regress_mat = [regress_mat occurance];
        predictor_var_names = [predictor_var_names; categories(i)];
    end
end

    
% use_movDir:sup = superior; inf = inferior; ant = anterior; post = posterior; contra = contralateral to dominant side
%ipsi = ipsilateral to dominant side
if use_movDir
    sign_data = annot.movDir;
    categories = unique(strrep(sign_data,' ',''));  % strrep used to remove spaces
    categories = categories(~strcmp(categories,''));          % Clear '' from categories
    categories = categories(~strcmpi(categories,'repeated')); %
    for i = 1:length(categories)
        occurance = double(strcmpi(sign_data, categories(i)));
        regress_mat = [regress_mat occurance];
        predictor_var_names = [predictor_var_names; categories(i)];
    end
end

if use_intMov
    sign_data = annot.intMov;
    categories = unique(strrep(sign_data,' ',''));  % strrep used to remove spaces
    categories = categories(~strcmpi(categories,'claw')); % has only a single occurance
    categories = categories(~strcmp(categories,''));          % Clear '' from categories
    for i = 1:length(categories)
        occurance = double(strcmpi(sign_data, categories(i)));
        regress_mat = [regress_mat occurance];
        predictor_var_names = [predictor_var_names; categories(i)];
    end
end    

% Specify the sign called for contact:
% Note that contact and uncontacting are classed together, since they were
% seen to have nearly identical responses
if use_contact
    occurance1 = double(strcmpi(annot.contact, 'no contact'));
    occurance2 = double(strcmpi(annot.contact, 'contact') | strcmpi(annot.contact, 'uncontacted'));
    regress_mat = [regress_mat occurance1 occurance2];
    predictor_var_names = [predictor_var_names; {'no_contact'}; {'contact'}];
end  




%% Clean Data.
% Only select data points used under data Tag.
% Only use predictors that are used throughout
regress_mat = regress_mat(Data_Tag,:);

row_means = mean(regress_mat,1); % average value of each predictor
unused = ((row_means == 1) | (row_means == 0));
unused(1) = false; % preserve offset column
regress_mat = regress_mat(:, ~unused);
predictor_var_names = predictor_var_names(~unused);



%% Normalize Predictors
% This assumes that all predictors in positive reals
if normalize
    regress_mat = gdivide(regress_mat,max(regress_mat));
end
predictor_var_names;

%% Optional Clearing of Sparse Variables 
if apply_sparse_thresh
    sparse_thresh = 10;
    cat_sizes = sum(regress_mat,1);
    sparse_cat = (cat_sizes < sparse_thresh);
    predictor_var_names(sparse_cat) = [];
    regress_mat(:,sparse_cat) = [];
end

%% Return Output
if isempty(out_type)
    output = regress_mat;
elseif strcmpi(out_type, 'names')
    output = [predictor_var_names];
else
    output = regress_mat;
end

end
