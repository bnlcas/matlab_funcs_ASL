function [ERP_Frequency_Table] = compare_Onset_with_eqParam_Table(ERPs, compare_param)
%% Creates a Table of the number of onset ERPs that are possible for comparison of each combination of 2 tags in compare_param)
% The comparison is assumed to be between Location with HS equalized or HS
% with location equalized. This is controled through the user input
% 'compare_param', as either 'handshape', or 'loc'.

%% Restrict to ERPs of good trials and lexical communication
Data_Tag = is_good_trial(ERPs) & strcmpi(ERPs.annot.filledLexTrans,'lexical');
Data_Tag = Data_Tag & ERPs.grouped_fill_annot.is_unique_erp;

is_onset = strcmpi(ERPs.annot.lexTrans,'lexical');
if strcmpi(compare_param, 'handshape')
    comp_data = ERPs.grouped_annot.handshape;
    eq_data_onset = ERPs.grouped_annot.loc;
    eq_data_ongoing = ERPs.filled_change_annot.loc;
    
elseif strcmpi(compare_param, 'loc')
    comp_data = ERPs.grouped_annot.loc;
    eq_data_onset = ERPs.grouped_annot.handshape;
    eq_data_ongoing = ERPs.filled_change_annot.handshape;
    
    
    comp_data = ERPs.grouped_fill_annot.loc;
    eq_data_onset = ERPs.grouped_fill_annot.handshape;
    eq_data_ongoing = ERPs.grouped_fill_annot.handshape;
else
    'ERROR!, enter either ''loc'' or ''handshape'', or change the damn code!'
end

%% GetList of unique and relevant equalize params for column values:
equalizers = unique(eq_data_onset(Data_Tag));
equalizers(strcmpi(equalizers, '')) = [];       % kill unwanted labels
equalizers(strcmpi(equalizers,'changing')) = [];
equalizers(strcmpi(equalizers, 'lax')) = [];

%% Restrict the list of equalizers to ones with > Thresh occurances
eq_size = zeros(length(equalizers),1);
eq_thresh = 40;
for i = 1:length(equalizers)
    eq_size(i) = sum(strcmpi(eq_data_onset(Data_Tag), equalizers(i)));
end
equalizers(eq_size < eq_thresh) = [];


%% Assemble List of Comparison Tags
comparators = unique(comp_data(Data_Tag));
comparators(strcmpi(comparators, '')) = [];
comparators(strcmpi(comparators, 'lax')) = [];
comparators(strcmpi(comparators, 'changing')) = [];

%% Restrict to compartors with >thresh instances
comp_size = zeros(length(comparators),1);
comp_thresh = 50;
for i = 1:length(comparators)
    comp_size(i) = sum(strcmpi(comp_data(Data_Tag), comparators(i)));
end
comparators(comp_size<comp_thresh) = [];


%% Run through Each Unique Combination of Comparators and find the number
% Of equalizable ERPs for this combination.
count = 1;
comparator_combs = [];
for i = 1:length(comparators)
    for j = (i+1):length(comparators)
        comparator_combs{count} = [comparators{i}, ' - ', comparators{j}];
        count = count+1;
    end
end
data_table.comparison = comparator_combs';


%% Find the Number of comparisons between each combination of HS for each location
for i = 1:length(equalizers);
    eq_name = equalizers{i};
    
    % Create separate categories of equalization data - one for onset erps,
    % another for erps with an ongoing tag
    is_eq_onset = Data_Tag & strcmpi(eq_data_onset, equalizers(i));
    is_eq_ongoing = Data_Tag & ~ is_eq_onset & strcmpi(eq_data_ongoing, equalizers(i));
    
    is_eq_onset = Data_Tag & strcmpi(ERPs.grouped_fill_annot.grouping_hs,'onset') & strcmpi(eq_data_onset, equalizers(i));
    is_eq_ongoing = Data_Tag & strcmpi(ERPs.grouped_fill_annot.grouping_hs,'sustained') & strcmpi(eq_data_onset, equalizers(i));
    
    count = 1;
    for k = 1:length(comparators)
        for j = (k+1):length(comparators)
            % Create is_Ongoing and is_onset booleans for hs(k) and hs(j)
            comp_1 = comparators(k);
            is_comp1_onset = Data_Tag & strcmpi(comp_data , comp_1);
             
            comp_2 = comparators(j);
            is_comp2_onset = Data_Tag & strcmpi(comp_data, comp_2);
                        
            Comp12_on_Eq = Equalize_Tag_Sizes(is_comp1_onset & is_eq_onset, is_comp2_onset & is_eq_onset) |  Equalize_Tag_Sizes(is_comp1_onset & is_eq_ongoing, is_comp2_onset & is_eq_ongoing);
            comp_comb_count(count) = sum(Comp12_on_Eq(:,1));
            count = count+1;
            if strcmpi(eq_name,'hand') & strcmpi(comp_1,'O') & strcmpi(comp_2,'S')
                 PlotECogGrid_Gen(ERPs, true, ERPs.ecog(:,:,Comp12_on_Eq(:,1)), ERPs.ecog(:,:,Comp12_on_Eq(:,2)))
            end
        end
    end
    if size(comp_comb_count,1) < size(comp_comb_count,2)
        comp_comb_count = comp_comb_count';
    end
    
    %% Correct Equalizer Names to be valid column header names:
    if strcmpi(eq_name,'5')
        eq_name = 'Five';
    elseif strcmpi(eq_name,'1')
        eq_name = 'One';
    elseif strcmpi(eq_name,'K/P')
        eq_name = 'KP';
    elseif strcmpi(eq_name, '5>')
        eq_name = 'Bent_5';
    end
    
    eval(['data_table.' eq_name,' = comp_comb_count;']);

end

ERP_Frequency_Table = struct2table(data_table);
