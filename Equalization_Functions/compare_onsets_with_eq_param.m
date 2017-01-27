function [comp_inds] = compare_onsets(ERPs, comp1, comp2, eq_param, comp_type)
%% This function takes two data tags comp1 and comp2 and equalizes returns
% indecies at which they occur such that both comp1 and comp2 tags have an
% equal number of onset and sustained co-responses to eq_param
% assumes that comp and eq are either Handshape of Location

%% Restrict to Good Data
is_good = is_good_trial(ERPs);
is_lex = strcmpi(ERPs.annot.filledLexTrans,'lexical');
Data_Tag = is_good & is_lex & ERPs.grouped_fill_annot.is_unique_erp;


%% Define Data Taggins
if strcmpi(comp_type, 'loc')
    comp_tag = ERPs.grouped_fill_annot.loc;
    comp_grouping = ERPs.grouped_fill_annot.grouping_loc;
    
    eq_tag = ERPs.grouped_fill_annot.handshape;
    eq_grouping = ERPs.grouped_fill_annot.grouping_hs;
elseif strcmpi(comp_type, 'handshape')
    comp_tag = ERPs.grouped_fill_annot.handshape;
    comp_grouping = ERPs.grouped_fill_annot.grouping_hs;
    
    eq_tag = ERPs.grouped_fill_annot.loc;
    eq_grouping = ERPs.grouped_fill_annot.grouping_loc;
end

%% Define Booleans is Comp1 and 2 such that they are on for unique onset ERPs with the comp tagging that are good channels

is_comp1 = strcmpi(comp_tag,comp1)...
    & strcmpi(comp_grouping,'onset')...
        & ERPs.grouped_fill_annot.is_unique_erp & Data_Tag;
    
is_comp2 = strcmpi(comp_tag,comp2)...
    & strcmpi(comp_grouping,'onset')...
        & ERPs.grouped_fill_annot.is_unique_erp & Data_Tag;
    
%% Define Boolean Eq Onset and Sustain ERPs
   
is_eq_onset = strcmpi(eq_tag,eq_param)...
    & strcmpi(eq_grouping,'onset');
is_eq_ongoing = strcmpi(eq_tag,eq_param)...
    & strcmpi(eq_grouping,'sustained');

%% Get Equalized Taggings

comp_inds = Equalize_Tag_Sizes(is_comp1 & is_eq_onset, is_comp2 & is_eq_onset)...
    |  Equalize_Tag_Sizes(is_comp1 & is_eq_ongoing, is_comp2 & is_eq_ongoing);
