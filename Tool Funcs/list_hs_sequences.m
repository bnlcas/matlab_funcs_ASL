function [] = list_hs_sequences(ERPs);
%%
%% Data Restriction
Data_Tag = is_good_trial(ERPs) & strcmpi(ERPs.annot.filledLexTrans,'lexical');
is_fs = Data_Tag & strcmpi(ERPs.grouped_fill_annot.loc,'fingerspelling');

is_fs(end) = [];

handshapes = ERPs.grouped_fill_annot.handshape;
hs_sequence = [handshapes(1:end-1), handshapes(2:end)];
hs_sequence_names = strcat(handshapes(1:end-1), handshapes(2:end)); % Needed to name each sequence

sequences = unique(hs_sequence_names(is_fs));
sequence_sizes = zeros(size(sequences));
for i = 1:length(sequences)
    sequence_sizes(i) = sum(strcmpi(hs_sequence_names(is_fs), sequences(i)));
end

a = 1;


[~,sort_order] = sort(sequence_sizes,'descend');
sequence_sizes = sequence_sizes(sort_order);
sequences = sequences(sort_order);

%% Partition E-R from T-R handshapes in FingerSpelling

is_ER = strcmpi(hs_sequence_names, 'ER') & is_fs;
is_TR = strcmpi(hs_sequence_names,'laxR') & is_fs;

inds_er = find(is_ER) + 1;
inds_tr = find(is_TR) + 1; % must shift the index over by one to capture the R

PlotECogGrid_Gen(ERPs, true, ERPs.ecog(:,:,inds_er), ERPs.ecog(:,:,inds_tr));
