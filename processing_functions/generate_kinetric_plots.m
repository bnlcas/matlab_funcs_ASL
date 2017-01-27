function [] = generate_kinetric_plots(ERPs)
%% Plot Graphics for Kinetic Features equalized for handshape
plot_full_grids = false;
plot_brain_maps = true;

twin = 170:230; % Relevant Area to measure difference


is_good = is_good_trial(ERPs);
is_lex = strcmpi(ERPs.annot.filledLexTrans, 'lexical');
Data_Tag = is_good & is_lex & ERPs.grouped_fill_annot.is_unique_erp;


ecog = ERPs.ecog;

%% Find Significant Channels that do not encode location
regress_grid_null = Regress_ECoG(ERPs, Data_Tag & get_ling_onset(ERPs), true);
sig_chans = get_sig_chans_ii(regress_grid_null);

%% Segment and find channels which are significant but have null location betas
is_chest_pre = false(size(sig_chans,1),1);
is_face_pre = false(size(sig_chans,1),1);
is_hand_pre = false(size(sig_chans,1),1);
is_fs_pre = false(size(sig_chans,1),1);
is_neut_pre = false(size(sig_chans,1),1);


sig_chan_list = find(sig_chans);
for i = 1:length(sig_chan_list)
    [~,peak_r2_ind] = max(squeeze(regress_grid_null(sig_chan_list(i),:,(end-1))));
    relevant_betas = squeeze(regress_grid_null(sig_chan_list(i),peak_r2_ind,21:25)); % Location Beta Weights at Peak R2 in a sig Chan

    is_chest_pre(sig_chan_list(i)) = ~(relevant_betas(1)==0);
    is_face_pre(sig_chan_list(i)) = ~(relevant_betas(2)==0);
    is_hand_pre(sig_chan_list(i)) = ~(relevant_betas(3)==0);
    is_fs_pre(sig_chan_list(i)) = ~(relevant_betas(4)==0);
    is_neut_pre(sig_chan_list(i)) = ~(relevant_betas(5)==0);
end
is_sig_non_loc = ~(is_face_pre | is_hand_pre | is_fs_pre | is_neut_pre) & sig_chans;

Data_Tag = is_good_trial(ERPs);
is_trans = Data_Tag & strcmpi(ERPs.annot.filledLexTrans,'transitional') & strcmpi(ERPs.grouped_fill_annot.grouping_hs,'onset');
is_lex = Data_Tag & strcmpi(ERPs.annot.filledLexTrans,'lexical') & strcmpi(ERPs.grouped_fill_annot.grouping_hs,'onset');

is_face = strcmpi(ERPs.grouped_fill_annot.loc,'face');
is_hand = strcmpi(ERPs.grouped_fill_annot.loc,'hand');
is_fs = strcmpi(ERPs.grouped_fill_annot.loc,'fingerspelling');
is_neut = strcmpi(ERPs.grouped_fill_annot.loc,'neutral');

%% Designate the Relevant Taggins of Lingusitic Data
is_ling = is_lex & (is_face | is_hand | is_fs | is_neut);   % Linguistic Data from the excluded locations
is_lex_lex = is_lex & (is_face | is_hand | is_neut);        % Lexical Data (fs excluded)
is_lex_fs = is_lex & is_fs;                                 % Fingerspellings only


%% Plot Trans vs Ling
comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag, 'handshape', is_trans, is_ling);

if plot_full_grids
    PlotECogGrid_Gen_Blackout(ERPs, true, ~is_sig_non_loc, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)))
    
    in_range = ERPs.time_axis <1000 & ERPs.time_axis > -1000;
    elect_trans = 179; elect_lex = 247;
%       elect_trans = 115; elect_lex = 202;
    figure; subplot(1,2,1)
    shadedErrorBar(ERPs.time_axis(in_range), mean(squeeze(ecog(elect_trans,in_range,comp_inds(:,1))),2), nansem(squeeze(ecog(elect_trans,in_range,comp_inds(:,1))),2),{'color', [0.8 0 0]},1)
    hold;
    shadedErrorBar(ERPs.time_axis(in_range), mean(squeeze(ecog(elect_trans,in_range,comp_inds(:,2))),2), nansem(squeeze(ecog(elect_trans,in_range,comp_inds(:,2))),2),{'color', [0 0 0.8]},1)
    xlabel('Time From Onset (ms)')
    ylabel('High Gamma Intensity')
    xlim([-1000 1000])
    title(['Comparison of Transitional (red) and Linguisitic (blue) ERPs in CH ' num2str(elect_trans)])
    plot([0 0], get(gca,'Ylim'),'k', 'LineWidth', 2)
    
    subplot(1,2,2)
    shadedErrorBar(ERPs.time_axis(in_range), mean(squeeze(ecog(elect_lex,in_range,comp_inds(:,1))),2), nansem(squeeze(ecog(elect_lex,in_range,comp_inds(:,1))),2),{'color', [0.8 0 0]},1)
    hold;
    shadedErrorBar(ERPs.time_axis(in_range), mean(squeeze(ecog(elect_lex,in_range,comp_inds(:,2))),2), nansem(squeeze(ecog(elect_lex,in_range,comp_inds(:,2))),2),{'color', [0 0 0.8]},1)
    plot([0 0], get(gca,'Ylim'),'k', 'LineWidth', 2)
    xlabel('Time From Onset (ms)')
    ylabel('High Gamma Intensity')
    xlim([-1000 1000])
    title(['Comparison of Transitional (red) and Linguisitic (blue) ERPs in CH ' num2str(elect_lex)])





end
%% Plot the relative difference between transitional and linguistic Data
if plot_brain_maps
    trans_ling_diff = mean(ERPs.ecog(:,twin,comp_inds(:,1)),3)-mean(ERPs.ecog(:,twin,comp_inds(:,2)),3);
    [~,max_ind] = max(abs(trans_ling_diff),[],2);
    for i = 1:256
        peak_diff(i) = trans_ling_diff(i,max_ind(i));
    end
    tmp = double(is_sig_non_loc).*peak_diff';
%     plot_brain_elecs_sandbox_greycbar(tmp)
%     title({'Maximum Difference in High Gamma Intensity of Linguistic and';'Transitional Gestures in Significant non-Location Electrodes'})


    %% Plot Trans-Lex Difference on the Brain:
    %comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag, 'handshape', is_trans, is_ling);

    twin = 170:230; 
    trans_mean_smth = mean(ERPs.ecog(:,twin,comp_inds(:,1)),3);
    lex_mean_smth = mean(ERPs.ecog(:,twin,comp_inds(:,2)),3);
    trans_lex_diff = trans_mean_smth - lex_mean_smth;

    [~,max_ind] = max(abs(trans_lex_diff),[],2);
    for i = 1:256
        peak_diff(i) = trans_lex_diff(i,max_ind(i));
    end
    tmp = -double(is_sig_non_loc).*peak_diff';
    lex_trans_diff = tmp(is_sig_non_loc);

%      plot_brain_elecs_sandbox_greycbar(tmp)
%      title({'Maximum Difference in High Gamma Intensity of Linguistic and';'Transitional Gestures in Significant non-Location Electrodes'})

    %% Plot Linguistic_Latency
%     [~,peak_time] = max(mean(ERPs.ecog(:,twin,is_lex_sig_loc),3),[],2);
%    [~,peak_time] = min(abs(diff(smooth_grid_time(mean(ERPs.ecog(:,twin,is_lex_sig_loc),3),20),1,2)),[],2);

   % peak_time = (peak_time-30)*10;
    peak_time = (max_ind - 30)*10;
    tmp = double(is_sig_non_loc).*peak_time;
    latency = tmp(is_sig_non_loc);
%     plot_brain_elecs_sandbox_greycbar(tmp)
%     title({'Latency Time of Peak High Gamma for Linguistic Gestures'; 'in Significant non-Location Electrodes'})
    %make_two_type_cbar(min(tmp), max(tmp))
end 
 
%% Plot Latency vs TransLex Difference:
fit_line_wparams(latency,lex_trans_diff)
xlabel('Latency of Peak Difference (ms)')
ylabel('Peak Difference in Linguistic and Transitional High Gamma Intensity')
title({'Plot of Peak Difference in Lexical and Transitional';'ERPs vs Latency Time of this Peak'})
hold;
plot([0 0],get(gca,'Ylim'),'k')
plot(get(gca,'XLim'),[0 0],'k')


%% Split Lex & Fingerspelling:
comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag, 'handshape', is_trans, is_lex_lex, is_lex_fs);
if plot_full_grids
    PlotECogGrid_Gen_Blackout(ERPs, true, ~is_sig_non_loc, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)),  ERPs.ecog(:,:,comp_inds(:,3)))
end
    
    

%% Addenda and Erratera:

is_uncontact = strcmpi(ERPs.annot.contact, 'uncontacted');
is_no_contact = strcmpi(ERPs.annot.contact, 'nocontact');
is_contact = strcmpi(ERPs.annot.contact,'contact');
is_contacting = is_contact | is_uncontact;

% Compare FS and LEX ERPs: with SVM or LDA
comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag, 'handshape', is_lex_fs, is_lex_lex);

twin = 190:210;
ecog_fs = ecog(is_sig_non_loc,twin,comp_inds(:,1));
ecog_lex = ecog(is_sig_non_loc,twin,comp_inds(:,2));
dat_fs_lex = [squeeze(mean(ecog_fs(:,190:210,:),2))'; squeeze(mean(ecog_lex(:,190:210,:),2))'];
%cat_fs_lex = [repmat({'fs'},436,1); repmat({'lex'},436,1)];
cat_fs_lex = [1*ones(size(comp_inds,1),1); 2*ones(size(comp_inds,1),1)];
SVM_model = fitcsvm(dat_fs_lex, cat_fs_lex);
%tmp(is_sig_non_loc) = abs(lda.Coeffs(1,2).Linear);
tmp(is_sig_non_loc) = abs(SVM_model.Beta);


end