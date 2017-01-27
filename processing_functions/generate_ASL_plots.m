function [] = generate_ASL_plots(ERPs)
%% This function generates a sequence of graphs for presentation purposes

%% controls:
colorlist = [0.8 0 0; 0 0 0.8;0.8 0.8 0; 0 0.8 0.8; 0 0.8 0; 0.8 0 0.8];    % Plots of colors listed as Red, Blue, Yellow, Teal, Green, Purple

%% Plot Graphs for RAW Data
plot_raw = false;       
    plot_full_grids = false;
    plot_select_electrodes = false;
    plot_sig_chans = false;

    
    plot_face_hand_in5 = false;
    plot_O_5_neut = false;
    plot_fs_comp = false;
    
%% Plot Kinetic Features:
plot_kinetic = false;
    plot_tl_diff_brain = true;

%% Plot map of Regression
plot_regression = true; 

%% Generate Heirarchical Clustering
plot_clustering = false; 

%% Generate Plots for SVM segmentation
plot_segmentation = false;   
    segment_loc = true;
    segment_hs = true;
    create_confusion_mats = true;

%% Restrict Data to Good Lexical ERPs
is_good = is_good_trial(ERPs);
is_lex = strcmpi(ERPs.annot.filledLexTrans, 'lexical');
Data_Tag = is_good & is_lex & ERPs.grouped_fill_annot.is_unique_erp;


ecog = ERPs.ecog;

%% Establish High Level Paramters for Regression, Clustering and Segmentation
if plot_raw | plot_segmentation | plot_regression | plot_clustering | plot_kinetic
    Data_Tag_hs_onset = Data_Tag & strcmpi(ERPs.grouped_fill_annot.grouping_hs,'onset') & ERPs.grouped_fill_annot.is_unique_erp;;    % used for handshape segmentation
    Data_Tag_loc_onset = Data_Tag & strcmpi(ERPs.grouped_fill_annot.grouping_loc,'onset') & ERPs.grouped_fill_annot.is_unique_erp;;  % used for location segmentation
    
    is_ling_onset = get_ling_onset(ERPs);
    Data_Tag_ling_onset = Data_Tag & is_ling_onset; % Main Data_Tag for Analysis of Clustering and Regression

    Regress_Mat_ling_onset = Regress_ECoG(ERPs, Data_Tag_ling_onset);
    sig_chans = get_sig_chans_ii(Regress_Mat_ling_onset);
end


%%
if plot_raw
%% Plot Face and Hand Locations w/5 Handshape;
if plot_face_hand_in5
    comp1 = 'face'; comp2 = 'hand'; eq_param = '5';
    tags = {comp1, comp2, eq_param};

%     is_hs_onset = strcmpi(ERPs.grouped_fill_annot.grouping_hs,'onset');
%     is_loc_onset = strcmpi(ERPs.grouped_fill_annot.grouping_loc,'onset');
%     is_hs_sustained = strcmpi(ERPs.grouped_fill_annot.grouping_hs,'sustained');
    [comp_inds] = compare_onsets_with_eq_param(ERPs, comp1, comp2, eq_param, 'handshape');
    
%     inds1 = strcmpi(ERPs.grouped_fill_annot.loc, comp1) & Data_Tag & strcmpi(ERPs.grouped_fill_annot.handshape,eq_param) & is_hs_onset & is_loc_onset;
%     inds2 = strcmpi(ERPs.grouped_fill_annot.loc, comp2) & Data_Tag & strcmpi(ERPs.grouped_fill_annot.handshape,eq_param) & is_hs_onset & is_loc_onset;
%     comp_inds = [inds1, inds2];
    
    %% Plot Full Grid
    if plot_full_grids
        PlotECogGrid_color_ctrl(ERPs, true, colorlist, ecog(:,:,comp_inds(:,1)), ecog(:,:,comp_inds(:,2)))
    end 
    %% Plot significant channels from regression
    if plot_sig_chans
        Plot_Sig_Chans(ERPs, Data_Tag, sig_chans, comp_inds(:,1), comp_inds(:,2))
    end
    
    %% Plot a select list of electrodes
    if plot_select_electrodes
        plot_electrodes = [110, 133, 147, 201];
        Plot_Select_Electrodes(plot_electrodes, ERPs, ecog, comp_inds, colorlist, tags)
    end
    
end

if plot_fs_comp
    is_fs = strcmpi(ERPs.grouped_fill_annot.loc,'fingerspelling');
    comp1 = 'O'; comp2 = 'S'; eq_param = 'fingerspelling'; 
    tags = {comp1, comp2, eq_param};
    
%    comp_inds = Equalize_Tag_Sizes(Data_Tag & is_fs & strcmpi(ERPs.grouped_fill_annot.handshape,comp1), Data_Tag & is_fs & strcmpi(ERPs.grouped_fill_annot.handshape,comp2));
    comp_inds = [Data_Tag & is_fs & strcmpi(ERPs.grouped_fill_annot.handshape,comp1), Data_Tag & is_fs & strcmpi(ERPs.grouped_fill_annot.handshape,comp2)]; % 
    if plot_full_grids
        PlotECogGrid_color_ctrl(ERPs, true, colorlist, ecog(:,:,comp_inds(:,1)), ecog(:,:,comp_inds(:,2)))
    end
    Plot_Sig_Chans(ERPs, Data_Tag, sig_chans, comp_inds(:,1), comp_inds(:,2))
    
    plot_electrodes = [78, 179, 218 232];
    Plot_Select_Electrodes(plot_electrodes, ERPs, ecog, comp_inds, colorlist, tags)

    comp1 = 'B'; comp2 = 'S';
    % comp_inds = Equalize_Tag_Sizes(Data_Tag & is_fs & strcmpi(ERPs.grouped_fill_annot.handshape,comp1), Data_Tag & is_fs & strcmpi(ERPs.grouped_fill_annot.handshape,comp2));
    comp_inds = [Data_Tag & is_fs & strcmpi(ERPs.grouped_fill_annot.handshape,comp1), Data_Tag & is_fs & strcmpi(ERPs.grouped_fill_annot.handshape,comp2)]; % 
    if plot_full_grids
        PlotECogGrid_color_ctrl(ERPs, true, colorlist, ecog(:,:,comp_inds(:,1)), ecog(:,:,comp_inds(:,2)))
    end
    Plot_Sig_Chans(ERPs, Data_Tag, sig_chans, comp_inds(:,1), comp_inds(:,2))


end
    %comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag_ling_onset,'handshape', strcmpi(ERPs.annot.handshape.
     

if plot_O_5_neut
    comp1 = 'O'; comp2 = 'S'; eq_param = 'neutral';
    %comp1 = 'O'; comp2 = '5'; eq_param = 'neutral';
    tags = {comp1, comp2, eq_param};

    [comp_inds] = compare_onsets_with_eq_param(ERPs, comp1, comp2, eq_param, 'handshape');
    if plot_full_grids
        PlotECogGrid_color_ctrl(ERPs, true, colorlist, ecog(:,:,comp_inds(:,1)), ecog(:,:,comp_inds(:,2)))
    end 
    if plot_select_electrodes
        plot_electrodes = [201, 232, 218, 203];
        plot_electrodes = [201, 232, 218, 223];
        Plot_Select_Electrodes(plot_electrodes, ERPs, ecog, comp_inds, colorlist, tags)
    end
        %% Plot significant channels from regression
    if plot_sig_chans
        Plot_Sig_Chans(ERPs, Data_Tag, sig_chans, comp_inds(:,1), comp_inds(:,2))
    end
    
end

end



%% plot regression
if plot_regression
    sig_thresh = 0.01/256
    
    % PLOT GRID of R^2
    tmp = double(Regress_Mat_ling_onset(:,:,end)<(sig_thresh)); %mask
    if plot_full_grids
        PlotECogGrid_Gen(ERPs,false, Regress_Mat_ling_onset(:,:,end-1), Regress_Mat_ling_onset(:,:,end-1).*(tmp./tmp))
    end
        %
    % PLOT BRAIN MAP
    tmp = double(sig_chans).*max(squeeze(Regress_Mat_ling_onset(:,:,(end-1))),[],2);
    plot_brain_elecs_sandbox_greycbar(tmp)
    title('Maximum R^2 in Significant Electrodes')
end

%% Plot Clustering
%% Need to Modify Cluster Function to accept variable inputs
if plot_clustering
    % HS-Loc
    % CLEAR FRINGE CASES:
    % 
    Data_Tag_cluster = Data_Tag_ling_onset & ~strcmpi(ERPs.grouped_fill_annot.handshape,'i') & ~strcmpi(ERPs.grouped_fill_annot.handshape,'k/p') & ~strcmpi(ERPs.grouped_fill_annot.handshape,'f') &~strcmpi(ERPs.grouped_fill_annot.loc,'chest')  &~strcmpi(ERPs.grouped_fill_annot.loc,'hand') &~strcmpi(ERPs.grouped_fill_annot.handshape,'5>');
    ClusterGestures_Comb_grouped_variable_comp(ERPs, Data_Tag_cluster, 6, sig_chans, 12, 'handshape', 'loc');
    % Loc-IntMov
    ClusterGestures_Comb_grouped_variable_comp(ERPs, Data_Tag_ling_onset, 6, sig_chans, 14, 'loc', 'intMov');
    % Int_Mov-HS
    ClusterGestures_Comb_grouped_variable_comp(ERPs, Data_Tag_ling_onset, 6, sig_chans, 14, 'handshape', 'intMov');
end

%% Plot Segmentation
if plot_segmentation
    Data_Tag_SVM = is_good & is_lex & strcmpi(ERPs.grouped_fill_annot.grouping_hs,'onset') & ERPs.grouped_fill_annot.is_unique_erp;
    % Timings:
    advance = 2;
    Analysis_frame = 10;
    frames = floor((size(ecog,2) - Analysis_frame)/advance); % clippout the last frame
    time_axis_svm = ERPs.time_axis(floor(((1:frames)-1)*advance + Analysis_frame/2));
    
    t_range = [-1800 1800];
    in_range = (time_axis_svm <= t_range(2)) & (time_axis_svm >=t_range(1));
    time_axis_svm = time_axis_svm(in_range);
    
    % Assumes_Existence of 
    if segment_hs
        Data_Tag_SVM = is_good & is_lex & strcmpi(ERPs.grouped_fill_annot.grouping_hs,'onset') & ERPs.grouped_fill_annot.is_unique_erp;
        [mean_accuracy_hs, std_error_hs] = Segment_Handshapes_Slide_KFold(ERPs, Data_Tag_SVM, sig_chans, 50, false, 0, 'handshape');
        [mean_accuracy_hs_rand, std_error_hs_rand] =  Segment_Handshapes_Slide_KFold(ERPs, Data_Tag_SVM, sig_chans, 50, true, 0, 'handshape');
        
        % Trim the Data to fit only the correct t_range
        mean_accuracy_hs = mean_accuracy_hs(in_range); mean_accuracy_hs_rand = mean_accuracy_hs_rand(in_range);
        std_error_hs = std_error_hs(in_range); std_error_hs_rand = std_error_hs_rand(in_range);

        figure; shadedErrorBar(time_axis_svm, mean_accuracy_hs, std_error_hs, {'color', [0 0.2 0.8]},1); hold;
        shadedErrorBar(time_axis_svm, mean_accuracy_hs_rand, std_error_hs_rand,{'color', [0.8 0.2 0]},1);
        xlabel('Time from Onset (ms)')
        ylabel('Classification Accuracy')
        title('Plot of Handshape SVM Classifier Accuracy over Time')
        plot([0 0], [0 0.5],'k','LineWidth',2)
        plot(t_range, [0.125 0.125],'k--','LineWidth',2)
        xlim(t_range)
        ylim([0 0.5])
    end
    
    if segment_loc
        Data_Tag_SVM = is_good & is_lex & strcmpi(ERPs.grouped_fill_annot.grouping_loc,'onset') & ERPs.grouped_fill_annot.is_unique_erp;
        [mean_accuracy_loc, std_error_loc] = Segment_Handshapes_Slide_KFold(ERPs, Data_Tag_SVM, sig_chans, 50, false, 0, 'loc');
        [mean_accuracy_loc_rand, std_error_loc_rand] =  Segment_Handshapes_Slide_KFold(ERPs, Data_Tag_SVM, sig_chans, 50, true, 0, 'loc');
        
          % Trim the Data to fit only the correct t_range
        mean_accuracy_loc = mean_accuracy_loc(in_range); mean_accuracy_loc_rand = mean_accuracy_loc_rand(in_range);
        std_error_loc = std_error_loc(in_range); std_error_loc_rand = std_error_loc_rand(in_range);

        figure; shadedErrorBar(time_axis_svm, mean_accuracy_loc, std_error_loc, [0 0.2 0.8],1); hold;
        shadedErrorBar(time_axis_svm, mean_accuracy_loc_rand, std_error_loc_rand,[0.8 0.2 0],1);
        xlabel('Time from Onset (ms)')
        ylabel('Classification Accuracy')
        title('Plot of Location SVM Classifier Accuracy over Time')
        plot([0 0], [0 0.8],'k','LineWidth',2)
        plot(t_range, [0.25 0.25],'k--','LineWidth',2)
        xlim(t_range)
        ylim([0 0.8])
    end
    
    %% Plot Conufsion Matricies
    if create_confusion_mats
        if segment_hs
            Data_Tag_SVM = is_good & is_lex & strcmpi(ERPs.grouped_fill_annot.grouping_hs,'onset') & ERPs.grouped_fill_annot.is_unique_erp;
            Conf_Mats_HS = Segment_Handshapes_Slide_Cleaned(ERPs, Data_Tag_SVM, sig_chans, 50, false, 0, 'handshape');
            Categories_HS = Get_Categories(ERPs, Data_Tag_SVM, 50, 'handshape');
            plot_confusion_mat_sequence(time_axis_svm, Categories_HS, Conf_Mats_HS) % confusion_mat_1) %, confusion_mat_2)

        end
        if segment_loc
            Data_Tag_SVM = is_good & is_lex & strcmpi(ERPs.grouped_fill_annot.grouping_loc,'onset') & ERPs.grouped_fill_annot.is_unique_erp;
            Conf_Mats_Loc = Segment_Handshapes_Slide_Cleaned(ERPs, Data_Tag_SVM, sig_chans, 50, false, 0, 'loc');
            Categories_Loc = Get_Categories(ERPs, Data_Tag_SVM, 50, 'loc');
            plot_confusion_mat_sequence(time_axis_svm, Categories_Loc, Conf_Mats_Loc) % confusion_mat_1) %, confusion_mat_2)
            a = 1;
        end
    end
end



%%*****************************************



%%


%%


%%******************************************


%% Plot Kinetic Features:
if plot_kinetic
    twin = 170:230; % Relevant Area to measure difference
        regress_grid_null = Regress_ECoG(ERPs, Data_Tag & get_ling_onset(ERPs), true);
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
    is_lex_sig_loc = is_lex & (is_face | is_hand | is_fs | is_neut);

    %% Plot FS vs Lex with Equal Dist
    is_lex_fs = is_lex & (is_fs);
    is_lex_lex = is_lex & (is_face | is_hand | is_neut);
    % equalize HS
    comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag, 'handshape', is_lex_lex, is_lex_fs);

    is_lex_lex_hs = false(size(is_lex_fs));
    is_lex_lex_hs(comp_inds(:,1)) = true;

    is_lex_fs_hs = false(size(is_lex_lex));
    is_lex_fs_hs(comp_inds(:,2)) = true;

    comp_dur = Equalize_Duration_Distribution(ERPs, Data_Tag, is_lex_lex_hs, is_lex_fs_hs);
    if plot_full_grids
        PlotECogGrid_Gen(ERPs, true, ERPs.ecog(:,:,comp_dur(:,1)), ERPs.ecog(:,:,comp_dur(:,2)))
    end

    %% Plot Trans vs Ling
    comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag, 'handshape', is_trans, is_lex_sig_loc);

    if plot_full_grids     
        PlotECogGrid_Gen_Blackout(ERPs, true, ~is_sig_non_loc, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)))
    end
    if plot_tl_diff_brain
        trans_ling_diff = mean(ERPs.ecog(:,twin,comp_inds(:,1)),3)-mean(ERPs.ecog(:,twin,comp_inds(:,2)),3);
        [~,max_ind] = max(abs(trans_ling_diff),[],2);
        for i = 1:256
            peak_diff(i) = trans_ling_diff(i,max_ind(i));
        end
        tmp = double(is_sig_non_loc).*peak_diff';
        plot_brain_elecs_sandbox_greycbar(tmp)
        title({'Maximum Difference in High Gamma Intensity of Linguistic and';'Transitional Gestures in Significant non-Location Electrodes'})
    end

    %% Split Lex & Fingerspelling:
    comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag, 'handshape', is_trans, is_lex & (is_face | is_hand | is_neut), is_lex & is_fs);

    if plot_full_grids
        PlotECogGrid_Gen_Blackout(ERPs, true, ~is_sig_non_loc, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)),  ERPs.ecog(:,:,comp_inds(:,3)))
    end



        %% Plot Trans-Lex Difference on the Brain:
        comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag, 'handshape', is_trans, is_lex_sig_loc);

        twin = 170:230; 
        trans_mean_smth = mean(ERPs.ecog(:,twin,comp_inds(:,1)),3);
        lex_mean_smth = mean(ERPs.ecog(:,twin,comp_inds(:,2)),3);
        trans_lex_diff = trans_mean_smth - lex_mean_smth;

        [~,max_ind] = max(abs(trans_lex_diff),[],2);
        for i = 1:256
            peak_diff(i) = trans_lex_diff(i,max_ind(i));
        end
        tmp = double(is_sig_non_loc).*peak_diff';
        
         plot_brain_elecs_sandbox_greycbar(tmp)
         title({'Maximum Difference in High Gamma Intensity of Linguistic and';'Transitional Gestures in Significant non-Location Electrodes'})

        %% Plot Linguistic_Latency
   %     [~,peak_time] = max(mean(ERPs.ecog(:,twin,is_lex_sig_loc),3),[],2);
    %    [~,peak_time] = min(abs(diff(smooth_grid_time(mean(ERPs.ecog(:,twin,is_lex_sig_loc),3),20),1,2)),[],2);

       % peak_time = (peak_time-30)*10;
        peak_time = (max_ind - 30)*10;
        tmp = double(is_sig_non_loc).*peak_time;
        plot_brain_elecs_sandbox_greycbar(tmp)
        title({'Latency Time of Peak High Gamma for Linguistic Gestures'; 'in Significant non-Location Electrodes'})
        %make_two_type_cbar(min(tmp), max(tmp))

comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag, 'handshape', is_trans, is_lex_lex, is_lex_fs);
if plot_full_grids
    PlotECogGrid_Gen_Blackout(ERPs, true, ~is_sig_non_loc, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)),  ERPs.ecog(:,:,comp_inds(:,3)))
end
    
    

    %% Addenda and Erratera:
    
if plot_contact_diff
    twin = 170:230; % Relevant Area to measure difference
    regress_grid_null = Regress_ECoG(ERPs, Data_Tag & get_ling_onset(ERPs), true);
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

    is_fs = strcmpi(ERPs.grouped_fill_annot.loc,'fingerspelling');

    is_uncontact = strcmpi(ERPs.annot.contact, 'uncontacted');
    is_no_contact = strcmpi(ERPs.annot.contact, 'nocontact');
    is_contact = strcmpi(ERPs.annot.contact,'contact');
    is_contacting = is_contact | is_uncontact;

    % Compare Contact, UnContact, No Contact (include FS)
    comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag_ling_onset, 'handshape', is_contact, is_no_contact,is_uncontact);
    PlotECogGrid_Gen_Blackout(ERPs, true, ~is_sig_non_loc, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)),  ERPs.ecog(:,:,comp_inds(:,3)))
    
    % Compare Contact, UnContact, No Contact (exclude FS)
    comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag_ling_onset, 'handshape', is_contact&~is_fs, is_no_contact&~is_fs,is_uncontact&~is_fs);
    PlotECogGrid_Gen_Blackout(ERPs, true, ~is_sig_non_loc, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)),  ERPs.ecog(:,:,comp_inds(:,3)))


    % Compare Contact, UnContact
    comp_inds_con_uncon = Equalize_Parameter_Distribution_N(ERPs, Data_Tag_ling_onset, 'handshape', is_contact, is_uncontact);
    PlotECogGrid_Gen_Blackout(ERPs, true, ~is_sig_non_loc, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)))

    % Compare Contact, No-Contact
    comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag_ling_onset, 'handshape', is_contact, is_no_contact);
    PlotECogGrid_Gen_Blackout(ERPs, true, ~is_sig_non_loc, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)))
 
    % Compare Contact, No-Contact (fs_excluded)
    comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag_ling_onset, 'handshape', is_contact&~is_fs, is_no_contact&~is_fs);
    PlotECogGrid_Gen_Blackout(ERPs, true, ~is_sig_non_loc, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)))
    
    % Compare UnContact, No-Contact (fs_excluded)
    comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag_ling_onset, 'handshape', is_uncontact&~is_fs, is_no_contact&~is_fs);
    PlotECogGrid_Gen_Blackout(ERPs, true, ~is_sig_non_loc, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)))

    %% find electrodes where the difference between contact and uncontact is insignificant in twin
    sig_ch = find(is_sig_non_loc);
    is_sig_cont_uncont = false(size(is_sig_non_loc));
    for i = 1:length(sig_ch)
        ch_dat_cont = squeeze(ERPs.ecog(sig_ch(i),twin, comp_inds_con_uncon(:,1)));
        ch_dat_uncont = squeeze(ERPs.ecog(sig_ch(i),twin, comp_inds_con_uncon(:,2)));
        sig_pts = 0;
        for j = 1:length(twin)
            ks_sig = kstest2(ch_dat_cont(j,:),ch_dat_uncont(j,:));
            if ks_sig
                sig_pts = sig_pts+1;
            end
        end
        
        is_sig_cont_uncont(sig_ch(i)) = (sig_pts>=10);
    end
    
    cont_uncont_same = is_sig_non_loc & ~is_sig_cont_uncont;
    
    comp_inds = Equalize_Parameter_Distribution_N(ERPs, Data_Tag_ling_onset, 'handshape', is_contacting & ~is_fs, is_no_contact & ~is_fs);
    PlotECogGrid_Gen_Blackout(ERPs, true, ~cont_uncont_same, ERPs.ecog(:,:,comp_inds(:,1)), ERPs.ecog(:,:,comp_inds(:,2)))
    
    mean_contact_diff = zeros(length(cont_uncont_same),1);
    for i = 1:length(cont_uncont_same)
        if cont_uncont_same(i)
            ch_dat_cont = squeeze(ERPs.ecog(i,twin, comp_inds(:,1)));
            ch_dat_nocont = squeeze(ERPs.ecog(i,twin, comp_inds(:,2)));
            mean_contact_diff(i) = max(mean(ch_dat_cont,2) - mean(ch_dat_nocont,2));
        end
    end
     plot_brain_elecs_sandbox_greycbar(mean_contact_diff)
     title({'Mean Difference in High Gamma Intensity of Contacting and';'Non-Contacting Gestures in Significant non-Location, non-uncontacting Electrodes'})

        
end

a = 1;
%         
%         
%         figure;
%         for i = 1:length(plot_electrodes)
%             % Create a subplot in the squarest grid
%             subplot(ceil(sqrt(length(plot_electrodes))),ceil(sqrt(length(plot_electrodes))),i)
%             data = squeeze(ecog(plot_electrodes(i),:,:));
%             time_axis = ERPs.time_axis;
%             for j = 1:size(comp_inds,2)
%                 shadedErrorBar(time_axis,mean(data(:,comp_inds(:,j)),2), nansem(data(:,comp_inds(:,j)),2),colorlist(j,:),1);
%                 hold on;
%             end
%             axis tight;
% %             set(gca,'YLim',[minbound maxbound]...
% %                 ,'XTickLabel',[],'YTickLabel',[],'YTick',[minbound maxbound]);
%             line([0 0],get(gca,'YLim'),'Color','k');
%             line(get(gca,'XLim'),[0 0],'Color','k');
%             xlabel('Time From Onset (ms)')
%             ylabel('High Gamma Intensity')
%             title(['Plot of High Gamma in Ch ', num2str(plot_electrodes(i)), ' for ', comp1, ' and ', comp2, ' Onset ERPs with equalized', eq_param, 'position'])
%         end      
%     end
% end
% 
% 
% 
% 
% % is_comp1 = strcmpi(ERPs.grouped_fill_annot.loc,comp1)...
% %     & strcmpi(ERPs.grouped_fill_annot.grouping_loc,'onset')...
% %         & ERPs.grouped_fill_annot.is_unique_erp & Data_Tag;
% %     
% is_comp2 = strcmpi(ERPs.grouped_fill_annot.loc,comp2)...
%     & strcmpi(ERPs.grouped_fill_annot.grouping_loc,'onset')...
%         & ERPs.grouped_fill_annot.is_unique_erp & Data_Tag;
%    
% is_eq_onset = strcmpi(ERPs.grouped_fill_annot,eq_param)...
%     & strcmpi(ERPs.grouped_fill_annot.grouping_hs,'onset');
% is_eq_ongoing = strcmpi(ERPs.grouped_fill_annot,eq_param)...
%     & strcmpi(ERPs.grouped_fill_annot.grouping_hs,'sustained');
% 
% Comp12_on_Eq = Equalize_Tag_Sizes(is_comp1 & is_eq_onset, is_comp2 & is_eq_onset)...
%     |  Equalize_Tag_Sizes(is_comp1 & is_eq_ongoing, is_comp2 & is_eq_ongoing);
% 
% PlotECogGrid_Gen(ERPs, true, ERPs.ecog(:,:,Comp12_on_Eq(:,1)), ERPs.ecog(:,:,Comp12_on_Eq(:,2)))
