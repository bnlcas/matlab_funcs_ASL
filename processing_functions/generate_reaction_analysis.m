function [] = generate_reaction_analysis(ERPs)
clear_unreason = @(x,y) x(~y);

%% Timing Flags
hist_kin_onset = false;  % Histogram reaction times
hist_yesno = false; % Histogram reaction times to yes/no (not intereting (n=17))
hist_lex_react = false;
hist_dup_react = false;
plot_brain_map = true;

%% ERPs plotting Flags:
plot_stim_locked = false; % (false plots Lex Movement Locked)

is_stim_onset = strcmpi(ERPs.annot.stimOnset,'S');

start_times = ERPs.annot.start_samp/10; % in ms
stop_times = ERPs.annot.end_samp/10;
is_ling = strcmpi(ERPs.annot.lexTrans,'lexical');
is_mov = strcmpi(ERPs.annot.movOnset,'M');
[is_ling_first, is_stim_matched] = Get_First_Subsequent_Tag(is_stim_onset, is_ling);
% Clear bad trials on stim or ling erps:
is_ling_first = is_ling_first & is_good_trial_block(ERPs);
is_stim_matched = is_stim_matched & is_good_trial_block(ERPs);

lexical_reaction_time = start_times(is_ling_first) - start_times(is_stim_matched);

freqs = ERPs.annot.frequency(is_stim_matched);
unreasonable_timings = lexical_reaction_time < 500 | lexical_reaction_time >3000;
freqs(unreasonable_timings) = [];
lexical_reaction_time(unreasonable_timings) = [];

handshapes = clear_unreason(ERPs.grouped_fill_annot.handshape(is_ling_first),unreasonable_timings);
% is psuedo
is_non_sign = (freqs == 0);
is_sign = ~is_non_sign;
% is duplicate
is_dup = clear_unreason(strcmpi(ERPs.annot.respType(is_ling_first),'dup'),unreasonable_timings);


%% Plot Reaction Histgram for Real/Pseudo Signs
if hist_lex_react
    figure; histogram(lexical_reaction_time(is_sign),25)
    hold on;
    histogram(lexical_reaction_time(is_non_sign),25)

    legend(['real sign (mean = ', num2str(mean(lexical_reaction_time(is_sign))) ,')'],...
        ['pseudo-sign (mean = ', num2str(mean(lexical_reaction_time(is_non_sign))),')'])
    ylabel('Number of Instances')
    xlabel('Reaction Time (ms)')
    title('Distribution of Linguistic Reaction Time to Stimuli')
end


%% Plot Reaction Hist for Real/Pseudo Duplication:
if hist_dup_react
    figure; histogram(lexical_reaction_time(is_sign & is_dup),25)
    hold on;
    histogram(lexical_reaction_time(is_non_sign & is_dup),25)
    legend(['real sign (mean = ', num2str(mean(lexical_reaction_time(is_sign & is_dup))),')'],...
        ['pseudo-sign (mean = ', num2str(mean(lexical_reaction_time(is_non_sign & is_dup))),')'])
    title('Distribution of Reaction Time for Duplications')
end

%% Plot Reaction Hist for YES/NO distinction
if hist_yesno
    is_yesno = clear_unreason((strcmpi(ERPs.annot.respType(is_ling_first),'yes') | strcmpi(ERPs.annot.respType(is_ling_first),'no')),unreasonable_timings);
    figure; histogram(lexical_reaction_time(is_sign & is_yesno),25)
    hold on;
    histogram(lexical_reaction_time(is_non_sign & is_yesno),25)
    legend(['real sign (mean = ', num2str(mean(lexical_reaction_time(is_sign & is_yesno))),')'],...
        ['pseudo-sign (mean = ', num2str(mean(lexical_reaction_time(is_non_sign & is_yesno))),')'])
    title('Distribution of Reaction Time for YES/NO Responses')
end



%% Compare ERPs to Duplicates restricted to the same set of duplicate words
if plot_stim_locked
    ecog = ERPs.ecog(:,:,is_stim_matched);
else
    ecog = ERPs.ecog(:,:,is_ling_first);
end
ecog(:,:,unreasonable_timings) = [];
%dups = clear_unreason(strcmpi(ERPs.annot.respType(is_ling_first),'dup'));

%PlotECogGrid_Gen(ERPs, true, ecog(:,:,is_dup & is_sign), ecog(:,:,is_dup & is_non_sign))


%% Compare Equal Handshapes
handshapes = clear_unreason(ERPs.grouped_fill_annot.handshape(is_ling_first), unreasonable_timings);
%[eq_inds] = Equalize_Parameter_Distribution(handshapes, true(size(handshapes)), is_dup & is_sign, is_dup & is_non_sign);
load('/Users/changlab/Documents/changrepo/matlab/analysis/ASL/Graphics/January Graphics/PseudoSign Analysis/Figures/matlab.mat');
eq_inds = eq_inds_ii;

%PlotECogGrid_Gen(ERPs, true, ecog(:,:,eq_inds(:,1)), ecog(:,:,eq_inds(:,2)))

%% Difference Between Real - Psuedo Signs
real_pseudo_diff = ecog(:,:,eq_inds(:,1)) - ecog(:,:,eq_inds(:,2));
real_pseudo_diff_mean = mean(real_pseudo_diff,3);
%PlotECogGrid_Gen(ERPs, true, real_pseudo_diff)

pmat = ttest_ecog_diff(ecog(:,:,eq_inds(:,1)), ecog(:,:,eq_inds(:,2)));
pthresh = 0.05;
consequtive_pts = 1;
sig_chans_diff = get_sig_chans_pvals(pmat(:,200:400), pthresh/200, consequtive_pts);

real_pseudo_diff_sig = real_pseudo_diff_mean.*double(pmat < pthresh);

[~, max_ind] = max(abs(real_pseudo_diff_sig),[],2); % largest difference (abs) must be resorted to get sign
 for i = 1:256
    peak_diff(i) = real_pseudo_diff_sig(i,max_ind(i));
 end
 
if plot_brain_map
    tmp = peak_diff'.*double(sig_chans_diff);
    plot_brain_elecs_sandbox_greycbar(tmp)
    title({'Maximum Significant Difference in High Gamma during the';'Production of Lexical - Pseudo Signs'})
end

mask = repmat(double(sig_chans_diff),1,size(ecog,2));
dat = real_pseudo_diff_mean.*mask;
%Frames_Psuedo = plot_brain_elecs_sandbox_greycbar_movie(ERPs, dat);
%% Project into PCA
ecog_dups = ecog(:,:,is_dup);
dat = squeeze(mean(ecog_dups(sig_chans_diff,200:300,:),2))';
pca_dat = dat - repmat(mean(dat,1),size(dat,1),1);
[pca_mat,~,~,~,exp] = pca(pca_dat);
pca_proj = pca_dat*pca_mat(:,1:3);

figure; plot(pca_proj(is_sign(is_dup),1), pca_proj(is_sign(is_dup),2),'b.', pca_proj(is_non_sign(is_dup),1),pca_proj(is_non_sign(is_dup),2),'r.','MarkerSize',12)
xlabel('PC 1'); ylabel('PC 2')
legend('Lexical Signs', 'Pseudo Signs')
title(['Plot of Real and Pseudo Signs in PC Space (' num2str(sum(exp(1:2)),2) '% Explained)'])
hold on;

figure; plot3(pca_proj(is_sign(is_dup),1), pca_proj(is_sign(is_dup),2),pca_proj(is_sign(is_dup),3),'b.', pca_proj(is_non_sign(is_dup),1),pca_proj(is_non_sign(is_dup),2),pca_proj(is_non_sign(is_dup),3),'r.')
xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3')
legend('Lexical Signs', 'Pseudo Signs')
title('Plot of Real and Pseudo Signs in PC Space (32 % Explained)')
hold on;

%% Pair Signs/Pseudo
proj_sign = pca_proj(is_sign(is_dup),:);
proj_nonsign = pca_proj(is_non_sign(is_dup),:);
wordnames = clear_unreason(ERPs.annot.word_name(is_stim_matched),unreasonable_timings);
words_dupsign = wordnames(is_sign(is_dup));
words_dupnonsign = wordnames(is_non_sign(is_dup));
figure; hold on;
%% Plot connection Vectors
for i = 1:length(words_dupnonsign)
    ind_sign = find(strcmpi(words_dupsign, words_dupnonsign(i)),1);
    if ~isempty(ind_sign)
        quiver(0, 0, (proj_nonsign(i,1) - proj_sign(ind_sign,1)), (proj_nonsign(i,2) - proj_sign(ind_sign,2)),'k')
        centroid_vec(i,1) = (proj_nonsign(i,1) - proj_sign(ind_sign,1));
        centroid_vec(i,2) = (proj_nonsign(i,2) - proj_sign(ind_sign,2));
    end
end
quiver(0, 0, mean(centroid_vec(:,1)), mean(centroid_vec(:,2)),'r','LineWidth',2)
xlabel('PC 1'); ylabel('PC 2')
title('Vector Plot of the Trajectory from Real to Pseudo Sign in PC Space')



%% Plot Connector Lines
for i = 1:length(words_dupnonsign)
    ind_sign = find(strcmpi(words_dupsign, words_dupnonsign(i)),1);
    if ~isempty(ind_sign)
        %plot3([proj_nonsign(i,1), proj_sign(ind_sign,1)], [proj_nonsign(i,2), proj_sign(ind_sign,2)], [proj_nonsign(i,3), proj_sign(ind_sign,3)],'k-.')
        plot([proj_nonsign(i,1), proj_sign(ind_sign,1)], [proj_nonsign(i,2), proj_sign(ind_sign,2)], 'k-.')
    end
end