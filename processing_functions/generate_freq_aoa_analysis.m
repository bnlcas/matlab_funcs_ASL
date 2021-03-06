function generate_freq_aoa_analysis(ERPs, Regress_Freq_Aoa_Stim, Regress_Freq_Aoa_Lex)
%% This function generates the analysis of the frequency and aoa regression
% Follow this Sequence:
% Data_Tag_Ling = get_ling_onset(ERPs) & is_good_trial(ERPs);
% this function needs to take the Regression matrcies generated by the
% function Regress_ECoG
% with the approprate changes made to the regression parameters to included
% in the function Assemble_predictor_mat_gen
% this would be challenging to automate, because there are a lot of options
% for features available in Assemble_predictor_mat_gen...
%


%% After adjusting parameters in Assemble_predictor_mat_gen:
% 
% % Data_Tag_lex = Data_Tag & strcmpi(ERPs.annot.respType,'dup') & strcmpi(ERPs.annot.lexTrans,'lexical') & (ERPs.annot.frequency ~= 0);
% % Data_Tag_Stim = false(length(Data_Tag),1);
% % stim_on_inds = find(strcmpi(ERPs.annot.stimOnset, 'S'));
% % lex_on_inds = find(Data_Tag_lex);
% % for i = 1:length(lex_on_inds)
% %   stim_inds = find(stim_on_inds < lex_on_inds(i));
% %   Data_Tag_Stim(stim_inds(end)) = true;
% % end
% % Data_Tag_Stim = Data_Tag_Stim & is_good_trial(ERPs);
% % 
% % Regress_Freq_Aoa_Stim = Regress_ECoG(ERPs, Data_Tag_Stim);
% % Regress_Freq_Aoa_Lex = Regress_ECoG(ERPs, Data_Tag_lex);
% %
% % generate_freq_aoa_analysis(ERPs, Regress_Freq_Aoa_Stim, Regress_Freq_Aoa_Lex)

%% Code
good_inds = find(is_good_trial(ERPs));
stim_inds = find(strcmpi(ERPs.annot.stimOnset,'S'));
lex_inds = find(strcmpi(ERPs.annot.lexTrans,'lexical'));
for i = 1:length(stim_inds)
lex_onset_inds(i) = lex_inds(find(lex_inds>stim_inds(i),1));
end
stim_inds = intersect(stim_inds, good_inds);
lex_onset_inds = intersect(lex_onset_inds, good_inds);
Data_Tag_Ling_Onset = false(size(ERPs.annot.stimOnset));
Data_Tag_Stim_Onset = false(size(ERPs.annot.stimOnset));
Data_Tag_Ling_Onset(lex_onset_inds) = true;
Data_Tag_Stim_Onset(stim_inds) = true;

% Regress_XXX_stim/lex = Regress_ECoG(ERPs, Data_Tag_Stim/Lex); 
% THEN ADJUST PARAMETERS IN Assemble_Predictor_Mat ... 

% % Regress_Freq_Aoa_Stim = Regress_ECoG(ERPs, Data_Tag_Stim_Onset);
% % Regress_Freq_Aoa_Lex = Regress_ECoG(ERPs, Data_Tag_Ling_Onset);
% %
% % generate_freq_aoa_analysis(ERPs, Regress_Freq_Aoa_Stim, Regress_Freq_Aoa_Lex)



%% Get Significant Channels
consequtive_sig_thresh = 5;
p_thresh = 0.05/400;
sig_chans_lex = get_sig_chans_ii(Regress_Freq_Aoa_Lex, p_thresh, consequtive_sig_thresh);
sig_chans_stim = get_sig_chans_ii(Regress_Freq_Aoa_Stim, p_thresh, consequtive_sig_thresh);



%% Generate lists of relevant timecourses:

freq_betas_stim = squeeze(Regress_Freq_Aoa_Stim(sig_chans_lex,:,2));
%aoa_betas_stim = squeeze(Regress_Freq_Aoa_Stim(sig_chans_lex,:,3));

freq_betas_lex = squeeze(Regress_Freq_Aoa_Lex(sig_chans_lex,:,2));
%aoa_betas_lex = squeeze(Regress_Freq_Aoa_Lex(sig_chans_lex,:,3));

%% Plot Average Beta Timecourses:
% Plot Stim AoA and Freqquency
figure; subplot(1,2,1);
shadedErrorBar(ERPs.time_axis,mean(freq_betas_stim,1),std(freq_betas_stim,[],1)/sqrt(size(freq_betas_stim,1)),{'color', [0 0 0.8]},1);
hold on;
plot([0 0], get(gca, 'YLim'), 'k', 'LineWidth', 2)
plot(get(gca, 'XLim'), [0 0], 'k', 'LineWidth', 2)
%shadedErrorBar(ERPs.time_axis,mean(aoa_betas_stim,1),std(aoa_betas_stim,[],1)/sqrt(size(aoa_betas_stim,1)),[0.8 0 0],1);

xlabel('Time from stimulus onset (ms)')
ylabel('Mean Beta Weight')
%title({'Mean Stimulus Locked Betas for Frequency (blue)';'and AoA (red) in Significant Channels'})
title({'Mean Stimulus Locked Betas for Frequency in Significant Channels'})


% Plot Lex Aoa and Frequency
subplot(1,2,2);
shadedErrorBar(ERPs.time_axis,mean(freq_betas_lex,1),std(freq_betas_lex,[],1)/sqrt(size(freq_betas_lex,1)),{'color', [0 0 0.8]},1);
hold on;
plot([0 0], get(gca, 'YLim'), 'k', 'LineWidth', 2)
plot(get(gca, 'XLim'), [0 0], 'k', 'LineWidth', 2)

%shadedErrorBar(ERPs.time_axis,mean(aoa_betas_lex,1),std(aoa_betas_lex,[],1)/sqrt(size(aoa_betas_lex,1)),[0.8 0 0],1);

xlabel('Time from movement onset (ms)')
ylabel('Mean Beta Weight')
%title({'Mean Movement Locked Betas for Frequency (blue)';'and AoA (red) in Significant Channels'})
title({'Mean Movement Locked Betas for Frequency in Significant Channels'})


%% Plot R^2 timecourse
r2_lex = squeeze(Regress_Freq_Aoa_Lex(sig_chans_lex,:,end-1));
r2_stim = squeeze(Regress_Freq_Aoa_Stim(sig_chans_lex,:,end-1));

figure;
subplot(1,2,1)
shadedErrorBar(ERPs.time_axis, mean(r2_stim,1),std(r2_stim,[],1)/sqrt(size(r2_stim,1)),'k',1)
hold on;
plot([0 0], get(gca, 'YLim'), 'k', 'LineWidth', 2)
xlabel('Time from stimulus onset (ms)')
ylabel('Mean R^2')
title({'Timecourse of Mean R^2 value in Significant Channels';'For Regression to Stimulus Locked ERPs'})


subplot(1,2,2)
shadedErrorBar(ERPs.time_axis, mean(r2_lex,1),std(r2_lex,[],1)/sqrt(size(r2_lex,1)),'k',1);
hold on;
plot([0 0], get(gca, 'YLim'), 'k', 'LineWidth', 2)
xlabel('Time from movement onset (ms)')
ylabel('Mean R^2')
title({'Timecourse of Mean R^2 value in Significant Channels';'For Regression to Movement Locked ERPs'})


%% Plot Max R2
plot_brain_elecs_sandbox_greycbar(double(sig_chans_stim).*max(squeeze(Regress_Freq_Aoa_Stim(:,:,end-1)),[],2))
title({'Maximum r^2 of Regression of Stimulus Locked ERPs';'to Frequency and Age of Acquisition'})
title({'Maximum r^2 of Regression of Stimulus Locked ERPs to Frequency'})
plot_brain_elecs_sandbox_greycbar(double(sig_chans_lex).*max(squeeze(Regress_Freq_Aoa_Lex(:,:,end-1)),[],2))
title({'Maximum r^2 of Regression of Movement Locked ERPs';'to Frequency and Age of Acquisition'})
title({'Maximum r^2 of Regression of Movement Locked ERPs Frequency'})




%% Plot Average Beta Timecourses:
% Plot Stim AoA and Freqquency
figure; subplot(1,2,1);
shadedErrorBar(ERPs.time_axis,mean(freq_betas_stim,1),std(freq_betas_stim,[],1)/sqrt(size(freq_betas_stim,1)),{'color', [0 0 0.8]},1);
hold on;
plot([0 0], get(gca, 'YLim'), 'k', 'LineWidth', 2)
plot(get(gca, 'XLim'), [0 0], 'k', 'LineWidth', 2)
%shadedErrorBar(ERPs.time_axis,mean(aoa_betas_stim,1),std(aoa_betas_stim,[],1)/sqrt(size(aoa_betas_stim,1)),[0.8 0 0],1);

xlabel('Time from stimulus onset (ms)')
ylabel('Mean Beta Weight')
%title({'Mean Stimulus Locked Betas for Frequency (blue)';'and AoA (red) in Significant Channels'})
title({'Mean Stimulus Locked Betas for Frequency in Significant Channels'})


% Plot Lex Aoa and Frequency
subplot(1,2,2);
shadedErrorBar(ERPs.time_axis,mean(freq_betas_lex,1),std(freq_betas_lex,[],1)/sqrt(size(freq_betas_lex,1)),{'color', [0 0 0.8]},1);
hold on;
plot([0 0], get(gca, 'YLim'), 'k', 'LineWidth', 2)
plot(get(gca, 'XLim'), [0 0], 'k', 'LineWidth', 2)

%shadedErrorBar(ERPs.time_axis,mean(aoa_betas_lex,1),std(aoa_betas_lex,[],1)/sqrt(size(aoa_betas_lex,1)),[0.8 0 0],1);

xlabel('Time from movement onset (ms)')
ylabel('Mean Beta Weight')
%title({'Mean Movement Locked Betas for Frequency (blue)';'and AoA (red) in Significant Channels'})
title({'Mean Movement Locked Betas for Frequency in Significant Channels'})





%% Afterward: if just frequency, plot r values:
if size(Regress_Freq_Aoa_Lex,3)==4
    r_val_lex = sqrt(r2_lex).*sign(freq_betas_stim);
    r_val_stim = sqrt(r2_stim).*sign(freq_betas_lex);
    figure; 
    subplot(1,2,1); hold on;
        shadedErrorBar(ERPs.time_axis, mean(r_val_stim,1),std(r_val_stim,[],1)/sqrt(size(r_val_stim,1)),'k',1)
        plot([0 0], get(gca, 'YLim'), 'k', 'LineWidth', 2)
        plot(get(gca,'XLim'), [0 0], 'k', 'LineWidth',2)
        xlabel('Time from stimulus onset (ms)')
        ylabel('Correlation Coefficient')
        title({'Timecourse of Mean r value in Significant Channels';'Between Stimulus Locked High Gamma and Lexical Frequency'})
    subplot(1,2,2); hold on;
        shadedErrorBar(ERPs.time_axis, mean(r_val_lex,1),std(r_val_lex,[],1)/sqrt(size(r_val_lex,1)),'k',1)
        plot([0 0], get(gca, 'YLim'), 'k', 'LineWidth', 2)
        plot(get(gca,'XLim'), [0 0], 'k', 'LineWidth',2)
        xlabel('Time from movement onset (ms)')
        ylabel('Correlation Coefficient')
        title({'Timecourse of Mean r value in Significant Channels';'Between Movement Locked High Gamma and Lexical Frequency'})
end
    
