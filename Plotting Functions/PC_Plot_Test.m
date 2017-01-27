function [] = PC_Plot_SandBox(ERPs, sig_chans)
Data_Tag = strcmpi(ERPs.annot.filledLexTrans, 'lexical');
Data_Tag = Data_Tag & is_good_trial(ERPs);

ecog = ERPs.ecog(sig_chans, :, Data_Tag);
dims = size(ecog);
ecog_reshape = reshape(ecog, [dims(1), dims(2)*dims(3)])'; % Transpose Necessary for PCA

%%%
% ecog_reshape = squeeze(mean(ecog,3))';
% ecog_reshape = squeeze(mean(ecog,2))';

[pc_mat, ~,~,~,exp] = pca(ecog_reshape);

ecog_reshape = reshape(ecog, [dims(1), dims(2)*dims(3)])';

%% Project Data Into PCs
num_comps = 2;
ecog_reshape = gsubtract(ecog_reshape, mean(ecog_reshape,1)); % Center ECoG_Reshape;
ecog_reshape_pca = (ecog_reshape*pc_mat)'; % Transpose back to rows by cols


ecog_pca = reshape(ecog_reshape_pca, dims);
ecog_pca = ecog_pca(1:num_comps,:,:); % Take first components of ecog_pca data
%ecog_pca = ecog_pca(2:3,:,:);



%% Compare Fingerspleeing and Lexical in PC space
class1 = (strcmpi(ERPs.annot.filledLexTrans, 'lexical') & strcmpi(ERPs.annot.respType,'fs'));
class2 = (strcmpi(ERPs.annot.filledLexTrans, 'lexical') & ~strcmpi(ERPs.annot.respType,'fs'));

class1 = class1(Data_Tag);
class2 = class2(Data_Tag);


%% Equalize For Handshape
eq_inds = Equalize_Parameter_Distribution(ERPs, Data_Tag, class1, class2);
ind1 = eq_inds(:,1);
ind2 = eq_inds(:,2);

data1_pca = ecog_pca(:,:,ind1);
data2_pca = ecog_pca(:,:,ind2);
%% Plot the Data

mean_1 = mean(data1_pca,3);
mean_2 = mean(data2_pca,3);

std_err_1 = nansem(data1_pca,3);
std_err_2 = nansem(data2_pca,3);



% figure;
% plot(mean_1(1,:),mean_1(2,:),'b', mean_2(1,:), mean_2(2,:),'r')
% xlabel('pc 1')
% ylabel('pc 2')
% hold
% plot([mean_1(1,1), mean_2(1,1)] , [mean_1(2,1), mean_2(2,1)],'yo', 'LineWidth', 2.5)
% plot([mean_1(1,201), mean_2(1,201)], [mean_1(2,201), mean_2(2, 201)], 'go', 'LineWidth', 3)
% legend('Fingerspelling', 'Lexical', 'Start Point', 'Event Onset')
% hold

 figure;
 plot([mean_1(1,1)], [mean_1(2,1)],'b', [mean_2(1,1)], [mean_2(2,1)],'r')
 xlabel('pc 1')
 ylabel('pc 2')
 hold
 plot([mean_1(1,1), mean_2(1,1)] , [mean_1(2,1), mean_2(2,1)],'yo', 'LineWidth', 3)
plot([mean_1(1,201), mean_2(1,201)], [mean_1(2,201), mean_2(2, 201)], 'go', 'LineWidth', 3)
legend('Fingerspelling', 'Lexical', 'Start Point', 'Event Onset')

%% Assemble Error Bar Regions
%figure; hold
t = (0:1/2560:1)'*2*pi;
x1 = [];
y1 = [];
x2 = [];
y2= [];
for i = 1:length(mean_1(1,:))
    plot([mean_1(1,i)],[mean_1(2,i)],'b.', [mean_2(1,i)], [mean_2(2,i)],'r.')
    
    x1 = [(mean_1(1,i)+std_err_1(1,i)*cos(t))];
    y1 = [mean_1(2,i)+std_err_1(2,i)*sin(t)];
    
    x2 = [mean_2(1,i)+std_err_2(1,i)*cos(t)];
    y2 = [mean_2(2,i)+std_err_2(2,i)*sin(t)];


    fill(x1,y1,'b', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)
    fill(x2,y2,'r', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)
    drawnow
end
plot([mean_1(1,1), mean_2(1,1)] , [mean_1(2,1), mean_2(2,1)],'yo', 'LineWidth', 3)
plot([mean_1(1,201), mean_2(1,201)], [mean_1(2,201), mean_2(2, 201)], 'go', 'LineWidth', 3)
title({'Parametric Plot of Mean Value of Lexical and Fingerspelling Gestures over Time on the'; 'First 2 Principle Components of all Linguistic ERPs in Significant Channels'})

end