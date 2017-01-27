function [] = PC_Plot_Peaks(ERPs, sig_chans)
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
num_comps = 65;
ecog_reshape = gsubtract(ecog_reshape, mean(ecog_reshape,1)); % Center ECoG_Reshape;
ecog_reshape_pca = (ecog_reshape*pc_mat)'; % Transpose back to rows by cols


ecog_pca = reshape(ecog_reshape_pca, dims);
ecog_pca = ecog_pca(1:num_comps,:,:); % Take first components of ecog_pca data

data_pca = squeeze(mean(ecog_pca(:,191:211,:),2));

%% Segment L Fingerspelling and L Handshape
class1 = strcmpi(ERPs.annot.respType,'fs');% & strcmpi(ERPs.annot.handshape,'O');
class2 = ~strcmpi(ERPs.annot.respType,'fs');% & strcmpi(ERPs.annot.handshape,'O');

class1 = class1(Data_Tag);
class2 = class2(Data_Tag);

data1 = data_pca(:,class1);
data2 = data_pca(:,class2);

figure;
for i = 1:64
    subplot(8,8,i)
    scatter(data1(i,:),data1(1,:)*0,'MarkerEdgeColor',[0 0 0.8])
        hold;
    scatter(mean(data1(i,:)),0,150,'d', 'MarkerFaceColor', [0 0 0.8], 'MarkerEdgeColor', [0 0 0.8])

    scatter(data2(i,:),data2(1,:)*0,'MarkerEdgeColor',[0.8 0 0])
    scatter(mean(data2(i,:)),0,150,'d', 'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', [0.8 0 0])
    xlabel(['PC ', num2str(i)]);
%    legend('Fingerspelling', 'Mean Fingerspelling', 'Lexical', 'Mean Lexical')
end


%% Plot Sub Categories of FS vs LEX
    
sub1 = strcmpi(ERPs.annot.handshape,'B');
sub2 = strcmpi(ERPs.annot.handshape,'C');
sub3 = strcmpi(ERPs.annot.handshape,'L');
sub4 = strcmpi(ERPs.annot.handshape,'O');
sub5 = strcmpi(ERPs.annot.handshape,'S');
sub6 = strcmpi(ERPs.annot.handshape,'Y');

figure; hold
num_subs = 6;
for i = 1:num_subs
    data_i = eval(['class1 & sub',num2str(i),'(Data_Tag)']);
    plot_data = data_pca(:,data_i);
    scatter(plot_data(1,:),plot_data(3,:),'MarkerEdgeColor',[0 1*i/num_subs 0.5])
end
for i = 1:num_subs
    data_i = eval(['class2 & sub',num2str(i),'(Data_Tag)']);
    plot_data = data_pca(:,data_i);
    scatter(plot_data(1,:),plot_data(3,:),'MarkerEdgeColor',[1 0.5*i/num_subs 0])
end

figure;
scatter(data1(1,:),data1(3,:),'MarkerEdgeColor',[0 0 0.8])


% class1 = strcmpi(ERPs.annot.respType,'fs') & strcmpi(ERPs.annot.handshape,'O');
% class2 = ~strcmpi(ERPs.annot.respType,'fs') & strcmpi(ERPs.annot.handshape,'O');
% 
handshapes = ERPs.annot.handshape;
handshapes = strrep(handshapes,' ', '');
tag1= strcmpi(ERPs.annot.respType,'fs');
tag2 = ~strcmpi(ERPs.annot.respType,'fs');

categories_fs = unique(handshapes(tag1));
categories_nfs = unique(handshapes(tag2));
categories_fs = categories_fs(~strcmpi(categories_fs,''));
categories_nfs = categories_nfs(~strcmpi(categories_nfs,''));

cat_count_fs = get_category_size(handshapes(tag1));
cat_count_nfs = get_category_size(handshapes(tag2));






%% define SubClasses of Data to display
class1 = strcmpi(ERPs.annot.respType,'fs') & strcmpi(ERPs.annot.handshape,'B');
class2 = strcmpi(ERPs.grouped_annot.loc,'neutral') & strcmpi(ERPs.grouped_annot.handshape,'B');
class3 = strcmpi(ERPs.grouped_annot.loc,'hand') & strcmpi(ERPs.grouped_annot.handshape,'B');
class4 = strcmpi(ERPs.grouped_annot.loc,'hand') & strcmpi(ERPs.grouped_annot.handshape,'5');
class5 = strcmpi(ERPs.grouped_annot.loc,'neutral') & strcmpi(ERPs.grouped_annot.handshape,'5');

varname = 'class'
num_classes = 5;

%% Clear for Data Tag
for i = 1:num_classes
    class_name = strcat(varname, num2str(i));
    eval(strcat(class_name, '=', class_name,'(Data_Tag);'));
end

%% Get Means
for i = 1:num_classes
    class_name = strcat(varname, num2str(i));
    mean_name = strcat('mean', num2str(i));
    eval(strcat(mean_name, ' = mean(data_pca(:,', class_name,'),2);'));
end
%% Get STDevs
for i = 1:num_classes
    class_name = strcat(varname, num2str(i));
    std_name = strcat('stdev', num2str(i));
    eval(strcat(std_name, ' = std(data_pca(:,', class_name,'),[],2);'));
end

%% Plot the Data
figure;
x1 = data_pca(1,class1);
y1 = data_pca(2,class1);
scatter(x1,y1,'MarkerFaceColor',[0 0 0.7])
hold
x2 = data_pca(1,class2);
y2 = data_pca(2,class2);
scatter(x2,y2,'MarkerFaceColor',[0.2 0.2 0.2])

x3 = data_pca(1,class3);
y3 = data_pca(2,class3);
scatter(x3,y3,'MarkerFaceColor',[0.4 0.4 0])

x4 = data_pca(1,class4);
y4 = data_pca(2,class4);
scatter(x4,y4,'MarkerFaceColor',[0.6 0.2 0])

x5 = data_pca(1,class5);
y5 = data_pca(2,class5);
scatter(x5,y5,'MarkerFaceColor',[0.6 0.4 0.4])

legend('Fingerspelling ''B''', 'Neutral ''B''', 'Hand ''B''', 'Hand ''5''', 'Neutral ''5''')

a=1
%% Assemble Error Bar Regions
% %figure; hold
% t = (0:1/2560:1)'*2*pi;
% x1 = [];
% y1 = [];
% x2 = [];
% y2= [];
% for i = 1:length(mean_1(1,:))
%     plot([mean_1(1,i)],[mean_1(2,i)],'b.', [mean_2(1,i)], [mean_2(2,i)],'r.')
%     
%     x1 = [(mean_1(1,i)+std_err_1(1,i)*cos(t))];
%     y1 = [mean_1(2,i)+std_err_1(2,i)*sin(t)];
%     
%     x2 = [mean_2(1,i)+std_err_2(1,i)*cos(t)];
%     y2 = [mean_2(2,i)+std_err_2(2,i)*sin(t)];
% 
% 
%     fill(x1,y1,'b', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)
%     fill(x2,y2,'r', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)
%     drawnow
% end
% plot([mean_1(1,1), mean_2(1,1)] , [mean_1(2,1), mean_2(2,1)],'yo', 'LineWidth', 3)
% plot([mean_1(1,201), mean_2(1,201)], [mean_1(2,201), mean_2(2, 201)], 'go', 'LineWidth', 3)
% title({'Parametric Plot of Mean Value of Lexical and Fingerspelling Gestures over Time on the'; 'First 2 Principle Components of all Linguistic ERPs in Significant Channels'})
% 
% end