
test_entropy = zeros(36,83);
for j = 1:167
    for i = 1:72
        test_entropy(i,j) = entropy(test((5*(i-1)+1):(5*(i-1)+1), (5*(j-1)+1):(5*(j)+1)));
    end
end
    
% 
% function [SVM_weights_FS] = script(ERPs, sig_chan_loc, Data_Tag, y_data)
% folds = 100;
% betas = [];
% cmap = cbrewer('seq','Reds',101);
% 
% %           scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
% %               elecSize,cmap(round(dat(i)*100) + 1,:),'filled');
% for i=1:folds
%     X_Data = squeeze(mean(ERPs.ecog(sig_chan_loc,191:211,Data_Tag),2))';
%     inds = randperm(length(y_data));
%     max_ind = floor(((folds - 1)/folds)*length(y_data)/2);
%     X_Data_fold = X_Data(inds(1:max_ind),:);
%     y_data_fold = y_data(inds(1:max_ind));
%     
%     
%     SVMModel = fitcsvm(X_Data_fold,y_data_fold,'KernelFunction','linear', 'ClassNames',logical([0,1]), 'BoxConstraint', 0.01, 'KernelScale',1);
%     betas(i,:) = SVMModel.Beta;
%     
%     X_Data_hold = X_Data(inds(max_ind:end),:);
%     y_data_hold = y_data(inds(max_ind:end),:);
%     y_data_predict = predict(SVMModel, X_Data_hold);
%     accuracy(i) = sum(~xor(y_data_hold, y_data_predict))/length(y_data_hold);
% end
% SVM_weights_FS = mean(betas,1).^2;
% mean(accuracy)