function class_acc = classification_accuracy(confusion_mat)
%% Function takes a N x N x Trials Confusion Matris and creates an array of the Accuracy Values of each TRIAL
trials = size(confusion_mat,3);
class_acc = zeros(1, trials);
for i = 1:trials
    mat_frame = squeeze(confusion_mat(:,:,i));
    class_acc(i) = trace(mat_frame)/sum(mat_frame(:));
end

end