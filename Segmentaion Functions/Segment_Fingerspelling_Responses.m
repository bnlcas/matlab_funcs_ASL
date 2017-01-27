function [confusion_mat] = Segment_Fingerspelling_Responses(ERPs)


SVMModels = [];
% select relevant fingerspellings to partition
relevant_chars = ['B' 'C' 'L' 'O' 'S' 'Y'];
relevant_chans = 200:220;
relevant_timepts = 200:250;
for i = 1:length(relevant_chars)
    %get ecog data on a character
    is_char = (strcmpi(ERPs.annot.handshape, relevant_chars(i)));
    ecog_char = squeeze(mean(ERPs.ecog(relevant_chans ,relevant_timepts ,is_char),2))';
    %ecog_char = reshape(ERPs.ecog(relevant_chans ,relevant_timepts ,is_char),length(relevant_chans),[])';

    
    
    %get ecog for all other characters
    not_char = false(size(ERPs.annot.handshape));
    for j = 1:length(relevant_chars)
        if ~strcmpi(relevant_chars(j),relevant_chars(i))
            not_char = [not_char | strcmpi(ERPs.annot.handshape, relevant_chars(j))];
        end
    end
    ecog_not_char =  squeeze(mean(ERPs.ecog(relevant_chans ,relevant_timepts ,not_char),2))';
    %ecog_not_char =  reshape(ERPs.ecog(relevant_chans ,relevant_timepts ,not_char),length(relevant_chans),[])';
    X_Data = [ecog_char; ecog_not_char];
    Y_Data = [true(size(ecog_char,1),1); false(size(ecog_not_char,1),1)];
    
    
    % Fit Support Vector Machine:
    SVMModel = fitcsvm(X_Data,Y_Data,'KernelFunction','linear');
    CVSVMModel = crossval(SVMModel,'Leaveout', 'on');
    TrainedModel = CVSVMModel.Trained{1}
    SVMModels{i} = TrainedModel;    
end


% % GET CONFUSION MATRIX
 confusion_mat = zeros(length(relevant_chars));
is_char = false(length(ERPs.annot.handshape),1);

for i = 1:length(relevant_chars)
    is_char = strcmpi(ERPs.annot.handshape, relevant_chars(i));
    ecog_char = squeeze(mean(ERPs.ecog(relevant_chans ,relevant_timepts ,is_char),2))';
    %ecog_char = reshape(ERPs.ecog(relevant_chans ,relevant_timepts ,is_char),length(relevant_chans),[])';
    for j = 1:sum(is_char)
        ecog_pt = ecog_char(j,:);
        scores = zeros(length(relevant_chars),1);
            for k = 1:length(relevant_chars)
                [~,score] = predict(SVMModels{k},ecog_pt);
                scores(k) = score(2);
            end
        [~,class_index] = max(scores);
        confusion_mat(i,class_index) = confusion_mat(i,class_index)+1;
    end
    %confusion_mat(i,:)./sum(is_char);
    
end
confusion_mat = gdivide(confusion_mat, sum(confusion_mat,2));

%%Plot Confusion Matrix:
figure;

%BORROWED FROM STACK EXCHANGE:
imagesc(confusion_mat);            %# Create a colored plot of the matrix values
%colormap(gray);  %# Change the colormap to gray (so higher values are
colorbar;                       
textStrings = num2str(confusion_mat(:),'%0.2f');  %# Create strings from the matrix values to 2 decimals
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:length(relevant_chars));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center');   
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(confusion_mat(:) < midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
ax = gca;
ax.XTickLabel = {'B' 'C' 'L' 'O' 'S' 'Y'};
ax.YTickLabel = {'B' 'C' 'L' 'O' 'S' 'Y'};
title({'Confusion Matrix of SVM Classifier with Linear Kernel for Common FingerSpellings';'Using Mean ECoG HG of first 0.5 s on Ch 200 to Ch 256'})

end




    
    
    