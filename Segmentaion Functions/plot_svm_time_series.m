function [] = plot_svm_time_series(ERPs, varargin) % confusion_mat_1) %, confusion_mat_2)
%% Function takes a variable number of frames x (n x n) confusion matries
% and plots the accuracy of these matricies over time

%% Set Time Axis for Frames
Analysis_frame = 10;                    % Average over 100ms window
time_advance = 5;                       % Advance the frame by 50ms each time
total_time_pts = size(ERPs.ecog,2);       % Standard for ECoG
frames = floor((total_time_pts - Analysis_frame)/time_advance) - 1; % clippout the last frame

time_axis_svm =  ERPs.time_axis(floor(((1:frames)-1)*time_advance + Analysis_frame/2)); % Takes the time closest to the midpoint of each analsis frame


num_conf_mats = length(varargin);
accuracy = zeros(1, frames);
figure; hold on;
%% For each Matrix Plot Accuracy
for i = 1:num_conf_mats
    confusion_mat = varargin{i};
    for k = 1:frames
        conf_frame = squeeze(confusion_mat(:,:,k));
            %if normalize
            %    conf_frame = gdivide(conf_frame, sum(conf_frame,2));
            %end
        accuracy(k) = trace(conf_frame)/sum(conf_frame(:));
    end
    plot(time_axis_svm, accuracy)
end

ax = gca;
ax.YLim(1)=0;
plot(ax.XLim, [1/size(conf_frame,1), 1/size(conf_frame,1)],'k--')
plot([0 0], ax.YLim, 'k','LineWidth',2)
title({'Time Evolution of Handshape Classification Accuracy';'On Channels'})
xlabel('Time From Handshape Onset (ms)')
ylabel('Classification Accuracy')
end
