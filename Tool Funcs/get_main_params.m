function [betas] = get_main_params(names, stat_grid)

secondary_param = true;
time_window = 150:350;

relevant_params = 1:50;
num_params = length(relevant_params);
%num_params = length(names);
betas = zeros(num_params,1);
dims = size(stat_grid);
for i = 1:dims(1)
    for j= time_window
        [~,ind]=max(abs(squeeze(stat_grid(i,j,relevant_params))));
        betas(ind) = betas(ind)+1;
        if secondary_param
            non_primary_params = relevant_params;
            non_primary_params(ind) = [];
            [~,ind]=max(abs(squeeze(stat_grid(i,j,non_primary_params))));
            betas(ind) = betas(ind)+1;
        end

    end
end

%Reformat Namings to replace '_' with ' '
for i = relevant_params
    mod = names{i};
    names{i} = strrep(mod,'_',' ');
end
figure;
barh(1:num_params, betas)
set(gca, 'YTick', 1:length(names), 'YTickLabel', names(relevant_params))
axis 'tight'

% ***Threshold Level:
% thresh = 1000;
% hold;
% plot([thresh thresh], [0.5, num_params+0.5],'k');

end