function [] = plot_respType_breakdown(ERPs)
% Creates a bar plot of the percentage of responses for each type of
% response to a stimuli
respTypes = ERPs.annot.respType(strcmpi(ERPs.annot.lexTrans,'lexical'));
resp_cats = unique(respTypes);
for i = 1:length(resp_cats)
    resp_freq(i) = sum(strcmpi(respTypes,resp_cats(i)));
end
resp_cats{5} = 'duplicate';
resp_cats{6} = 'fingerspelling';
resp_freq = resp_freq.*(100/length(respTypes));
[~,order] = sort(resp_freq,'descend');


figure; bar(resp_freq(order))
ylabel('Percentage of Responses')
xticklabel_rotate(1:9,45,resp_cats(order))
title('Frequency of Responses to ASL Stimuli')
