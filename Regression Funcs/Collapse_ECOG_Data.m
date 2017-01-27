function [collapsed_ecog] = Collapse_ECOG_Data(collapser, block_on_trigger, ecog)
%% Divide regression matrix into smaller blocks, divided when the collapser input is triggered
% this could also be stimOnset triggered by {'S'} or trans/lex...

% block trigger is one if any of the blockontrigger tags appear
block_trigger = false(length(collapser),1); 
for i = 1:length(block_on_trigger)
    block_trigger = (block_trigger | strcmpi(collapser, block_on_trigger(i)));
end

collapsed_ecog = ecog(:,:,block_trigger);
end