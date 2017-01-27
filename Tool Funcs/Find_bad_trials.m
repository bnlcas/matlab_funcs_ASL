function badTrials = Find_bad_trials(datarootdir,subj,block,stimOnsetTime,timeLim)
% Written by Matt Sept 2014
% Changed by MAx on 9/16/2014

bad = load([datarootdir '/data/raw_data/' subj '/' subj '_' num2str(block) '/Artifacts/badTimeSegments.mat']);

badTrials = [];

for j = 1:length(stimOnsetTime)
    checkRange = [stimOnsetTime(j) + timeLim(1) ...
        stimOnsetTime(j) + timeLim(2)];
    if ~isempty(bad.badTimeSegments)
%         if any(bad.badTimeSegments(:,1) >= checkRange(1) & bad.badTimeSegments(:,2) <= checkRange(2) | ...
%                 bad.badTimeSegments(:,1) <= checkRange(1) & bad.badTimeSegments(:,2) >= checkRange(2))
%             badTrials = [badTrials j];
%         end
        if any(bad.badTimeSegments(:,1) >= checkRange(1) & bad.badTimeSegments(:,1) <= checkRange(2) | ...
                bad.badTimeSegments(:,2) <= checkRange(2) & bad.badTimeSegments(:,2) >= checkRange(1) | ...
                bad.badTimeSegments(:,1) <= checkRange(1) & bad.badTimeSegments(:,2) >= checkRange(2));
            badTrials = [badTrials j];
        end
    end
end
end