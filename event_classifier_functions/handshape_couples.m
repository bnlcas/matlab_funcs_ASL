function [] = handshape_couples(ERPs, Data_Tag)
%% Function creates a listing the indecies of handshape pairings
movement_inds = find(strcmpi(ERPs.annot.movOnset,'M'));

%% Process handshape list
handshapes = lower(ERPs.annot.handshape);
handshapes = strrep(handshapes,' ','');


exclude = ~Data_Tag | strcmp(handshapes,'lax') | strcmpi(handshapes,'changing');


er_inds = [];
tr_inds = [];

pairings = {''};
count = 1;
for i = 1:length(movement_inds);
   
    if i == length(movement_inds)
        stop_pt = length(handshapes)-1;
    else
        stop_pt = movement_inds(i+1);
    end
    inds = movement_inds(i):stop_pt;
    
    handshapes_inds = handshapes(inds);
    block = exclude(inds);
    handshapes_inds(block) = []; % Kill non-linguisting ERPs
    handshapes_inds(strcmpi(handshapes_inds,'')) = [];
    
    if (length(handshapes_inds)-1) > 0
        for j = 1:(length(handshapes_inds)-1)
            pairings(count) = strcat(handshapes_inds(j),handshapes_inds(j+1));
            
            if strcmp(pairings(count),'er')
                ind = find(strcmp(handshapes(inds),'r'),1);
                
                er_inds = [er_inds, ind+movement_inds(i)];
            end
            if strcmp(pairings(count),'tr')
                ind = find(strcmp(handshapes(inds),'r'),1);
                tr_inds = [tr_inds, ind+movement_inds(i)];
            end
            count = count+1;
        end
    end
end



pairings;