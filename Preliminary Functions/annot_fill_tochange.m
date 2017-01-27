function [filled_annot] = annot_fill_tochange(annot)
%% This function takes the Handshape, Location and Internal Movement Annotations
% and fills in empty values based on the entries of preceding values using
% the following rules:

% Internal Movements Annotations are only filled in for the rows of the 
% single Linguistic Gesture in which they occured (This owes to the
% sparsity of their encoding)
% Location and Handshapes rows that empty recieve the nearest non-empty
% annotation with a lower row number, Unless there are no annotations prior
% to this occurance, or unless the nearest non-empty annotation is
% 'changing' in which case leave blank.

filled_annot = annot;
rows = length(filled_annot.loc);

%% Fill Internal Movement Annotations:
filled_annot.intMov = strrep(filled_annot.intMov ,' ',''); % clean blanks
intMov_annots = find(~strcmpi(filled_annot.intMov,'')); % Listing of non-empty intMov
for i = 1:length(intMov_annots)
    intmov_ind = filled_annot.intMov(intMov_annots(i));
    gesture = filled_annot.filledLexTrans(intMov_annots(i)); % Gesture of current annotation - used to mark new gesture boundary
    index = intMov_annots(i)+1; % point after filled index

    while strcmpi(filled_annot.intMov(index),'') & strcmpi(filled_annot.filledLexTrans(index),gesture) & index<=rows
        filled_annot.intMov(index) = intmov_ind;
        index = index+1;
    end
    
end    

%% Fill Locations
filled_annot.loc = strrep(filled_annot.loc ,' ',''); % clean blanks

loc_annots = find(~strcmpi(filled_annot.loc,'') & ~strcmpi(filled_annot.loc,'changing')); % indecies of all definite location annotaion
for i = 1:length(loc_annots)
    loc_ind = filled_annot.loc(loc_annots(i));
    index = loc_annots(i) + 1;
    if index < rows
    while strcmpi(filled_annot.loc(index),'') & (index < rows)
        filled_annot.loc(index) = loc_ind;
        index = index + 1;
    end
    end
end

%% Fill Handshapes
filled_annot.handshape = strrep(filled_annot.handshape ,' ',''); % clean blanks

handshape_annots = find(~strcmpi(filled_annot.handshape,'') & ~strcmpi(filled_annot.handshape,'changing'));
for i = 1:length(handshape_annots)
    handshape_ind = filled_annot.handshape(handshape_annots(i));
    index = handshape_annots(i) + 1;
    if index < rows
    while strcmpi(filled_annot.handshape(index),'') & (index < rows)
        filled_annot.handshape(index) = handshape_ind;
        index = index+1;
    end
    end
end

filled_annot;
