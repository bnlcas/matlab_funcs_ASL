function hs_index = get_handshape_sequence(ERPs)
%% Function takes the non-grouped handshape and finds the order in which
% handshapes occur during linguistic signs

handshapes = strrep(ERPs.annot.handshape,' ','');
linguistic_start_index = find(strcmpi(ERPs.annot.lexTrans,'lexical'));

hs_index = zeros(size(handshapes));

for i = 1:(length(linguistic_start_index)-1)
    j = linguistic_start_index(i);
    sequence_num = 1;
    while j<linguistic_start_index(i+1) & strcmpi(ERPs.annot.filledLexTrans(j),'lexical')
        if ~strcmpi(handshapes(j),'') & ~strcmpi(handshapes(j),'changing')
            hs_index(j) = sequence_num;
            sequence_num = sequence_num + 1;
        end
        j = j+1;
    end
end
%% End Case
    j = linguistic_start_index(end);
    sequence_num = 1;
    while j<length(handshapes) & strcmpi(ERPs.annot.filledLexTrans(j),'lexical')
        if ~strcmpi(handshapes(j),'') & ~strcmpi(handshapes(j),'changing')
            hs_index(j) = sequence_num;
            sequence_num = sequence_num + 1;
        end
        j = j+1;
    end

    
end


