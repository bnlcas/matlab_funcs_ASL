function index = get_event_sequence(ERPs, event_label)
%% Function takes a boolean label whether and returns the the order in which
% labellings (boolean 1) occur during linguistic signs

linguistic_start_index = find(strcmpi(ERPs.annot.lexTrans,'lexical'));

index = zeros(size(event_label));

for i = 1:(length(linguistic_start_index)-1)
    j = linguistic_start_index(i);
    sequence_num = 1;
    while j<linguistic_start_index(i+1) & strcmpi(ERPs.annot.filledLexTrans(j),'lexical')
        if event_label(j)
            index(j) = sequence_num;
            sequence_num = sequence_num + 1;
        end
        j = j+1;
    end
end
%% End Case
    j = linguistic_start_index(end);
    sequence_num = 1;
    while j<length(event_label) & strcmpi(ERPs.annot.filledLexTrans(j),'lexical')
        if event_label(j)
            index(j) = sequence_num;
            sequence_num = sequence_num + 1;
        end
        j = j+1;
    end

    
end


