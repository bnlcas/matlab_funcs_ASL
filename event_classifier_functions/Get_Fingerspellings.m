%%Get Letters
function FS = Get_Fingerspellings(ERPs)
resp_type = ERPs.annot.respType;
HandShape = ERPs.annot.handshape; %is either length 1 or 

is_fingerspell = strcmpi(resp_type, 'fs'); % is a fingerspelling
% a single_letter is either 1 char
is_letter = false(size(HandShape));
for i = 1:length(HandShape)
    word = char(HandShape(i));
    is_one_char = (length(word) == 1);
    if (length(word) == 2)
        is_i_end = strcmpi(word(2), 'I');
    else
        is_i_end = false;
    end
    is_kp = strcmpi(word, 'K/P');
    is_letter(i) = ((is_one_char | is_i_end) | is_kp); 
end
FS = HandShape(is_letter & is_fingerspell);

letters = unique(FS);
occurances = zeros(size(letters));
for i = 1:length(letters);
    occurances(i) = sum(strcmpi(HandShape, letters(i)));
end
figure;
bar(1:length(letters), occurances)
set(gca, 'XTick', 1:length(letters), 'XTickLabel', letters)
axis 'tight'

