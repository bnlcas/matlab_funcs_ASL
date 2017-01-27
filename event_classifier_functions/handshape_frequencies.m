function [handshape_freq_fs] = handshape_frequencies(ERPs)

%phoneme recognition over time
% frequency of handshape occurrences as duplicates, fingerspelling, and comment categories:

handshape = ERPs.annot.handshape;
handshape = strrep(handshape ,' ','');
handshape = strrep(handshape ,' ','');
categories = unique(handshape );
categories = categories(~strcmp(categories,''));    % Clear blanks from categoriessparse_category = false(size(categories));          % Will remove categories with too few instances
category_size_thresh = 5;                           % Min of Three instances per category
for i = 1:length(categories)
	  cat_count = sum(strcmpi(handshape, categories(i)));
		    if cat_count < category_size_thresh
			 sparse_category(i) = true;
		    end
end
categories = categories(~sparse_category);

% make category frequency plot:
handshape_freq = zeros(size(categories));
for i = 1:length(categories)
	handshape_freq(i) = sum(strcmpi(handshape, categories(i)));
end
handshape_freq = handshape_freq./sum(handshape_freq);

is_fs = strcmpi(ERPs.annot.respType,'fs');
handshape_fs = ERPs.annot.handshape(is_fs);
handshape_fs = strrep(handshape_fs ,' ','');
handshape_fs = strrep(handshape_fs ,' ','');
handshape_freq_fs = zeros(size(categories));
for i = 1:length(categories)
	handshape_freq_fs(i) = sum(strcmpi(handshape_fs, categories(i)));
end
handshape_freq_fs = handshape_freq_fs./sum(handshape_freq_fs);

is_dup = strcmpi(ERPs.annot.respType, 'dup');
handshape_dup = ERPs.annot.handshape(is_dup);
handshape_dup = strrep(handshape_dup ,' ','');
handshape_dup = strrep(handshape_dup ,' ','');
handshape_freq_dup = zeros(size(categories));
for i = 1:length(categories)
	handshape_freq_dup(i) = sum(strcmpi(handshape_dup, categories(i)));
end
handshape_freq_dup = handshape_freq_dup./sum(handshape_freq_dup);


is_comment = strcmpi(ERPs.annot.respType, 'comment');
handshape_comment = ERPs.annot.handshape(is_comment);
handshape_comment = strrep(handshape_comment ,' ','');
handshape_comment = strrep(handshape_comment ,' ','');
handshape_freq_comment = zeros(size(categories));
for i = 1:length(categories)
	handshape_freq_comment(i) = sum(strcmpi(handshape_comment, categories(i)));
end
handshape_freq_comment = handshape_freq_comment./sum(handshape_freq_comment);
figure;
%X = [handshape_freq_fs, handshape_freq_dup, handshape_freq_comment];
%X = [handshape_freq, handshape_freq_fs];
X = [(handshape_freq_dup+handshape_freq_comment), handshape_freq_fs];
bar(X)
drawnow
axis tight;
set(gca,'xtick', 1:length(categories), 'xticklabel',categories);
ylabel('Lexical Frequency of HandShape')
%legend('FingerSpellings','Duplication', 'Comments')
legend('Duplications & Comments', 'FingerSpellings')
title('Freqneucy of Handshapes used in Different Response Types')
