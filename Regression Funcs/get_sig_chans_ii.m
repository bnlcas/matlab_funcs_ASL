function is_sig_chan = get_sig_chans_ii(regress_grid, varargin)
%% Takes the chans x timepts x vars+2 regress_grid and returns a chans x 1 boolean of channels 
% a certain number of consequtive sample point that are p significant
% if enter optional settings of p_thresh

if isempty(varargin)
    consequtive_sig_thresh = 10; % requires n time points to be significant
    p_thresh = 0.01/256; %/(length(sig_win));   % ARBITRARY
end
if length(varargin) == 2
    p_thresh = varargin{1};
    consequtive_sig_thresh = varargin{2};
elseif length(varargin) == 1
    p_thresh = varargin{1};
    consequtive_sig_thresh = 10;
end
%     consequtive_sig_thresh = varargin{1};
%     p_thresh = 0.05/400;


pval_grid = squeeze(regress_grid(:,:,end));
is_sig = pval_grid < p_thresh;
is_sig_chan = false(size(is_sig,1),1);

for i = 1:size(is_sig,1)
    consequtive_sigs = 0;
    for j = 1:size(is_sig,2)
        if is_sig(i,j)
            consequtive_sigs = consequtive_sigs + 1;
            if consequtive_sigs >= consequtive_sig_thresh
                is_sig_chan(i) = true;
            end
        else
            consequtive_sigs = 0;
        end
    end

end

a=1;