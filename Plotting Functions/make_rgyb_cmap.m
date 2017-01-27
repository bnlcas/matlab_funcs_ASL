function map = make_rgyb_cmap(dat, gray_val)
%% This function takes a 1d array of data and creates
% a nx3 set of RGB triples defined such that negative values are blue,
% positive values are red, and values below grey_thresh (the fraction of
% the normalized intensity) are set to be gray
map = ones(length(dat),3)*gray_val;
if size(dat,2) > size(dat,1)
    dat = dat';
end

%% get the max and min ranges of the spectrum
rval = max([max(dat), 0]);
%bval = min([min(dat),0]);
bval = min(dat);
norm = max(abs([rval, bval]));

%% Normalize Data
dat_raw = dat;
dat = dat/norm; rval = rval/norm; bval = bval/norm;



%% Set Blue values:
is_blue = dat < 0; blue_dat = abs(dat(is_blue));
blue_ch = blue_dat*(abs(bval)-gray_val)/abs(bval) + gray_val;
red_ch = gray_val - gray_val*blue_dat/abs(bval);
green_ch = red_ch;
map(is_blue,:) = [red_ch, green_ch, blue_ch];

%% Set Red Values:

is_red = dat > 0; red_dat = abs(dat(is_red));
%red_ch = red_dat*(abs(rval)-gray_val)/abs(rval) + gray_val;
red_ch = red_dat*(1-gray_val)/(rval - bval) + (rval*gray_val-bval)/(rval-bval);
blue_ch = gray_val*red_dat/(abs(rval)-bval) - gray_val*rval/(rval-bval);
green_ch = blue_ch;
map(is_red,:) = [red_ch, green_ch, blue_ch];

plot_cbar = true;
if plot_cbar
    figure; imagesc(flipud(dat_raw))
    colormap(map)
    colorbar
end
