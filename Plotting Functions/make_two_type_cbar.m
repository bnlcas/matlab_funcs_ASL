function [] = make_two_type_cbar(bmax, rmax, null_gray)
%% Function generates a color bar for data that ranges from Negative to Positve Values
% The negative is in blue and the largest bound is given as bmax, while
% positive is in red and is given in terms of rmax

gray_middle = false; % Paramter to determine if the 0 intensity map is gray or black

saturation = max(abs([rmax, bmax]));    % The highest color level
grid = [bmax:(rmax-bmax)/400:rmax]';
grid = grid/saturation;

if ~gray_middle
    gray_thresh = 0;
    map = zeros(length(grid),3)/3;
    map((grid < -gray_thresh),:) = [zeros(sum(grid<-gray_thresh),1), zeros(sum(grid<-gray_thresh),1), abs(grid((grid < -gray_thresh)))];
    map((grid > gray_thresh),:) = [grid(grid > gray_thresh), zeros(sum(grid>gray_thresh),1), zeros(sum(grid>gray_thresh),1)];
else
    gray_thresh = 0.25;
    %gray_sat = min(2*gray_thresh,1); % 

    map = ones(length(grid),3)*gray_thresh;
    grid_blue = grid < -gray_thresh;
    grid_red = grid > gray_thresh;
    grid_grey = ~grid_blue & ~ grid_red;


    map(grid_blue,:) = [gray_thresh*(0:(sum(grid_blue)-1))'/(sum(grid_blue)-1), gray_thresh*(0:(sum(grid_blue)-1))'/(sum(grid_blue)-1), abs(grid((grid < -gray_thresh)))];
    map(grid_red,:) = [grid((grid > gray_thresh)), sort(gray_thresh*(0:(sum(grid_red)-1))','descend')/(sum(grid_red)-1), sort(gray_thresh*(0:(sum(grid_red)-1))'/(sum(grid_red)-1),'descend')];
    %map(grid_grey,:) = map(grid_grey,:) + 0.1*repmat(1+(cos(grid(grid_grey)*pi/gray_thresh)),1,3);
    map(grid_grey,:) = map(grid_grey,:) + 0.2*repmat(1+(cos(grid(grid_grey)*pi/(gray_thresh))),1,3);
end
test_cbar = true;
if test_cbar
    %gray_sat = min(2*gray_thresh,1); % 

    map = ones(length(grid),3)*null_gray;
    grid_blue = grid < 0;
    grid_red = grid > 0;

    map(grid_blue,:) = [null_gray*(0:(sum(grid_blue)-1))'/(sum(grid_blue)-1), null_gray*(0:(sum(grid_blue)-1))'/(sum(grid_blue)-1), (abs(bmax):-((abs(bmax)-null_gray)/(sum(grid_blue)-1)):null_gray)'];
    map(grid_red,:) = [(null_gray:((abs(rmax)-null_gray)/(sum(grid_red)-1)):rmax)', sort(null_gray*(0:(sum(grid_red)-1))','descend')/(sum(grid_red)-1), sort(null_gray*(0:(sum(grid_red)-1))'/(sum(grid_red)-1),'descend')];
    map(grid_red,:) = [(null_gray:((abs(rmax)-null_gray)/(sum(grid_red)-1)):rmax)', zeros(sum(grid_red),1), zeros(sum(grid_red),1)];

end



figure; imagesc(flipud(grid))
colormap(map)
colorbar
