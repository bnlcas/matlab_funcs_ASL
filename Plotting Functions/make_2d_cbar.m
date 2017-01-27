function [] = make_2d_cbar()
gridsize = 100;
maxval = 0.316469; maxval = maxval/(gridsize);
color_grid = ones(gridsize); % create 100 x 100 matrix

for i = 1:gridsize
    for j = 1:gridsize
        color_grid(j,i) = maxval*(i + j/gridsize);
    end
end

map = [];
for i = 1:gridsize
    for j = 1:gridsize
        map = [map; [i/gridsize 0 j/gridsize]];
    end
end
figure; imagesc(flipud(color_grid'))
colormap(map)
