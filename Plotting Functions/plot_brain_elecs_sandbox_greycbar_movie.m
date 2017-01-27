function [F] = plot_brain_elecs_sandbox_greycbar_movie(ERPs, dat)
%% Function takes a 256xtimepts array and maps each timepoint as a frame in a movie

flatten = @(x) x(:);
%% Get Colormap:
%Settings in subfunction
%color_res = 500; %gray_level = 0.5;

% negative_rgb = [0.8 0.8 0];
% positive_rgb = [0 0.8 0.8]; % Traditional 

%negative_rgb = [0.8 0.25 0.4]; % Cyan to Purple (pseudo -real)
%positive_rgb = [0 0.8 0.8];
negative_rgb = [1 0.32 0.5]; % Cyan to Purple (pseudo -real)
positive_rgb = [0 1 1];

 %negative_rgb = [1 0 0];
 %positive_rgb = [0 0 1]; % Red to Blue (Traditional/ Translex)

[cmap, data_axis] = make_approx_cbar_var(dat, negative_rgb, positive_rgb);

sphere_size = 1.55;
 [sphx,sphy,sphz] = sphere(60);
 sphx = sphx*sphere_size; sphy = sphy*sphere_size; sphz = sphere_size*sphz;


%% Normalize Data
dat_raw = dat;
mask = (dat ~= 0); % is True for significant channels 

dat = dat/(max(abs(flatten(dat(mask)))));



subj = 'CH';
hem = 'lh';

roi = {'superiorfrontal','caudalmiddlefrontal','parstriangularis','parsopercularis','superiortemporal','precentral','postcentral','superiorparietal','supramarginal'};
clrs = linspace(-1,0,length(roi));

roi_label_flag = 0;
plot_elecs_flag = 1;
elec_label_flag = 0;
elecSize = 60;
offset = 5;


%%Changes!!!
rootdir = '/Users/changlab/Documents/changrepo/matlab/analysis/ASL';

%%

load([rootdir '/MRI/Meshes/' subj '_' hem '_pial.mat']);
load([rootdir '/MRI/elecs/hd_grid.mat']);

if strcmpi(hem,'lh')
    offset = offset * -1;
end

ctab_fid = fopen([rootdir '/MRI/aparc.annot.ctab']);
ctab = textscan(ctab_fid,'%d%s%d%d%d%d');
roi_names = ctab{2};
% roi = roi_names(2:end);
% clrs = linspace(-1,1,length(roi));

rndColorOrder = randperm(length(clrs));
clrs = clrs(rndColorOrder);
clrs(find(clrs == 0)) = 0.25;

fid = fopen([rootdir '/MRI/' hem '.aparc.annot.dpv']);
vert = textscan(fid,'%d%d%d%d%d');
vert = vert{5};

%cmap = cbrewer('seq','Reds',101);





for frame = 1:length(ERPs.time_axis)
    fig = figure;
    dat_frame = dat(:,frame);
    time_pt = ERPs.time_axis(frame);

ctmr_gauss_plot(cortex,[0 0 0],0,hem);
%view([-100,20])

kids = get(gca,'Children');


        
        
        

if roi_label_flag
    for i = 1:length(roi)
        roi_idx = find(strcmpi(roi_names,roi{i}));
        
        roi_verts{i} = find(vert == roi_idx);
        kids(2).FaceVertexCData(roi_verts{i},:) = repmat(clrs(i),length(roi_verts{i}),1);
    end
end
if plot_elecs_flag
    for i = 1:size(elecmatrix,1)
%         scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
%             elecSize,'y','filled');
         if mask(i)
%           scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
%              elecSize,cmap(round(dat(i)*100) + 1,:),'filled');
% 
                % [~,c_ind] = min(abs(data_axis - dat_frame(i)));
%                  scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
%                      elecSize, cmap(c_ind,:),'filled'); %   colorscaling3(i) colorscaling2(i)],'filled');
                % transp = abs(dat(i)/max(dat(:)));
                if dat_frame(i) >= 0
                    cval = cmap(length(data_axis),:);    
                else
                    cval = cmap(1,:);   
                end
                transp = abs(dat_frame(i)/max(abs(dat(:))));
%                    surf(sphx+elecmatrix(i,1)+offset,sphy+elecmatrix(i,2),sphz+elecmatrix(i,3), ...
%                        'FaceColor', cval,'FaceAlpha', transp , 'LineStyle', 'none','FaceLighting','none', 'BackFaceLighting', 'unlit')
                  surf(sphx+elecmatrix(i,1)+offset,sphy+elecmatrix(i,2),sphz+elecmatrix(i,3), ...
                       'FaceColor', cval,'FaceAlpha', transp , 'LineStyle', 'none','FaceLighting','none', 'BackFaceLighting', 'unlit')
        scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
                elecSize, [0.35 0.35 0.35]); %   colorscaling3(i) colorscaling2(i)],'filled');



         else
               scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
                elecSize, [0.35 0.35 0.35]); %   colorscaling3(i) colorscaling2(i)],'filled');

          end
             
        hold on;
        if elec_label_flag
            text(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
                num2str(i),'Color','y');
        end
    end
end

%title(['Difference in Mean High Gamma of Linguistic and Non-Linguistic ERPs at Time = ' num2str(time_pt), ' (ms)'])
%title(['Mean High Gamma of Linguistic ERPs at Time = ' num2str(time_pt), ' (ms)'])
%title(['Difference in Mean High Gamma of Real and Psuedo Signs at Time = ' num2str(time_pt), ' (ms)'])
title(['time = ' num2str(time_pt), ' (ms)'])
view([270 0])
create_movie_time_bar(time_pt)
frames(frame) = getframe(fig);
close(fig)
end

F = frames;

gen_movie = false;
movie_name = 'Citizen Galacticus.avi';
movie_directory = '/Users/changlab/Documents/changrepo/matlab/analysis/ASL/Graphics/Movies/Linguistic ERP';
if gen_movie
    cd(movie_directory)
    movie2avi(frames,movie_name)
end
    
