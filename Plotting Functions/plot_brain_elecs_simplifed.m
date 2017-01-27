function [] = plot_brain_elecs_simplifed(electrodes_positions)
%% This is a simplified function designed to plot a grid of spherically
% shaped electrodes on an existing plot of a brain, on the basis of the
% input data 'electrodes_positions'
% 
% it is assumed that electrode positions is a (num_electrodes x 3) matrix,
% where the first, second and third columns correspond to the x, y and z coordinates respectively.
%
% it is possible to change the transparency and color of the electrodes,
% but this functionality has not been inculded here.

%% - Ben Lucas 4/2/2016

num_elecs = size(electrodes_positions,1);
%% Electrode Color Paramters
default_color = [0.8 0.2 0.1];
default_transparency = 0.4;

cmap = repmat(default_color, num_elecs, 1);         % Electrode by electrode color values
transp = repmat(default_transparency, num_elecs, 1);

%% Electrode shape parameters:
sphere_size = 1.65;            % Scaling for the size of the sphere
[sphx,sphy,sphz] = sphere(60); % Sphere with 60 pt mesh
sphx = sphx*sphere_size;       % Scale the sphere
sphy = sphy*sphere_size;
sphz = sphere_size*sphz;

offset = [-5 0 0];  % Offset term to prevent electrodes from collding with brain
% 

%% Plot the electrodes:
for i = 1:num_elecs
    hold on;
    surf(sphx+electrodes_positions(i,1)+offset(1),sphy+electrodes_positions(i,2)+offset(2),...
        sphz+electrodes_positions(i,3)+offset(3), ...
            'FaceColor', cmap(i,:),'FaceAlpha', transp(i) ,...
            'LineStyle', 'none','FaceLighting','none', 'BackFaceLighting', 'unlit')

end

end