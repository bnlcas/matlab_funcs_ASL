function [] = Plot_Sig_Chans(ERPs, Data_Tag, sig_chan, varargin )
%Temp Function to Plot the Mean ECoG of Significant Channels as given by
%256 element boolean array for weather a channel is significant

%inputs are Data tag which selects for certain classes of data to examine
%(lex vs transitional)
% Plots the differences given by a series of 256 element booleans listed in
% sig chan
% 

%% Chanels and time Parameters of plot
chans = sig_chan;
chan_numbers = find(chans);
twin = 101:301;
time_axis = ERPs.time_axis(twin);


%% Assemble Inputs
num_inputs = length(varargin);

means_mat = zeros(sum(chans),length(twin), num_inputs);
sems_mat = zeros(sum(chans),length(twin), num_inputs);
for i = 1:num_inputs
    ecog = ERPs.ecog(chans,twin,(Data_Tag & varargin{i}));
    means_mat(:,:,i) = mean(ecog,3);
end
for i = 1:num_inputs
    ecog = ERPs.ecog(chans,twin,(Data_Tag & varargin{i}));
    sems_mat(:,:,i) = nansem(ecog,3);
end

maxbound = max(max(max(means_mat + sems_mat)));
minbound = min(min(min(means_mat - sems_mat)));

%% Plot on significant channels
% Subplot dimensions are variable due to style considerations
cols = 8;
rows = 8;
grd = reshape(1:(cols*rows),rows,cols);
figure;
%colorlist = ['r', 'b', 'g', 'c', 'y'];
colorlist = [0.8 0 0; 0 0 0.8;0.8 0.8 0; 0 0.8 0.8; 0 0.8 0; 0.8 0 0.8];    % Plots of colors listed as Red, Blue, Yellow, Teal, Green, Purple


for i = 1:sum(chans)
    p = plotGridPosition(grd(i),cols*rows, cols);
    subplot('position',p);
%    subplot(rows,cols,i);
    for j = 1:num_inputs
        shadedErrorBar(ERPs.time_axis(twin), squeeze(means_mat(i,:,j)), squeeze(sems_mat(i,:,j)),{'color', colorlist(j,:)},1);
        hold on;
    end
    axis tight;
    set(gca,'YLim',[minbound maxbound]...
            ,'XTickLabel',[],'YTickLabel',[],'YTick',[minbound maxbound]);
        line([0 0],get(gca,'YLim'),'Color','k'); 
        line(get(gca,'XLim'),[0 0],'Color','k');
    text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim')) * 0.1)),num2str(chan_numbers(i)));
     
end
