% This script plots the main ERP plot, topoplots for N1 and P2 and
% colorbars for visual conditions.

%% set up preliminaries
clear all
close all
clc

%% File paths
% Fieldtrip is necessary for the match_str function which finds the channel
% indices. The indices are also provided (see electrode section).
ft_path = '';  % your folder for the fieldtrip toolbox
addpath(ft_path);
ft_default.spmversion = 'spm12';
ft_defaults;

%% Bounded line
% boundedline.m is necessary for the ERP plots. Kelly Kearney (2022). boundedline.m (https://github.com/kakearney/boundedline-pkg),
% GitHub. Retrieved May 4, 2022.
% Add it to the path here.
addpath(genpath(''))

%% Colourmap
% cbrewer : colorbrewer schemes for Matlab is necessary for the topoplots
%Charles (2022). cbrewer : colorbrewer schemes for Matlab (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab),
% MATLAB Central File Exchange. Retrieved May 4, 2022.
set(0,'DefaultFigureColormap',brewermap(64,'*RdBu'));
 
%% Load data
conditionlist = {
    'visual_active'
    'visual_passive'
    };

vua = importdata(sprintf('%s.mat', conditionlist{1}));
vup = importdata(sprintf('%s.mat', conditionlist{2}));
%% Set up variables and parameters
participants = [1, 3:13, 15:17, 19:28, 30];
num_pp = length(participants);

get_mean = @(data, channels) squeeze(mean(data.individual(:, channels, :), [1 2]));
get_SD = @(data, channels) squeeze(std(mean(data.individual(:,channels,:),2)))/sqrt(30);

chans = match_str(vua.label, {'Oz', 'O1', 'O2'});
% chans = [15, 16, 17]; % Use if you haven't got fieldtrip in the path.

%% Create the ERP Plot
plot_conds = {'vua', 'vup'};
line_col = {'k','r'};

close all
h = figure(1);

% Plot the ERPs with shaded areas showing standard deviation
[line1, hp1] = boundedline(vua.time, get_mean(eval(plot_conds{1}), chans), get_SD(eval(plot_conds{1}), chans),...
    line_col{1}, 'alpha');
xlim([-0.2, 0.6])
[line2, hp2] = boundedline(vua.time, get_mean(eval(plot_conds{2}), chans), get_SD(eval(plot_conds{2}), chans),...
    line_col{2}, 'alpha');
xlim([-0.2, 0.6])
line1.LineWidth = 1;
line2.LineWidth = 1;

% Display green shaded areas showing the time window for analysis
N1bars = [0.184 0.208];
hold on
patch([N1bars(1) N1bars(1), N1bars(2) N1bars(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], 'g','LineStyle','none')
alpha(0.3)

N2bars = [0.232 0.256];
hold on
patch([N2bars(1) N2bars(1), N2bars(2) N2bars(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], 'g','LineStyle','none')
alpha(0.3)

set(gca,'children',flipud(get(gca,'children')))
set(gca,'YDir','reverse')
xline(0.1)
yline(0)
xticks(-0.2:0.1:0.8)
xticklabels([-300, -200, -100, 0, 100, 200, 300, 400, 500, 600, 700])
legend([line1, line2], 'Active', 'Passive')
legend('boxoff')
legend('Location', 'southeast')
set(gca,'FontSize', 13)
title('Visual (Oz, O1, O2)', 'FontSize', 13)
xlabel('Time(ms)')
ylabel('Mean Amplitude (μV)')
h.Position = [100 100 585 360];
set(gca,'fontname','Helvitica')

% Save
% f = gcf;
% exportgraphics(f,'Visual ERP.png','ContentType','vector', 'resolution', '300')

%% N1 Topoplots
% Get the time index for N1
vis_ERPs = zeros(2, 375);
vis_ERPs(1,:) = squeeze(mean(vua.individual(:, chans, :), [1 2]))';
vis_ERPs(2,:) = squeeze(mean(vup.individual(:, chans, :), [1 2]))';

% Take the max value across all ERPs and get an index for the time
minval = zeros(1,height(vis_ERPs));
minidx = zeros(1,height(vis_ERPs));

for cond_idx = 1:height(vis_ERPs)
    [minval(cond_idx), minidx(cond_idx)] = min(vis_ERPs(cond_idx, :), [], 'all', 'linear');
end
[~, idx] = min(minval);

centre_idx_N1 = minidx(idx);
t_range_N1 = centre_idx_N1 - 3: centre_idx_N1 +3;
t_range_N1 = vua.time(t_range_N1);

h2 = figure(2);
cfg                     = [];
cfg.parameter           = 'individual';
cfg.xlim                = [t_range_N1(1), t_range_N1(end)];
cfg.zlim                = [-1.5, 2];
cfg.colorbar            = 'no';
cfg.layout              = 'EEG1020.lay';
cfg.comment             = 'no';
cfg.markersize          = 3;
cfg.markerfontsize      = 15;
ft_topoplotER(cfg, vua)
title('Active', 'FontSize', 13)
h2.Position = [100 100 240 120];
% f = gcf;
% exportgraphics(f,'Active N1 Topo.png','ContentType','vector', 'resolution', '300')

h3 = figure(3);
ft_topoplotER(cfg, vup)
title('Passive', 'FontSize', 13)
h3.Position = [100 100 240 120];
% f = gcf;
% exportgraphics(f,'Passive N1 Topo.png','ContentType','vector', 'resolution', '300')

%% P2 Topoplots

% Get the time index for P2
vis_ERPs = zeros(2, 375);
vis_ERPs(1,:) = squeeze(mean(vua.individual(:, chans, :), [1 2]))';
vis_ERPs(2,:) = squeeze(mean(vup.individual(:, chans, :), [1 2]))';

% Take the max value across all ERPs and get an index for the time
maxval = zeros(1,height(vis_ERPs));
maxidx = zeros(1,height(vis_ERPs));

for cond_idx = 1:height(vis_ERPs)
    [maxval(cond_idx), maxidx(cond_idx)] = max(vis_ERPs(cond_idx, :), [], 'all', 'linear');
end
[~, idx] = max(maxval);

centre_idx_P2 = maxidx(idx);
t_range_P2 = centre_idx_P2 - 3: centre_idx_P2 +3;
t_range_P2 = vua.time(t_range_P2);

h4 = figure(4);
cfg                     = [];
cfg.parameter           = 'individual';
cfg.xlim                = [t_range_P2(1), t_range_P2(end)];
cfg.zlim                = [-2, 8];
cfg.colorbar            = 'no';
cfg.layout              = 'EEG1020.lay';
cfg.comment             = 'no';
cfg.markersize          = 3;
cfg.markerfontsize      = 15;
ft_topoplotER(cfg, vua)
title('Active', 'FontSize', 13)
h4.Position = [100 100 240 120];

% f = gcf;
% exportgraphics(f,'Active P2 Topo.png','ContentType','vector', 'resolution', '300')

h5 = figure(5);
ft_topoplotER(cfg, vup)
title('Passive', 'FontSize', 13)
h5.Position = [100 100 240 120];

% f = gcf;
% exportgraphics(f,'Passive P2 Topo.png','ContentType','vector', 'resolution', '300')

%% N1 colorbar
% close all
h6 = figure(6);
ax = axes;
c = colorbar(ax);
ax.Visible = 'off';
caxis([-1.5, 2])
box('off')
c.Ticks = [-1.5, 2];
c.Location = 'west';
set(gca,'FontSize', 13)
h6.Position = [100 100 140 140];

% f = gcf;
% exportgraphics(f,'N1 Colorbar.png','ContentType','vector', 'resolution', '300')

%% P2 colorbar
close all
h7 = figure(7);
ax = axes;
c = colorbar(ax);
ax.Visible = 'off';
caxis([-2 8])
box('off')
c.Ticks = [-2 8];
c.Location = 'west';
set(gca,'FontSize', 13)
h7.Position = [100 100 140 140];

% f = gcf;
% exportgraphics(f,'P2 Colorbar.png','ContentType','vector', 'resolution', '300')