% This script loads in group data and extracts N1 peak values for auditory
% and visual conditions.
clear
close all
clc

%% File paths
% Fieldtrip is necessary for the match_str function which finds the channel
% indices. The indices are also provided (see electrode section).
ft_path = '';  % your folder for the fieldtrip toolbox
addpath(ft_path);
ft_default.spmversion = 'spm12';
ft_defaults;

%% Cohen's D
% computeCohen_d(x1, x2, varargin) is necessary to calculate Cohen's D in
% the stats. Ruggero G. Bettinardi (2022). computeCohen_d(x1, x2, varargin)
%(https://www.mathworks.com/matlabcentral/fileexchange/62957-computecohen_d-x1-x2-varargin),
% MATLAB Central File Exchange. Retrieved May 4, 2022.

% Add the path here if you want to use it.
% addpath('') % Add to compute Cohen's d
%% Load the data
conditionlist = {
    'visual_active'
    'visual_passive'
    'auditory_active'
    'auditory_passive'
    };

vua = importdata(sprintf('%s.mat', conditionlist{1}));
vup = importdata(sprintf('%s.mat', conditionlist{2}));
aua = importdata(sprintf('%s.mat', conditionlist{3}));
aup = importdata(sprintf('%s.mat', conditionlist{4}));
%% Set up variables and parameters
participants = [1, 3:13, 15:17, 19:28, 30];
num_participants = length(participants);
vis_cond = {
    'vua'
    'vup'
    };
aud_cond = {
    'aua'
    'aup'
    };
% Empty variables to contain the peak values
num_conditions = length(vis_cond);
vis_peaks_N1 = zeros((num_participants), 3);
vis_peaks_N1(:,1) = participants;
aud_peaks_N1 = vis_peaks_N1;

% Define a function which takes a mean value for 1 participant across the selected channels and time window.
get_mean = @(data, chans, pp, t_start, t_end) squeeze(mean(data.individual(pp, chans, t_start:t_end), [2,3]));

% Electrodes
vis_electrodes = {'Oz', 'O1', 'O2'};
vis_chans = match_str(vua.label, vis_electrodes);
% vis_chans = [15, 16, 17]; % Use if you haven't got fieldtrip in the path.

target_electrodes = {'Cz', 'C3', 'C4'};
aud_chans = match_str(vua.label, target_electrodes);
% aud_chans = [8, 23, 24];

%% Visual
% Avg the ERP across all PPs in each condition
vis_ERPs = zeros(2, 375);
vis_ERPs(1,:) = squeeze(mean(vua.individual(:, vis_chans, :), [1 2]))';
vis_ERPs(2,:) = squeeze(mean(vup.individual(:, vis_chans, :), [1 2]))';

% Take the min value across all ERPs and get an index for the time
minval = zeros(1,height(vis_ERPs));
minidx = zeros(1,height(vis_ERPs));

for cond_idx = 1:height(vis_ERPs)
    [minval(cond_idx), minidx(cond_idx)] = min(vis_ERPs(cond_idx, :), [], 'all', 'linear');
end
[~, idx] = min(minval);

% Find the outer indices of the time window
centre_idx = minidx(idx);
centre_value_vis = vua.time(centre_idx);
t_start_vis = vua.time(centre_idx - 3);
t_end_vis = vua.time(centre_idx + 3);
t_range_vis = centre_idx - 3: centre_idx +3;

% Plot the result as a sanity check
figure(1)
for cond_idx = 1:length(vis_cond)
    plot(vua.time, vis_ERPs(cond_idx, :))
    hold on
end
% Plot vertical lines on the plot to show the time window.
xline(centre_value_vis, 'r')
xline(t_start_vis, 'k')
xline(t_end_vis, 'k')

% Print the values
for pp_idx = 1:length(participants)
    for cond_idx = 1:length(vis_cond)
        vis_peaks_N1(pp_idx, 2) = get_mean(vua, vis_chans, pp_idx, centre_idx-3, centre_idx+3);
        vis_peaks_N1(pp_idx, 3) = get_mean(vup, vis_chans, pp_idx, centre_idx-3, centre_idx+3);
    end
end

% writematrix(vis_peaks_N1, 'vis_peaks_N1.csv')
% save('vis_peaks_N1.mat', 'vis_peaks_N1')

%% Auditory
% Avg the ERP across all PPs in each condition
aud_ERPs = zeros(4, 375);
aud_ERPs(1,:) = squeeze(mean(aua.individual(:, aud_chans, :), [1 2]))';
aud_ERPs(2,:) = squeeze(mean(aup.individual(:, aud_chans, :), [1 2]))';

% Take the min value across all ERPs and get an index for the time
minval = zeros(1,height(aud_ERPs));
minidx = zeros(1,height(aud_ERPs));

for cond_idx = 1:height(aud_ERPs)
    [minval(cond_idx), minidx(cond_idx)] = min(aud_ERPs(cond_idx, :), [], 'all', 'linear');
end
[~, idx] = min(minval);

centre_idx = minidx(idx);
centre_value_aud = vua.time(centre_idx);
t_start_aud = vua.time(centre_idx - 3);
t_end_aud = vua.time(centre_idx + 3);
t_range_aud = centre_idx - 3: centre_idx +3;

% Sanity check
figure(2)
for cond_idx = 1:length(aud_cond)
    plot(aua.time, aud_ERPs(cond_idx, :))
    hold on
end
xline(centre_value_aud, 'r')
xline(t_start_aud, 'k')
xline(t_end_aud, 'k')

for pp_idx = 1:length(participants)
    for cond_idx = 1:length(aud_cond)
        aud_peaks_N1(pp_idx, 2) = get_mean(aua, vis_chans, pp_idx, centre_idx-3, centre_idx+3);
        aud_peaks_N1(pp_idx, 3) = get_mean(aup, vis_chans, pp_idx, centre_idx-3, centre_idx+3);
    end
end

% writematrix(aud_peaks_N1, 'aud_peaks_N1.csv')
% save('aud_peaks_N1.mat', 'aud_peaks_N1')
%% Statistics
active_vis = vis_peaks_N1(:,2);
passive_vis = vis_peaks_N1(:,3);
[vs.H_vis, vs.P_vis, vs.CI, vs.stats] = ttest(active_vis, passive_vis);
vs.d_vis = computeCohen_d(active_vis, passive_vis, 'paired');
vs.mean_act_vis = mean(active_vis);
vs.mean_pas_vis = mean(passive_vis);
vs.sd_act_vis = std(active_vis);
vs.sd_pas_vis = std(passive_vis);
vs.time_range_vis = (vua.time(t_range_vis)-0.1)*1000;
vs.t = vs.stats.tstat;
vs.df = vs.stats.df;
vs

% save('vis_stats_N1.mat', 'vs')

active_aud = aud_peaks_N1(:,2);
passive_aud = aud_peaks_N1(:,3);
[as.H_aud, as.P_aud, as.CI, as.stats] = ttest(active_aud, passive_aud);
as.d_aud = computeCohen_d(active_aud, passive_aud, 'paired');
as.mean_act_aud = mean(active_aud);
as.mean_pas_aud = mean(passive_aud);
as.sd_act_aud = std(active_aud);
as.sd_pas_aud = std(passive_aud);
as.time_range_aud = (aua.time(t_range_aud)-0.1)*1000;
as.t = as.stats.tstat;
as.df = as.stats.df;
as

% save('aud_stats_N1.mat', 'as')

%% Same stats but without participant 25, whose behavioural data were not saved.

vis_peaks_N1_no25 = vis_peaks_N1;
vis_peaks_N1_no25(22, :) = []; % Remove participant 25
active_vis = vis_peaks_N1_no25(:,2);
passive_vis = vis_peaks_N1_no25(:,3);
[vs.H_vis, vs.P_vis, vs.CI, vs.stats] = ttest(active_vis, passive_vis);
vs.d_vis = computeCohen_d(active_vis, passive_vis, 'paired');
vs.mean_act_vis = mean(active_vis);
vs.mean_pas_vis = mean(passive_vis);
vs.sd_act_vis = std(active_vis);
vs.sd_pas_vis = std(passive_vis);
vs.time_range_vis = (vua.time(t_range_vis)-0.1)*1000;
vs.t = vs.stats.tstat;
vs.df = vs.stats.df;
vs

% save('vis_stats_N1_no25.mat', 'vs')

aud_peaks_N1_no25 = aud_peaks_N1;
aud_peaks_N1_no25(22,:) = []; % Remove participant 25
active_aud = aud_peaks_N1_no25(:,2);
passive_aud = aud_peaks_N1_no25(:,3);
[as.H_aud, as.P_aud, as.CI, as.stats] = ttest(active_aud, passive_aud);
as.d_aud = computeCohen_d(active_aud, passive_aud, 'paired');
as.mean_act_aud = mean(active_aud);
as.mean_pas_aud = mean(passive_aud);
as.sd_act_aud = std(active_aud);
as.sd_pas_aud = std(passive_aud);
as.time_range_aud = (aua.time(t_range_aud)-0.1)*1000;
as.t = as.stats.tstat;
as.df = as.stats.df;
as

% save('aud_stats_N1_no25.mat', 'as')
