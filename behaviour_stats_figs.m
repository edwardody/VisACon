% Edward Ody May 2022
% This script: 
% 1. loads total number of trials and proportion of 'second
% brighter' responses. 
% 2. Fits psychometric functions to the data and
% plots the results with one figure per participant showing both active and
% passive.
% 3. Extracts thresholds
% 4. Performs t-tests

%%
clear all
close all
%% Psignifit 4
% Psigniift 4 Matlab is required for this script. It can be downloaded from
% https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/neuronale-informationsverarbeitung/research/software/psignifit/
% Add the path where it located here:
addpath('')
%% Set general parameters and variables
ppnum = [1, 3:24, 26:28, 30]; % Participants to be included in the analysis

% Psignifit options
options = struct;
options.sigmoidName    = 'logistic'; % logistic sigmoid
options.expType        = 'equalAsymptote';
options.threshPC       = 0.5;

stimulus_intensities = (1:5);

% Variables to store the ouput
Vis_Thresh = zeros(27,3);
Vis_Thresh(:,1) = ppnum;
Aud_Thresh = Vis_Thresh;

% Plotting options
plotOptions.lineWidth = 2;
plotOptions.plotAsymptote = false;
plotOptions.plotThresh = false;
plotOptions.dataSize = 10;
linecol = [0 0 0; 1 0 0];
x_axis_lab = {'Brightness', 'Loudness'};
y_axis_lab = {'Proportion of ''2nd brighter'' responses','Proportion of ''2nd louder'' responses'};

plot_conds = {
    'Active'
    'Passive'
    };

%% Load files
load('visual_second.mat')
load('visual_total.mat')
load('auditory_second.mat')
load('auditory_total.mat')
%% Visual
close all
for pp_idx = 1:length(ppnum)
    % Set up the figure properties
    figure(pp_idx)
    % Plot the curves
    for cond_idx = 1:length(plot_conds)
        % Set up plotting conditions:
        results_vec = [];
        threshold_vec = [];
        slope_vec = [];
        
        plotOptions.lineColor = linecol(cond_idx, :);
        plotOptions.dataColor = linecol(cond_idx, :);
        plotOptions.xLabel = x_axis_lab{1};
        plotOptions.yLabel = y_axis_lab{1};

        number_of_second = vis_sec(:, cond_idx, pp_idx);
        number_of_trials = vis_total(:, cond_idx, pp_idx);
        
        data = [stimulus_intensities; number_of_second'; number_of_trials']';
        
        % Run the analysis
        result = psignifit(data,options);
        result.Slope = getSlopePC(result, 0.5);
        plotPsych(result, plotOptions);
        threshold = result.Fit(1) %Display the thresh for sanity check
        hold on
        
        % Print the result
        Vis_Thresh(pp_idx, cond_idx + 1) = result.Fit(1);
    end
    title(sprintf('Participant %d', ppnum(pp_idx)))
end

%% Auditory
close all
for pp_idx = 1:length(ppnum)
    figure(pp_idx)
    for cond_idx = 1:length(plot_conds)
        results_vec = [];
        threshold_vec = [];
        slope_vec = [];
        
        plotOptions.lineColor = linecol(cond_idx, :);
        plotOptions.dataColor = linecol(cond_idx, :);
        plotOptions.xLabel = x_axis_lab{2};
        plotOptions.yLabel = y_axis_lab{2};
        
        number_of_second = aud_sec(:, cond_idx, pp_idx);
        number_of_trials = aud_total(:, cond_idx, pp_idx);
        
        data = [stimulus_intensities; number_of_second'; number_of_trials']';
        
        result = psignifit(data,options);
        result.Slope = getSlopePC(result, 0.5);
        plotPsych(result, plotOptions);
        threshold = result.Fit(1)
        Aud_Thresh(pp_idx, cond_idx + 1) = result.Fit(1);
        hold on
    end
    title(sprintf('Participant %d', ppnum(pp_idx)))
end

%% Remove 14 and 18 due to nonfonforming thresholds
Vis_Thresh(Vis_Thresh(:, 1) == 14, :) = [];
Vis_Thresh(Vis_Thresh(:, 1) == 18, :) = [];
Aud_Thresh(Aud_Thresh(:, 1) == 14, :) = [];
Aud_Thresh(Aud_Thresh(:, 1) == 18, :) = [];
%% Perform t-tests
active_vis = Vis_Thresh(:,2);
passive_vis = Vis_Thresh(:, 3);
[vis.H_vis, vis.P_vis, vis.CI, vis.stats] = ttest(active_vis, passive_vis);
vis.d_vis = computeCohen_d(active_vis, passive_vis, 'paired');
vis.mean_act_vis = mean(active_vis);
vis.mean_pas_vis = mean(passive_vis);
vis.sd_act_vis = std(active_vis);
vis.sd_pas_vis = std(passive_vis);
vis.t = vis.stats.tstat;
vis.df = vis.stats.df;
vis

active_aud = Aud_Thresh(:,2);
passive_aud = Aud_Thresh(:, 3);
[aud.H_aud, aud.P_aud, aud.CI, aud.stats] = ttest(active_aud, passive_aud);
aud.d_aud = computeCohen_d(active_aud, passive_aud, 'paired');
aud.mean_act_aud = mean(active_aud);
aud.mean_pas_aud = mean(passive_aud);
aud.sd_act_aud = std(active_aud);
aud.sd_pas_aud = std(passive_aud);
aud.t = aud.stats.tstat;
aud.df = aud.stats.df;
aud