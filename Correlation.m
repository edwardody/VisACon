% This script loads all the N1 and P2 peaks and thresholds for correlation analysis.

%% set up preliminaries
clear all
close all
clc

%% Set up variables and parameters
ppnum = [1, 3:13, 15:17, 19:24, 26:28, 30]; % Participants to be included in the analysis

%% Set up file paths
% ERP Path- add the directory where you have stored ERP peak values.
addpath('')

% Behaviour Path- add the directory where you have stored thresholds.
addpath('')

%% Load the data
load('auditory_N1.mat')
load('visual_N1.mat')
load('auditory_P2.mat')
load('visual_P2.mat')

% Behaviour is missing pp25 so remove this from neural before further
% analysis.
aud_peaks_N1(22, :) = [];
aud_peaks_P2(22, :) = [];
vis_peaks_N1(22, :) = [];
vis_peaks_P2(22, :) = [];

load('auditory_thresholds.mat')
load('visual_thresholds.mat')

%% Calculate suppression
aud_N1_sup = aud_peaks_N1(:, 2) - aud_peaks_N1(:, 3);
aud_P2_sup = aud_peaks_P2(:, 2) - aud_peaks_P2(:, 3);
vis_N1_sup = vis_peaks_N1(:, 2) - vis_peaks_N1(:, 3);
vis_P2_sup = vis_peaks_P2(:, 2) - vis_peaks_P2(:, 3);

vis_beh_sup = Vis_Thresh(:, 2) - Vis_Thresh(:, 3);
aud_beh_sup = Aud_Thresh(:, 2) - Aud_Thresh(:, 3);
%% Correlation and plots
% Auditory, N1, threshold
figure(1)
[a1.r, a1.p] = corrcoef(aud_N1_sup, aud_beh_sup);
scatter(aud_N1_sup, aud_beh_sup)

% Auditory, P2, threshold
figure(2)
[a2.r, a2.p] =corrcoef(aud_P2_sup, aud_beh_sup);
scatter(aud_P2_sup, aud_beh_sup)

% Visual, N1, threshold
figure(3)
[v1.r, v1.p] = corrcoef(vis_N1_sup, vis_beh_sup);
scatter(vis_N1_sup, vis_beh_sup)

% Visual, P2, threshold
figure(4)
[v2.r, v2.p] = corrcoef(vis_P2_sup, vis_beh_sup);
scatter(vis_P2_sup, vis_beh_sup)

