%% RAW_HEATMAP_COMBINED.m
%
% Created by: Manan Bhatt
%
% Date: 01/28/2025
%
% Version: 0.1.0
% *Requirements*: 
% 1. Original Individual Muscle heatmaps
%
% *Description*: This code loads a processed Motor Unit (MU) data matrix from 
% a .mat file, extracts a specific range of data, and computes a summed 
% heatmap across multiple layers. The resulting heatmap is displayed as an 
% image and saved both as a .mat file and a normalized .png image for visualization.

clc; clear; close all;

% Load the processed .mat file
data = load('Original Individual EMG .mat file Address');
MUs = data.MUs; % Extract the matrix
filtered = MUs(:,100:550,100:600);
[numLayers, rows, cols] = size(filtered);

% Sum all 5 layers
summedHeatmap = sum(filtered, 1); % Sum across layers
summedHeatmap = squeeze(summedHeatmap); % Remove singleton dimension

% Display the summed heatmap
figure;
imagesc(summedHeatmap);
colormap hot;
colorbar;
title('Summed Heatmap of All 5 Layers');
axis equal tight;

% Save summed heatmap as a .mat file
save('Destination address as .mat', 'summedHeatmap');
disp('Saved summed heatmap as Summed_Heatmap.mat');

% Normalize for saving as PNG (avoid black images)
minVal = min(summedHeatmap(:));
maxVal = max(summedHeatmap(:));
normalizedSummed = (summedHeatmap - minVal) / (maxVal - minVal); % Scale to [0,1]

% Save as PNG
imwrite(normalizedSummed, 'Destination address as .png');
disp('Saved summed heatmap as Summed_Heatmap.png');