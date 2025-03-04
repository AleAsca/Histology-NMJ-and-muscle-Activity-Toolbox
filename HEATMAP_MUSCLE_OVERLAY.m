%% HEATMAP_MUSCLE_OVERLAY.m
%
% Created by: Manan Bhatt
%
% Date: 01/28/2025
%
% Version: 0.1.0
% *Requirements*: 
% 1. Original Muscle fig file with NMJs
% 2. RAW summed EMG Heatmap obtained from "RAW_HEATMAP_COMBINED.m"
%
% *Description*: This code overlays a rotated heatmap onto the first Z-layer 
% of a muscle structure extracted from a saved figure. It aligns the heatmap 
% using PCA, ensures visibility with transparency adjustments, and places fixed 
% electrode positions. The updated figure is saved with enhanced visualization.

clc; clear; close all;

%% Step 1: Load Original Figure (Keep Everything Else Unchanged)
figFile = "FigFile Address";
heatmapData = load("Raw Summed EMG data");
figHandle = openfig(figFile, 'invisible'); % Open original figure in the background

% Get current figure axis
axesHandle = findobj(figHandle, 'Type', 'axes');

%% Step 2: Identify First Z-Layer for Heatmap Placement
allLines = findobj(figHandle, 'Type', 'line');
% **Change NMJ Color to Blue**
for i = 1:length(allLines)
    set(allLines(i), 'Color', [0, 0, 0]); % Change NMJ color to blue
    set(allLines(i), 'LineWidth', 1.5); % Make NMJs more visible
end

muscleX = [];
muscleY = [];
muscleZ = [];

for i = 1:length(allLines)
    xData = get(allLines(i), 'XData');
    yData = get(allLines(i), 'YData');
    zData = get(allLines(i), 'ZData');
    
    if all(zData == zData(1))  
        muscleX = [muscleX, xData];
        muscleY = [muscleY, yData];
        muscleZ = [muscleZ, zData];
    end
end

firstLayerZ = max(muscleZ);  

%% Step 3: PCA for Heatmap Alignment
musclePoints = [muscleX(:), muscleY(:)];
[coeff, score, ~] = pca(musclePoints);

principalVector = coeff(:,1);
muscleCenter = mean(musclePoints, 1);
muscleLength = max(score(:,1)) - min(score(:,1));

startPoint = muscleCenter - 0.3 * muscleLength * principalVector';
endPoint = startPoint + 0.6 * muscleLength * principalVector';

%% Step 4: Load Heatmap & Define Placement
summedHeatmap = heatmapData.summedHeatmap;
summedHeatmap = rot90(summedHeatmap, -1);  % Rotate heatmap 180 degrees

[heatmapHeight, heatmapWidth] = size(summedHeatmap);
heatmapAspectRatio = heatmapWidth / heatmapHeight;

newWidth = norm(endPoint - startPoint);
newHeight = newWidth / heatmapAspectRatio;

perpendicularVector = [-principalVector(2), principalVector(1)];

corner1 = startPoint - (newHeight / 2) * perpendicularVector;
corner2 = startPoint + (newHeight / 2) * perpendicularVector;
corner3 = endPoint + (newHeight / 2) * perpendicularVector;
corner4 = endPoint - (newHeight / 2) * perpendicularVector;

[xGrid, yGrid] = meshgrid(linspace(corner1(1), corner3(1), heatmapWidth), ...
                          linspace(corner1(2), corner3(2), heatmapHeight));

heatmapZ = firstLayerZ * ones(size(xGrid));

%% Step 5: Overlay Heatmap on First Layer
figure(figHandle);
hold on;

h = surface(xGrid, yGrid, heatmapZ, ...
        'CData', flipud(summedHeatmap), ...
        'FaceColor', 'texturemap', ...
        'EdgeColor', 'none', ...
        'FaceLighting', 'none');

% **Fix 1: Ensure Heatmap Visibility**
alpha(h, 0.5); % Make heatmap fully opaque

% **Fix 2: Apply Correct Colormap**
colormap hot;
caxis([min(summedHeatmap(:)), max(summedHeatmap(:))]); % Ensure correct color scaling
colorbar;

%% Step 6: Place Electrodes
numElectrodesX = 8;
numElectrodesY = 8;

electrodeX = linspace(corner1(1), corner3(1), numElectrodesX);
electrodeY = linspace(corner1(2), corner3(2), numElectrodesY);

[Xe, Ye] = meshgrid(electrodeX, electrodeY);
electrodePositions = [Xe(:), Ye(:)];
electrodeZ = firstLayerZ * ones(size(electrodePositions, 1), 1);

scatter3(electrodePositions(:,1), electrodePositions(:,2), electrodeZ, ...
         80, 'green', 'filled', 'MarkerEdgeColor', 'k');

%% Step 7: Save the Updated Figure
savefig(figHandle, 'Destination Fig File Address');
saveas(figHandle, 'Destination PNG Address.png');

disp('Updated figure saved with visible heatmap and electrodes.');