%% Heatmap_NMJ_Overlay_90_Shift.m
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
% *Description*: This code extracts muscle structure from a saved figure, 
% aligns a heatmap using PCA, and overlays it at a computed Z-layer after 
% rotating it 180 degrees. Electrodes are positioned without modification, 
% and the final visualization is saved with a fixed electrode layout.

%clc; clear; close all;

%% **Step 1: Load the .fig file**
figFile = "Figfile Address";
heatmapData = load("Raw Summed EMG data");
figHandle = openfig(figFile, 'invisible'); % Open in background

% Extract all axes and lines from the figure
axesHandle = findobj(figHandle, 'Type', 'axes');
allLines = findobj(figHandle, 'Type', 'line');

%% **Step 2: Extract X, Y, and Z data from the muscle structure**
muscleX = [];
muscleY = [];
muscleZ = [];

for i = 1:length(allLines)
    xData = get(allLines(i), 'XData');
    yData = get(allLines(i), 'YData');
    zData = get(allLines(i), 'ZData');
    
    if all(zData == zData(1))  % Store only the first detected layer (constant Z)
        muscleX = [muscleX, xData];
        muscleY = [muscleY, yData];
        muscleZ = [muscleZ, zData]; 
    end
end

% **Step 3: Identify the Z-layer spacing**
uniqueZ = unique(muscleZ);  % Unique Z-coordinates in sorted order
layerSpacing = mean(diff(uniqueZ));  % Compute average spacing between layers
disp(['Estimated layer spacing: ', num2str(layerSpacing), ' mm']);

% **Identify the topmost Z-layer**
topLayerZ = max(uniqueZ); 

% **Determine the Z-position for the heatmap**
heatmapZ = topLayerZ + layerSpacing;
disp(['Placing heatmap at Z = ', num2str(heatmapZ)]);

%% **Step 4: Compute PCA to Align the Heatmap**
musclePoints = [muscleX(:), muscleY(:)];
[coeff, score, ~] = pca(musclePoints);

principalVector = coeff(:,1);
muscleCenter = mean(musclePoints, 1);
muscleLength = max(score(:,1)) - min(score(:,1));

startPoint = muscleCenter - 0.3 * muscleLength * principalVector';
endPoint = startPoint + 0.6 * muscleLength * principalVector';

%% **Step 5: Load the Heatmap**
summedHeatmap = heatmapData.summedHeatmap;

[heatmapHeight, heatmapWidth] = size(summedHeatmap);
heatmapAspectRatio = heatmapWidth / heatmapHeight;

% **Compute heatmap corners (size unchanged)**
newWidth = norm(endPoint - startPoint);
newHeight = newWidth / heatmapAspectRatio;

% Compute perpendicular vector
perpendicularVector = [-principalVector(2), principalVector(1)];

corner1 = startPoint - (newHeight / 2) * perpendicularVector;
corner2 = startPoint + (newHeight / 2) * perpendicularVector;
corner3 = endPoint + (newHeight / 2) * perpendicularVector;
corner4 = endPoint - (newHeight / 2) * perpendicularVector;

% **Rotate the heatmap by 180 degrees (only the image, not placement)**
xGrid = linspace(corner1(1), corner3(1), heatmapWidth);  
yGrid = linspace(corner1(2), corner3(2), heatmapHeight);  
rotatedHeatmap = flipud(rot90(summedHeatmap, -1));  % Rotate heatmap 180 degrees


%% **Step 6: Generate Electrodes (Keep Positions Unchanged)**
numElectrodesX = 8; 
numElectrodesY = 8;

electrodeX = linspace(corner1(1), corner3(1), numElectrodesX); % **No change**
electrodeY = linspace(corner1(2), corner3(2), numElectrodesY); % **No change**
[Xe, Ye] = meshgrid(electrodeX, electrodeY);
electrodePositions = [Xe(:), Ye(:)];

% Assign electrodes to the heatmap layer
electrodeZ = heatmapZ * ones(size(electrodePositions, 1), 1);

%% **Step 7: Create the 3D Visualization**
figure;
copyobj(allLines, gca);  % Copy original muscle structures
hold on;

% **Overlay the rotated heatmap only at heatmapZ**
surface(xGrid, yGrid, heatmapZ * ones(heatmapHeight, heatmapWidth), ...
        'CData', rotatedHeatmap, 'FaceColor', 'texturemap', ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5); % Keep image orientation rotated

% **Plot Electrodes at the Heatmap Level (No Change in Position)**
scatter3(electrodePositions(:,1), electrodePositions(:,2), electrodeZ, ...
         50, 'bo', 'filled', 'MarkerEdgeColor', 'k'); % Blue electrodes

colormap hot;
colorbar;

% **Format the figure**
title('3D Muscle with 180-Degree Rotated Heatmap & Fixed Electrodes');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)'); 

view(3); % 3D View
grid on;
hold off;

% Save the final figure
savefig(figHandle, 'Destination Fig File Address');
saveas(figHandle, 'Destination PNG Address.png');
disp('Saved overlayed figure with 180-degree rotated heatmap and fixed electrodes.');