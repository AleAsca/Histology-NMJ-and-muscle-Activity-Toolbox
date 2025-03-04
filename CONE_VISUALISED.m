%% CONE_VISUALISED.m
%
% Created by: Manan Bhatt
%
% Date: 01/28/2025
%
% Version: 0.1.0
% *Requirements*: 
% 1. Muscle fig file with NMJs and electrode positions obtained from "Heatmap_NMJ_Overlay_90_Shift.m"
%
% *Description*: This code extracts NMJ and electrode positions from a saved 
% figure and identifies the electrode with the highest NMJ count. It models 
% the NMJ distribution using an inverted cone, visualises the structure in 3D, 
% and saves the resulting figure for further analysis.

clc; clear; close all;

%% Step 1: Load the Saved Figure
figFile = "Figfile Address";
figHandle = openfig(figFile, 'invisible'); % Load saved figure

% Extract all axes and graphical objects
axesHandle = findobj(figHandle, 'Type', 'axes');
allLines = findobj(figHandle, 'Type', 'line');  % NMJs
allScatter = findobj(figHandle, 'Type', 'scatter');  % Electrodes

%% Step 2: Extract NMJ Positions
muscleX = [];
muscleY = [];
muscleZ = [];

for i = 1:length(allLines)
    xData = get(allLines(i), 'XData');
    yData = get(allLines(i), 'YData');
    zData = get(allLines(i), 'ZData');

    if all(zData == zData(1))  % Ensure it's a 2D NMJ layer
        muscleX = [muscleX, xData];
        muscleY = [muscleY, yData];
        muscleZ = [muscleZ, zData];
    end
end

% **Step 3: Identify the Topmost and Bottommost Z-Layers**
topLayerZ = max(muscleZ);  % Topmost muscle Z-layer
bottomLayerZ = min(muscleZ);  % Bottommost muscle layer
coneHeight = topLayerZ - bottomLayerZ;  % Total height of the cone

%% Step 4: Extract Electrode Positions
electrodeX = get(allScatter, 'XData')';
electrodeY = get(allScatter, 'YData')';
electrodeZ = get(allScatter, 'ZData')';

% Ensure all electrode data are column vectors
electrodePositions = [electrodeX, electrodeY, electrodeZ];
numElectrodes = size(electrodePositions, 1);

%% Step 5: Compute Base Radius for Each Cone
distances = pdist2(electrodePositions(:,1:2), electrodePositions(:,1:2)); % Compute 2D distances
distances(distances == 0) = Inf; % Ignore self-distance
baseRadius = min(distances, [], 2) / 2; % Half of the nearest electrode distance

% **Define how much the cone expands as it goes deeper**
expansionFactor = max(baseRadius) * 4; % Adjust this to control widening speed

%% Step 6: Find the Electrode with the Largest NMJ Count
nmjCounts = zeros(numElectrodes, 1);  % Store the number of NMJs per electrode
for i = 1:numElectrodes
    % Define Electrode Position
    electrodeX_i = electrodePositions(i, 1);
    electrodeY_i = electrodePositions(i, 2);

    % Find NMJs inside the **Inverted Cone**
    insideCone = false(size(muscleX));  % Logical mask

    for j = 1:length(muscleX)
        nmjZ_j = muscleZ(j);
        if nmjZ_j > topLayerZ || nmjZ_j < bottomLayerZ
            continue; % Skip if NMJ is out of range
        end
        
        % **Compute cone radius at NMJ's Z-level**
        depthRatio = (topLayerZ - nmjZ_j) / coneHeight;  % Normalize depth (0 = top, 1 = bottom)
        coneRadiusAtZ = baseRadius(i) + (depthRatio * expansionFactor);  % Expanding radius
        
        % **Check if NMJ is inside the cone**
        nmjDistanceXY = sqrt((muscleX(j) - electrodeX_i)^2 + (muscleY(j) - electrodeY_i)^2);
        if nmjDistanceXY <= coneRadiusAtZ
            insideCone(j) = true;
        end
    end

    nmjCounts(i) = sum(insideCone);  % Store NMJ count
end

% **Find the electrode with the highest NMJ count**
[~, maxNMJIndex] = max(nmjCounts);
electrodeX_max = electrodeX(maxNMJIndex);
electrodeY_max = electrodeY(maxNMJIndex);
electrodeZ_max = electrodeZ(maxNMJIndex);
baseRadius_max = baseRadius(maxNMJIndex);

%% **Step 7: Generate 3D Mesh for the Cone**
theta = linspace(0, 2*pi, 30);  % Angular resolution
zValues = linspace(topLayerZ, bottomLayerZ, 20);  % Vertical layers

% Compute expanding radius at each depth (MATCHING THE ORIGINAL CONE EQUATION)
coneRadiusZ = baseRadius_max + ((topLayerZ - zValues) / coneHeight) * expansionFactor;

% Convert to Cartesian coordinates
[Theta, Z_mesh] = meshgrid(theta, zValues);
X_mesh = coneRadiusZ(:) .* cos(Theta);
Y_mesh = coneRadiusZ(:) .* sin(Theta);

% Translate to electrode position
X_mesh = X_mesh + electrodeX_max;
Y_mesh = Y_mesh + electrodeY_max;

%% **Step 8: Find NMJs Inside This Cone**
insideConeMax = false(size(muscleX));  % Logical mask for NMJs in selected electrode's cone

for j = 1:length(muscleX)
    nmjZ_j = muscleZ(j);
    if nmjZ_j > topLayerZ || nmjZ_j < bottomLayerZ
        continue; % Skip if NMJ is out of range
    end

    % Compute cone radius at this depth
    depthRatio = (topLayerZ - nmjZ_j) / coneHeight;
    coneRadiusAtZ = baseRadius_max + (depthRatio * expansionFactor);

    % Check if NMJ is inside the cone
    nmjDistanceXY = sqrt((muscleX(j) - electrodeX_max)^2 + (muscleY(j) - electrodeY_max)^2);
    if nmjDistanceXY <= coneRadiusAtZ
        insideConeMax(j) = true;
    end
end

%% **Step 9: Plot the Cone & NMJs**
figure;
hold on;
grid on;

% Plot all NMJs (Gray Small Points)
scatter3(muscleX, muscleY, muscleZ, 10, [0.6, 0.6, 0.6], 'filled', 'MarkerFaceAlpha', 0.3);

% Plot NMJs inside the selected electrode's cone (Red)
scatter3(muscleX(insideConeMax), muscleY(insideConeMax), muscleZ(insideConeMax), ...
         30, 'r', 'filled'); 

% Plot the selected electrode (Blue)
scatter3(electrodeX_max, electrodeY_max, electrodeZ_max, 150, 'b', 'filled', 'MarkerEdgeColor', 'k');

% Plot the cone as a wireframe
surf(X_mesh, Y_mesh, Z_mesh, 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'FaceColor', 'cyan');

% Customize View
title('Visualization of Largest NMJ Electrode Cone');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
view(3);
legend({'All NMJs', 'NMJs in Cone', 'Electrode', 'Cone'}, 'Location', 'best');
hold off;

%% **Step 10: Save the Visualization**
savefig('Saved_FIG Address');
saveas(gcf, 'Saved_PNG Address .png');
disp('Saved largest NMJ electrode cone visualization.');