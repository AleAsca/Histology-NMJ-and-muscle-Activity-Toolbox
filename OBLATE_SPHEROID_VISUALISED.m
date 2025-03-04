%% OBLATE_SPHEROID_VISUALISED.m
%
% Created by: Manan Bhatt
%
% Date: 01/28/2025
%
% Version: 0.1.0
% *Requirements*: 
% 1. Muscle fig file with NMJs and electrode positions obtained from "Heatmap_NMJ_Overlay_90_Shift.m"
%
% *Description*: This code extracts neuromuscular junction (NMJ) and electrode 
% positions from a saved figure, identifies the electrode with the highest NMJ count, 
% and models its influence using a partial oblate spheroid. The spheroid is visualised 
% in 3D, highlighting NMJs within its region, and the results are saved as figures.

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

%% **Step 4: Set User-Controlled Spheroid Depth**
spheroidDepth = 0.7 * (topLayerZ - bottomLayerZ);  % Default to 50% depth (adjustable)
spheroidBottomZ = topLayerZ - spheroidDepth;  % Lowest extent of the oblate spheroid

%% Step 5: Extract Electrode Positions
electrodeX = get(allScatter, 'XData')';
electrodeY = get(allScatter, 'YData')';
electrodeZ = get(allScatter, 'ZData')';

% Ensure all electrode data are column vectors
electrodePositions = [electrodeX, electrodeY, electrodeZ];
numElectrodes = size(electrodePositions, 1);

%% Step 6: Compute Base Radius for Each Electrode's Spheroid
distances = pdist2(electrodePositions(:,1:2), electrodePositions(:,1:2)); % Compute 2D distances
distances(distances == 0) = Inf; % Ignore self-distance
baseRadius = min(distances, [], 2) / 2; % Half of the nearest electrode distance

%% Step 7: Find the Electrode with the Largest NMJ Count
nmjCounts = zeros(numElectrodes, 1);  % Store NMJ count for each electrode

for i = 1:numElectrodes
    % Define Electrode Position
    electrodeX_i = electrodePositions(i, 1);
    electrodeY_i = electrodePositions(i, 2);
    electrodeZ_i = electrodePositions(i, 3);

    % Find NMJs inside the **Partial Oblate Spheroid**
    insideSpheroid = false(size(muscleX));  % Logical mask

    for j = 1:length(muscleX)
        nmjZ_j = muscleZ(j);
        if nmjZ_j > topLayerZ || nmjZ_j < spheroidBottomZ
            continue; % Skip if NMJ is out of range
        end

        % **Compute spheroid radius at NMJ's Z-level**
        depthRatio = (topLayerZ - nmjZ_j) / spheroidDepth;
        spheroidRadiusAtZ = baseRadius(i) * sqrt(1 - depthRatio^2);  % Oblate spheroid equation

        % **Check if NMJ is inside spheroid**
        nmjDistanceXY = sqrt((muscleX(j) - electrodeX_i)^2 + (muscleY(j) - electrodeY_i)^2);
        if nmjDistanceXY <= spheroidRadiusAtZ
            insideSpheroid(j) = true;
        end
    end

    nmjCounts(i) = sum(insideSpheroid);  % Store NMJ count
end

% **Find the electrode with the highest NMJ count**
[~, maxNMJIndex] = max(nmjCounts);
electrodeX_max = electrodeX(maxNMJIndex);
electrodeY_max = electrodeY(maxNMJIndex);
electrodeZ_max = electrodeZ(maxNMJIndex);
baseRadius_max = baseRadius(maxNMJIndex);

%% **Step 8: Generate 3D Mesh for the Partial Oblate Spheroid**
theta = linspace(0, 2*pi, 30);  % Angular resolution
zValues = linspace(topLayerZ, spheroidBottomZ, 15);  % Vertical layers

% Compute shrinking radius at each depth
spheroidRadiusZ = baseRadius_max * sqrt(1 - ((topLayerZ - zValues) / spheroidDepth).^2);

% Convert to Cartesian coordinates
[Theta, Z_mesh] = meshgrid(theta, zValues);
X_mesh = spheroidRadiusZ(:) .* cos(Theta);
Y_mesh = spheroidRadiusZ(:) .* sin(Theta);

% Translate to electrode position
X_mesh = X_mesh + electrodeX_max;
Y_mesh = Y_mesh + electrodeY_max;

%% **Step 9: Find NMJs Inside This Spheroid**
insideSpheroidMax = false(size(muscleX));  % Logical mask for NMJs in selected electrode's spheroid

for j = 1:length(muscleX)
    nmjZ_j = muscleZ(j);
    if nmjZ_j > topLayerZ || nmjZ_j < spheroidBottomZ
        continue; % Skip if NMJ is out of range
    end

    % Compute spheroid radius at this depth
    depthRatio = (topLayerZ - nmjZ_j) / spheroidDepth;
    spheroidRadiusAtZ = baseRadius_max * sqrt(1 - depthRatio^2);

    % Check if NMJ is inside spheroid
    nmjDistanceXY = sqrt((muscleX(j) - electrodeX_max)^2 + (muscleY(j) - electrodeY_max)^2);
    if nmjDistanceXY <= spheroidRadiusAtZ
        insideSpheroidMax(j) = true;
    end
end

%% **Step 10: Plot the Spheroid & NMJs**
figure;
hold on;
grid on;

% Plot all NMJs (Gray Small Points)
scatter3(muscleX, muscleY, muscleZ, 10, [0.6, 0.6, 0.6], 'filled', 'MarkerFaceAlpha', 0.3);

% Plot NMJs inside the selected electrode's spheroid (Red)
scatter3(muscleX(insideSpheroidMax), muscleY(insideSpheroidMax), muscleZ(insideSpheroidMax), ...
         30, 'r', 'filled'); 

% Plot the selected electrode (Blue)
scatter3(electrodeX_max, electrodeY_max, electrodeZ_max, 150, 'b', 'filled', 'MarkerEdgeColor', 'k');

% Plot the spheroid as a wireframe
surf(X_mesh, Y_mesh, Z_mesh, 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'FaceColor', 'cyan');

% Customize View
title('Visualization of Largest NMJ Electrode Partial Oblate Spheroid');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
view(3);
legend({'All NMJs', 'NMJs in Spheroid', 'Electrode', 'Spheroid'}, 'Location', 'best');
hold off;

%% **Step 11: Save the Visualization**
savefig('Saved_FIG Address');
saveas(gcf, 'Saved_PNG Address .png');
disp('Saved largest NMJ electrode spheroid visualization.');