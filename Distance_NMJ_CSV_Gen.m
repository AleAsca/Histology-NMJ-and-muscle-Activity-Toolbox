%% Distance_NMJ_CSV_Gen.m
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
% figure and computes the average NMJ distance for each electrode using a 
% cylindrical region. The results, including electrode-specific NMJ distances, 
% are saved to a CSV file for further analysis.

clc; clear; close all;

%% Step 1: Load the Saved Figure
figFile = "Figfile Address";
figHandle = openfig(figFile, 'invisible'); % Load saved figure

% Extract all axes and graphical objects
axesHandle = findobj(figHandle, 'Type', 'axes');
allLines = findobj(axesHandle, 'Type', 'line');  % NMJs
allScatter = findobj(axesHandle, 'Type', 'scatter');  % Electrodes

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

% **Step 3: Identify the Topmost Z-Layer**
lastLayerZ = max(muscleZ);  % Topmost muscle Z-layer

%% Step 4: Extract Electrode Positions
electrodeX = get(allScatter, 'XData')';
electrodeY = get(allScatter, 'YData')';
electrodeZ = get(allScatter, 'ZData')';

% Ensure all electrode data are column vectors
electrodePositions = [electrodeX, electrodeY, electrodeZ];
numElectrodes = size(electrodePositions, 1);

%% Step 5: Compute Cylinder Radius
distances = pdist2(electrodePositions(:,1:2), electrodePositions(:,1:2)); % Compute 2D distances
distances(distances == 0) = Inf; % Ignore self-distance
cylinderRadius = min(distances, [], 2) / 2; % Half of the nearest electrode distance

%% Step 6: Loop Through Each Electrode and Compute NMJ Distances
avgNMJDistances = nan(numElectrodes, 1);  % Ensure the correct number of rows

for i = 1:numElectrodes
    % Define Cylinder Region
    electrodeX_i = electrodePositions(i, 1);
    electrodeY_i = electrodePositions(i, 2);
    radius_i = cylinderRadius(i);
    
    % Find NMJs inside the cylinder (2D radius)
    distancesNMJ = sqrt((muscleX - electrodeX_i).^2 + (muscleY - electrodeY_i).^2);
    insideCylinder = distancesNMJ <= radius_i;  % Logical mask

    % Compute 3D distances for valid NMJs
    nmjX = muscleX(insideCylinder);
    nmjY = muscleY(insideCylinder);
    nmjZ = muscleZ(insideCylinder);

    if isempty(nmjX)
        avgNMJDistances(i) = NaN;  % No NMJs found in this cylinder
    else
        nmjDistances = sqrt((nmjX - electrodeX_i).^2 + (nmjY - electrodeY_i).^2 + (nmjZ - electrodeZ(i)).^2);
        avgNMJDistances(i) = mean(nmjDistances);
    end
end

%% **Step 7: Save Results to CSV**
csvData = table(electrodeX, electrodeY, electrodeZ, cylinderRadius, avgNMJDistances, ...
                'VariableNames', {'Electrode_X', 'Electrode_Y', 'Electrode_Z', 'Cylinder_Radius', 'Avg_NMJ_Distance'});

csvFileName = 'Destination Address';
writetable(csvData, csvFileName);

disp(['Saved NMJ distances for electrodes to ', csvFileName]);