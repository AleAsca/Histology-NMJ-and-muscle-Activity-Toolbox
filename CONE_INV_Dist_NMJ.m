%% CONE_INV_Dist_NMJ.m
%
% Created by: Manan Bhatt
%
% Date: 01/28/2025
%
% Version: 0.1.0
% *Requirements*: 
% 1. Muscle fig file with NMJs and electrode positions obtained from "Heatmap_NMJ_Overlay_90_Shift.m"
% 2. Distance CSV file obtained from "Distance_NMJ_CSV_Gen.m"
%
% *Description*: This code extracts NMJ and electrode positions from a saved 
% figure and computes NMJ contributions for each electrode using an inverted 
% cone model. It calculates inverse NMJ distances, saves the results to a CSV 
% file, and visualizes electrode contributions in 3D with size-scaled markers.

%clc; clear; close all;

%% Step 1: Load the Saved Figure
figFile = "Figfile Address";
csvFileName = "CSV_File Address";
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

%% Step 6: Loop Through Each Electrode and Compute NMJ Distances Inside the Cone
avgNMJDistances = nan(numElectrodes, 1);  % Ensure the correct number of rows
nmjCounts = zeros(numElectrodes, 1);  % Store the number of NMJs per electrode

for i = 1:numElectrodes
    % Define Electrode Position
    electrodeX_i = electrodePositions(i, 1);
    electrodeY_i = electrodePositions(i, 2);
    
    % Find NMJs inside the **Inverted Cone**
    insideCone = false(size(muscleX));  % Logical mask for NMJs inside the cone

    for j = 1:length(muscleX)
        % Compute expanding cone radius at NMJ's Z-level
        nmjZ_j = muscleZ(j);
        if nmjZ_j > topLayerZ || nmjZ_j < bottomLayerZ
            continue; % Skip if NMJ is out of range
        end
        
        % **Calculate radius of the cone at depth Z_j**
        depthRatio = (topLayerZ - nmjZ_j) / coneHeight;  % Normalize depth (0 = top, 1 = bottom)
        coneRadiusAtZ = baseRadius(i) + (depthRatio * expansionFactor);  % Expanding radius
        
        % **Check if NMJ is inside cone radius at this depth**
        nmjDistanceXY = sqrt((muscleX(j) - electrodeX_i)^2 + (muscleY(j) - electrodeY_i)^2);
        if nmjDistanceXY <= coneRadiusAtZ
            insideCone(j) = true;
        end
    end

    % Compute 3D distances for NMJs inside the cone
    nmjX = muscleX(insideCone);
    nmjY = muscleY(insideCone);
    nmjZ = muscleZ(insideCone);
    nmjCounts(i) = sum(insideCone);  % Store NMJ count

    if isempty(nmjX)
        avgNMJDistances(i) = NaN;  % No NMJs found in this cone
    else
        nmjDistances = sqrt((nmjX - electrodeX_i).^2 + (nmjY - electrodeY_i).^2 + (nmjZ - topLayerZ).^2);
        nmjDistances = 1./(nmjDistances);
        avgNMJDistances(i) = sum(nmjDistances);
        %avgNMJDistances(i) = mean(nmjDistances);
        %avgNMJDistances(i) = median(nmjDistances);
    end
end

%% **Step 7: Save Results to CSV**
csvData = table(electrodeX, electrodeY, electrodeZ, baseRadius, avgNMJDistances, nmjCounts, ...
                'VariableNames', {'Electrode_X', 'Electrode_Y', 'Electrode_Z', 'Base_Radius', 'Avg_NMJ_Distance', 'NMJ_Count'});


writetable(csvData, csvFileName);
disp(['Saved NMJ distances for electrodes (cone model) to ', csvFileName]);

%% **Step 8: Create a 3D Visualization**
figure;
hold on;
grid on;

% Plot NMJs (Gray Small Points)
scatter3(muscleX, muscleY, muscleZ, 10, [0.6, 0.6, 0.6], 'filled', 'MarkerFaceAlpha', 0.3); 

% **Scale Electrode Size Based on NMJ Count**
minSize = 50;  % Smallest electrode marker size
maxSize = 300; % Largest electrode marker size
scaledSize = minSize + (nmjCounts - min(nmjCounts)) / (max(nmjCounts) - min(nmjCounts)) * (maxSize - minSize);

% Plot Electrodes (Size = NMJ Contribution)
scatter3(electrodeX, electrodeY, electrodeZ, scaledSize, nmjCounts, 'filled'); 

% Customize Colorbar
colormap turbo;  % Change colormap (e.g., 'jet', 'parula', 'cool')
colorbar;
caxis([min(nmjCounts), max(nmjCounts)]); % Scale based on NMJ counts

% Labels & Formatting
title('Electrode NMJ Contribution (Inverted Cone)');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
view(3); % 3D View
legend({'NMJs', 'Electrodes (Size = NMJ Count)'}, 'Location', 'best');

hold off;

%% **Step 9: Save the Plot**
savefig('Saved_FIG Address');
saveas(gcf, 'Saved_PNG Address .png');
disp('Saved NMJ contribution visualization (cone model).');

