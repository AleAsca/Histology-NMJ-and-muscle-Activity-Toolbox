%% ALL_INV_Dist_NMJ.m
%
% Created by: Manan Bhatt
%
% Date: 01/28/2025
%
% Version: 0.1.0
% *Requirements*: 
% 1. Original Muscle fig file with NMJs
% 2. Distance CSV file obtained from "Distance_NMJ_CSV_Gen.m"
%
% *Description*: This code extracts NMJ and electrode positions from a saved 
% figure and computes the sum of inverse distances from each electrode to all 
% NMJs, applying an exponential decay factor. The results are saved to a CSV 
% file, and a 3D visualization is generated with electrode sizes scaled by NMJ 
% contribution.

clc; clear; close all;

%% **Step 1: Load the Saved Figure**
figFile = "Figfile Address";
csvFileName = 'CSV_File Address';
figHandle = openfig(figFile, 'invisible'); % Load saved figure

% Extract all axes and graphical objects
axesHandle = findobj(figHandle, 'Type', 'axes');
allLines = findobj(axesHandle, 'Type', 'line');  % NMJs
allScatter = findobj(axesHandle, 'Type', 'scatter');  % Electrodes

%% **Step 2: Extract NMJ Positions**
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

numNMJs = length(muscleX);  % Total NMJs

%% **Step 3: Extract Electrode Positions**
electrodeX = get(allScatter, 'XData')';
electrodeY = get(allScatter, 'YData')';
electrodeZ = get(allScatter, 'ZData')';

% Ensure all electrode data are column vectors
electrodePositions = [electrodeX, electrodeY, electrodeZ];
numElectrodes = size(electrodePositions, 1); % Should be 64

%% **Step 4: Compute Sum of Inverse Distances**
sumInverseDistances = nan(numElectrodes, 1); % Initialize result vector
damp= -0.015; 

for i = 1:numElectrodes
    % Compute distances from electrode i to all NMJs
    distancesNMJ = sqrt((muscleX - electrodeX(i)).^2 + ...
                        (muscleY - electrodeY(i)).^2 + ...
                        (muscleZ - electrodeZ(i)).^2);
    
    % Avoid division by zero (if NMJ exactly at electrode position)
    distancesNMJ(distancesNMJ == 0) = NaN; 
    
    % Compute sum of inverse distances
    sumInverseDistances(i) = sum(exp(damp * distancesNMJ))

end

%% **Step 5: Save Results to CSV**
csvData = table(electrodeX, electrodeY, electrodeZ, sumInverseDistances, ...
                'VariableNames', {'Electrode_X', 'Electrode_Y', 'Electrode_Z', 'Avg_NMJ_Distance'});


writetable(csvData, csvFileName);
disp(['Saved sum of inverse NMJ distances for electrodes to ', csvFileName]);

%% **Step 6: Create a 3D Visualization**
figure;
hold on;
grid on;

% Plot NMJs (Gray Small Points)
scatter3(muscleX, muscleY, muscleZ, 10, [0.6, 0.6, 0.6], 'filled', 'MarkerFaceAlpha', 0.3); 

% Scale Electrode Size Based on Sum of Inverse Distances
minSize = 50;  % Smallest marker size
maxSize = 300; % Largest marker size
scaledSize = minSize + (sumInverseDistances - min(sumInverseDistances)) / ...
                        (max(sumInverseDistances) - min(sumInverseDistances)) * (maxSize - minSize);

% Plot Electrodes (Size = NMJ Contribution)
scatter3(electrodeX, electrodeY, electrodeZ, scaledSize, sumInverseDistances, 'filled'); 

% Customize Colorbar
colormap turbo;  % Change colormap (e.g., 'jet', 'parula', 'cool')
colorbar;
caxis([min(sumInverseDistances), max(sumInverseDistances)]); % Scale based on values

% Labels & Formatting
title('Electrode NMJ Contribution (Sum of Inverse Distances)');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
view(3); % 3D View
legend({'NMJs', 'Electrodes (Size = Contribution)'}, 'Location', 'best');

hold off;

%% **Step 7: Save the Plot**
savefig('Saved_FIG Address');
saveas(gcf, 'Saved_PNG Address .png');
disp('Saved NMJ contribution visualization using sum of inverse distances.');