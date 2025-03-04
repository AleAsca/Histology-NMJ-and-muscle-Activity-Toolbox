%% LOSSY_CONE_NMJ_Electrode_Point.m
%
% Created by: Manan Bhatt
%
% Date: 01/28/2025
%
% Version: 0.1.0
% *Requirements*: 
% 1. Muscle fig file with NMJs and electrode positions obtained from "Heatmap_NMJ_Overlay_90_Shift.m"
%
%
% *Description*: This code extracts NMJ and electrode positions from a saved 
% figure, assigns NMJs to electrodes using an inverted cone model, and applies 
% a lossy tissue attenuation model. The results, including NMJ-electrode mappings, 
% are saved to a CSV file and visualized in 3D with size-scaled electrodes.

clc; clear; close all;

%% **Step 1: Load the Saved Figure**
figFile = "Figfile Address";
figHandle = openfig(figFile, 'invisible');

% Extract all axes and graphical objects
axesHandle = findobj(figHandle, 'Type', 'axes');
allLines = findobj(figHandle, 'Type', 'line');  % NMJs
allScatter = findobj(figHandle, 'Type', 'scatter');  % Electrodes

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

% **Generate NMJ IDs**
numNMJs = length(muscleX);
nmjIDs = "NMJ_" + string(1:numNMJs);  % Assign unique NMJ IDs

% **Step 3: Identify Top and Bottom Z-Layers**
topLayerZ = max(muscleZ);
bottomLayerZ = min(muscleZ);
coneHeight = topLayerZ - bottomLayerZ;

%% **Step 4: Extract Electrode Positions**
electrodeX = get(allScatter, 'XData')';
electrodeY = get(allScatter, 'YData')';
electrodeZ = get(allScatter, 'ZData')';

% Ensure all electrode data are column vectors
electrodePositions = [electrodeX, electrodeY, electrodeZ];
numElectrodes = size(electrodePositions, 1);

%% **Step 5: Compute Base Radius for Each Cone**
distances = pdist2(electrodePositions(:,1:2), electrodePositions(:,1:2)); % Compute 2D distances
distances(distances == 0) = Inf; % Ignore self-distance
baseRadius = min(distances, [], 2) / 2; % Half of the nearest electrode distance

% **Define how much the cone expands as it goes deeper**
expansionFactor = max(baseRadius) * 4;

%% **Step 6: Loop Through Each Electrode and Compute NMJ Assignments**
damp = -0.015;  % Tissue attenuation factor

lossyNMJSum = nan(numElectrodes, 1);
nmjCounts = zeros(numElectrodes, 1);
nmjListPerElectrode = repmat({''}, numElectrodes, 1);  % Store NMJ IDs
electrodeListPerNMJ = repmat({''}, numNMJs, 1);  % Store Electrode IDs per NMJ

for i = 1:numElectrodes
    electrodeX_i = electrodePositions(i, 1);
    electrodeY_i = electrodePositions(i, 2);
    
    % **Compute Electrode Row & Column ID**
    colIdx = ceil(i / 8);  % Swap: Compute column index first
    rowIdx = 9 - mod(i - 1, 8) - 1;  % Flip row indexing (top row should be R1, bottom row should be R8)

    electrodeID = "R" + string(colIdx) + "_C" + string(rowIdx);
    % **Find NMJs inside the Inverted Cone**
    insideCone = false(size(muscleX));

    for j = 1:length(muscleX)
        nmjZ_j = muscleZ(j);
        if nmjZ_j > topLayerZ || nmjZ_j < bottomLayerZ
            continue; % Skip if NMJ is out of range
        end

        depthRatio = (topLayerZ - nmjZ_j) / coneHeight;
        coneRadiusAtZ = baseRadius(i) + (depthRatio * expansionFactor);
        
        nmjDistanceXY = sqrt((muscleX(j) - electrodeX_i)^2 + (muscleY(j) - electrodeY_i)^2);
        if nmjDistanceXY <= coneRadiusAtZ
            insideCone(j) = true;
            electrodeListPerNMJ{j} = strjoin(unique([electrodeListPerNMJ{j}, electrodeID]), ', ');
        end
    end

    % **Store NMJ Info**
    nmjCounts(i) = sum(insideCone);
    linkedNMJs = nmjIDs(insideCone);
    nmjListPerElectrode{i} = strjoin(linkedNMJs, ', ');  

    % **Compute Lossy Sum**
    nmjX = muscleX(insideCone);
    nmjY = muscleY(insideCone);
    nmjZ = muscleZ(insideCone);
    if isempty(nmjX)
        lossyNMJSum(i) = NaN;
    else
        nmjDistances = sqrt((nmjX - electrodeX_i).^2 + (nmjY - electrodeY_i).^2 + (nmjZ - topLayerZ).^2);
        lossyNMJSum(i) = sum(exp(damp * nmjDistances));
    end
end

%% **Step 7: Save CSV**
csvData = table(...
    strcat("R", string(ceil((1:numElectrodes)' / 8)), "_C", string(mod((1:numElectrodes)'-1, 8) + 1)), ...
    electrodeX, electrodeY, electrodeZ, baseRadius, lossyNMJSum, nmjCounts, nmjListPerElectrode, ...
    'VariableNames', {'Electrode_ID', 'Electrode_X', 'Electrode_Y', 'Electrode_Z', 'Base_Radius', 'Avg_NMJ_Distance', 'NMJ_Count', 'NMJ_List'});

writetable(csvData, 'Destination CSV Address');
disp('Saved CSV with Electrode IDs and NMJ Assignments.');

%% **Step 8: Create Improved 3D Visualization**
figure;
hold on;
grid on;

% **Plot NMJs (Larger and Different Color)**
scatter3(muscleX, muscleY, muscleZ, 25, 'm', 'filled', 'MarkerFaceAlpha', 0.5); 

% **Scale Electrode Size Based on NMJ Count**
minSize = 50;
maxSize = 300;
scaledSize = minSize + (nmjCounts - min(nmjCounts)) / (max(nmjCounts) - min(nmjCounts)) * (maxSize - minSize);

% **Plot Electrodes**
scatter3(electrodeX, electrodeY, electrodeZ, scaledSize, nmjCounts, 'filled'); 

% **Label NMJs with Electrode ID**
for j = 1:numNMJs
    text(muscleX(j), muscleY(j), muscleZ(j), char(electrodeListPerNMJ{j}), 'FontSize', 8, 'Color', 'k');
end

% **Customize Plot**
colormap turbo;
colorbar;
caxis([min(nmjCounts), max(nmjCounts)]);
title('Electrode NMJ Contribution (Lossy Tissue Model - Cone)');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
view(3);
legend({'NMJs', 'Electrodes (Size = NMJ Count)'}, 'Location', 'best');
hold off;

%% **Save the Plot**
savefig('Saved_FIG Address');
saveas(gcf, 'Saved_PNG Address .png');
disp('Saved improved 3D visualization.');