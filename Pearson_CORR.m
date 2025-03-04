%% Pearson_CORR.m
%
% Created by: Manan Bhatt
%
% Date: 01/28/2025
%
% Version: 0.1.0
% *Requirements*: 
% 1. RAW summed EMG Heatmap obtained from "RAW_HEATMAP_COMBINED.m"
% 2. Individual Distance CSV files from all the shapes to be viewed from :
% "Electrode_NMJ_Distances.csv", "Electrode_NMJ_SumInverseDistances.csv", "Electrode_NMJ_Distances_Spheroid.csv", "Electrode_NMJ_Distances_Cone.csv", "Electrode_NMJ_Distances_Lossy_Cone_With_NMJs.csv"
%
% *Description*: This code rotates a heatmap 180 degrees, extracts average 
% activation at electrode positions, and computes Pearson correlations with 
% NMJ distances from various models. Results are visualized with scatter plots, 
% and electrode-NMJ data is saved to CSV.

%clc; clear; close all;

%% **Step 1: Load & Rotate the Original Heatmap**
heatmapData = load('Raw Summed EMG data');
summedHeatmap = heatmapData.summedHeatmap; % Heatmap matrix

% **Rotate the heatmap by 180 degrees**
summedHeatmap = flipud(rot90(summedHeatmap, -1));  

% Get heatmap dimensions
[heatmapHeight, heatmapWidth] = size(summedHeatmap);

%% **Step 2: Generate the 8×8 Electrode Mesh (Row-wise from Top-Left)**
numElectrodes = 8; % 8 per row/column (7 segments)
regionSize = 40; % 40x40 region around each electrode

% Define evenly spaced electrode positions (row-wise, starting from top-left)
xElectrode = round(linspace(1, heatmapWidth, numElectrodes));
yElectrode = round(linspace(1, heatmapHeight, numElectrodes));
[X_mesh, Y_mesh] = meshgrid(xElectrode, yElectrode);

% Convert meshgrid to column vectors (64 electrodes, **ordered row-wise**)
electrodeX = X_mesh(:);
electrodeY = Y_mesh(:);

%% **Step 3: Compute 40×40 Average Intensity for Each Electrode (Row-wise)**
heatmapActivation = zeros(size(electrodeX));

for i = 1:length(electrodeX)
    % Define region bounds (ensure it doesn't go out of bounds)
    xMin = max(1, electrodeX(i) - floor(regionSize/2));
    xMax = min(heatmapWidth, electrodeX(i) + floor(regionSize/2));
    yMin = max(1, electrodeY(i) - floor(regionSize/2));
    yMax = min(heatmapHeight, electrodeY(i) + floor(regionSize/2));
    
    % Extract region and compute average intensity from **rotated heatmap**
    region = summedHeatmap(yMin:yMax, xMin:xMax);
    heatmapActivation(i) = mean(region(:), 'omitnan'); % Avoid NaN issues
end

%% **Step 4: Load NMJ Distance Data from All Three Models**
models = {'Cylinder', 'Cone', 'Spheroid', 'All', 'Lossy Cone'};
csvFiles = {...
    'Cylinder CSV Address', ...
    'Cone CSV Address', ...
    'Spheroid CSV Address' ...
    'All NMJs CSV Address' ...
    'Lossy Cone CSV Address' ...
};

correlationResults = table('Size', [3 3], 'VariableTypes', {'string', 'double', 'double'}, ...
                           'VariableNames', {'Model', 'Pearson_R', 'P_Value'});

%% **Step 5: Compute Pearson Correlation for Each Model and Plot Separately**
for i = 1:length(models)
    % Load CSV file
    nmjData = readtable(csvFiles{i});
    
    % Ensure electrodes are sorted row-wise
    nmjData = sortrows(nmjData, {'Electrode_Y', 'Electrode_X'}); 

    % Extract NMJ distances in the correct order
    avgNMJDistances = nmjData.Avg_NMJ_Distance;
    
    % Remove NaN values
    validIdx = ~isnan(avgNMJDistances) & ~isnan(heatmapActivation);
    avgNMJDistances = avgNMJDistances(validIdx);
    heatmapActivationValid = heatmapActivation(validIdx);

    % Compute Pearson correlation
    [R, pValue] = corrcoef(avgNMJDistances, heatmapActivationValid);
    R_value = R(1,2); % Extract correlation coefficient
    p_value = pValue(1,2);
    
    % Store results in a table
    correlationResults.Model(i) = models{i};
    correlationResults.Pearson_R(i) = R_value;
    correlationResults.P_Value(i) = p_value;
    
    % **Display results in console**
    disp(['Model: ', models{i}]);
    disp(['   Pearson Correlation (R) = ', num2str(R_value, '%.4f')]);
    disp(['   P-Value = ', num2str(p_value, '%.4f')]);

    %% **Step 6: Create a Separate Figure for Each Model with Clickable Electrode ID**
    figure;
    scatter(avgNMJDistances, heatmapActivationValid, 80, 'b', 'filled'); 
    hold on;
    
    % **Compute and Plot Linear Fit**
    p = polyfit(avgNMJDistances, heatmapActivationValid, 1);
    xFit = linspace(min(avgNMJDistances), max(avgNMJDistances), 100);
    yFit = polyval(p, xFit);
    plot(xFit, yFit, 'r', 'LineWidth', 2); % Trendline
    
    % **Labels & Formatting**
    xlabel('Lossy Tissue NMJ Distance');
    ylabel('Heatmap Activation Value (40x40 Avg)');
    title([models{i}, ': NMJ Distance vs. Heatmap Activation (R = ', num2str(R_value, '%.2f'), ')']);
    legend({'Electrodes', 'Trendline'}, 'Location', 'best');
    grid on;
    hold off;
    
    %% **Enable Clickable Data Cursor for Electrode ID**
    dcm = datacursormode(gcf);
    set(dcm, 'UpdateFcn', @(obj, event_obj) displayElectrodeID(event_obj, nmjData));
    
    %% **Save Each Figure with Clickable IDs**
    savefig(['Destination Address/NMJ_vs_Heatmap_', models{i}, '_Rotated.fig']);
    saveas(gcf, ['Destination Address/NMJ_vs_Heatmap_', models{i}, '_Rotated.png']);
    disp(['Saved NMJ distance vs. heatmap activation plot for ', models{i}]);

    %% **Step 7: Save Each Figure Separately**
    savefig(['Destination Address/NMJ_vs_Heatmap_', models{i}, '_Rotated.fig']);
    saveas(gcf, ['Destination Address/NMJ_vs_Heatmap_', models{i}, '_Rotated.png']);
    disp(['Saved NMJ distance vs. heatmap activation plot for ', models{i}]);
end

%% **Step 8: Save Correlation Results**
writetable(correlationResults, 'Destination Address/Correlation_Results_Rotated.csv');
disp('Saved correlation results as CSV.');

%% **Step 9: Display Table with Correlation Results**
disp('Pearson Correlation Results (Rotated Heatmap):');
disp(correlationResults);


%% **Custom Data Cursor Function (Electrode ID & NMJ List Display)**
function txt = displayElectrodeID(event_obj, nmjData)
    % Get clicked data point index
    pos = get(event_obj, 'Position');
    x_clicked = pos(1);  % NMJ Distance
    y_clicked = pos(2);  % Heatmap Activation
    
    % Find the closest electrode match
    [~, idx] = min(abs(nmjData.Avg_NMJ_Distance - x_clicked));  % Find closest NMJ
    
    % Extract the Electrode ID (Row and Column)
    numElectrodes = 8;  % 8x8 Grid
    rowIdx = ceil(idx / numElectrodes);  % Compute row index
    colIdx = mod(idx - 1, numElectrodes) + 1;  % Compute column index
    electrodeID = "R" + string(rowIdx) + "_C" + string(colIdx);
    
    % Extract NMJs associated with this electrode
    nmjsLinked = nmjData.NMJ_List(idx);  % Get the NMJs for this electrode
    
    % Display Data in the Tooltip
    txt = {['NMJ Distance: ', num2str(x_clicked, '%.2f')], ...
           ['Heatmap Activation: ', num2str(y_clicked, '%.2f')], ...
           ['Electrode ID: ', electrodeID], ...
           ['Linked NMJs: ', nmjsLinked]};
end


%% **Step 10: Save Electrode Data with NMJs in CSV**

% Read NMJ distance data
nmjCSV = 'Address/ Electrode_NMJ_Distances_Lossy_Cone.csv';
nmjData = readtable(nmjCSV);

% Ensure NMJ data is sorted correctly
nmjData = sortrows(nmjData, {'Electrode_Y', 'Electrode_X'});

% Extract NMJ Distance
correctedNMJDistances = nmjData.Avg_NMJ_Distance;

% Generate a column for NMJ List (default is empty)
nmjListPerElectrode = strings(size(electrodeX));

% Assign NMJs to each electrode
for i = 1:length(electrodeX)
    % Find NMJs belonging to this electrode
    matchIdx = find(nmjData.Electrode_X == electrodeX(i) & nmjData.Electrode_Y == electrodeY(i));
    
    if ~isempty(matchIdx)
        nmjListPerElectrode(i) = strjoin(string(nmjData.NMJ_ID(matchIdx)), ', '); % Join NMJ IDs
    else
        nmjListPerElectrode(i) = "None"; % No NMJs for this electrode
    end
end

% Generate Electrode ID ('Row_Col')
electrodeID = strings(size(electrodeX));
for i = 1:length(electrodeX)
    rowIdx = ceil(i / numElectrodes);
    colIdx = mod(i - 1, numElectrodes) + 1;
    electrodeID(i) = "R" + string(rowIdx) + "_C" + string(colIdx);
end

% Create Table
electrodeTable = table(electrodeID, electrodeX, electrodeY, heatmapActivation, correctedNMJDistances, nmjListPerElectrode, ...
                       'VariableNames', {'Electrode_ID', 'Electrode_X', 'Electrode_Y', ...
                                         'Heatmap_Activation', 'NMJ_Distance', 'NMJ_List'});

% Save CSV
csvFileName = 'Address/Electrode_Activation_NMJ_Distances_With_NMJs.csv';
writetable(electrodeTable, csvFileName);
disp(['Saved corrected electrode data with NMJ assignments to ', csvFileName]);