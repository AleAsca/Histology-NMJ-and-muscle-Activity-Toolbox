function points = Muscle3DModel(filename,MeshC,MeshA)
    %% Muscle3DModel.m
    %
    % Created by: Alessandro Ascani Orsini
    %
    % Date: 4/17/2024
    %
    % Version: 0.0.1
    %
    % *Description*: This function builds a 3D model of the muscle

    %% Settings
    fig = openfig(filename);  % Open the figure
    axesHandles = findall(fig, 'type', 'axes'); % Access the axes in the figure
    % Preallocate data
    points = [];

    %% Correct missing variables
    if ~exist('MeshC','var') || max(size(MeshC)~=[1,3])
        MeshC = [0,0,0];
    end
    if ~exist('MeshA','var')
        MeshA = 0.05;
    end

    %% Iterate over all axes
    for ax = axesHandles'
        surfObjects = findall(ax, 'Type', 'surface'); % get surf objects
        for surfObj = surfObjects'
            cData = get(surfObj,'CData'); % get the CData
            xData = get(surfObj, 'XData');
            yData = get(surfObj, 'YData');
            zData = get(surfObj, 'ZData');
            contourn = edge(~isnan(cData),'Sobel'); % get the contourn
            points = [points;
                        xData(contourn),yData(contourn),zData(contourn)];
        end
    end
    points = points(~isnan(points(:,1)), :); % Remove NaN entries from points
    K = convhull(points(:,1),points(:,2),points(:,3)); % get the convolutional hull
    hold on
    freezeColors;
    trisurf(K, points(:,1), points(:,2), points(:,3), ...
                    'FaceColor', MeshC, 'FaceAlpha', MeshA, 'EdgeColor', 'none','HandleVisibility', 'off');
    hold off
    %close(fig); % Close the figure
end
