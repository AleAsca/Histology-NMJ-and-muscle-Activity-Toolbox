function [points, w, l] = extract3DPointsFromFig(filename)
    %% extract3DPointsFromFig.m
    %
    % Created by: Alessandro Ascani Orsini
    %
    % Date: 4/16/2024
    %
    % Version: 1.0.0
    %
    % *Description*: This function extracts the NMJ points from a 3D muscle
    % reconstruction as well as its dimensions (width along the y axis and
    % length along the x axis).

    %% Settings
    fig = openfig(filename, 'invisible');  % Open the figure invisibly
    axesHandles = findall(fig, 'type', 'axes'); % Access the axes in the figure
    % Preallocate data
    points = [];
    w = 0;
    l = 0;
    maxNonNanHeight = 0;  % Initialize max non-NaN height

    %% Iterate over all axes
    for ax = axesHandles'
        % Get all children of type 'scatter' or 'line'
        scatterObjects = findall(ax, 'Type', 'scatter');
        lineObjects = findall(ax, 'Type', 'line');
        surfObjects = findall(ax, 'Type', 'surface');
        
        % Extract points from scatter objects
        for scatterObj = scatterObjects'
            xData = get(scatterObj, 'XData');
            yData = get(scatterObj, 'YData');
            zData = get(scatterObj, 'ZData');
            points = [points; [xData(:), yData(:), zData(:)]];
        end
        
        % Extract points from line objects
        for lineObj = lineObjects'
            xData = get(lineObj, 'XData');
            yData = get(lineObj, 'YData');
            zData = get(lineObj, 'ZData');
            points = [points; [xData(:), yData(:), zData(:)]];
        end
        
        % Detect width and length of muscle
        for surfObj = surfObjects'
            cData = get(surfObj,'CData');
            validIndices = ~isnan(cData);
            
            tempw = max(sum(validIndices,1));
            templ = max(sum(validIndices,2));
            if tempw>w
                w=tempw;
            end
            if templ>l
                l=templ;
            end
        end
    end
    points = points(~isnan(points(:,1)), :); % Remove NaN entries from points
    close(fig); % Close the figure
end
