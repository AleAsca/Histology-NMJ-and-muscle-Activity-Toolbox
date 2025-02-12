function points = ShowSide(filename)
    %% ShowSide.m
    %
    % Created by: Alessandro Ascani Orsini
    %
    % Date: 7/15/2024
    %
    % Version: 0.0.1
    %
    % *Description*: This function shows the side of each slide

    %% Settings
    fig = openfig(filename);  % Open the figure
    axesHandles = findall(fig, 'type', 'axes'); % Access the axes in the figure
    % Preallocate data
    points = [];
    
    %% Iterate over all axes
    hold on
    for ax = axesHandles'
        surfObjects = findall(ax, 'Type', 'surface'); % get surf objects
        for surfObj = surfObjects'
            cData = get(surfObj,'CData'); % get the CData
            xData = get(surfObj, 'XData');
            yData = get(surfObj, 'YData');
            zData = get(surfObj, 'ZData');
            contourn = edge(~isnan(cData),'Sobel'); % get the contourn
            points = [xData(contourn),yData(contourn),zData(contourn)];
            plot3(points(:,1), points(:,2), points(:,3), 'k-', 'LineWidth', 2);
        end
    end
    hold off
end
