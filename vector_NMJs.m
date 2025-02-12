function [pv, ZXangle] = vector_NMJs(filename,vectorScale)
    %% vector_NMJs.m
    %
    % Created by: Alessandro Ascani Orsini
    %
    % Date: 1/30/2025
    %
    % Version: 0.0.1
    %
    % *Description*: This function plots the PCA vector of the points in
    % fig filename. It then calculates the angle on the ZX plane.

    %% Preset vectorScale
    if ~exist('vectorScale','var')|| isempty(vectorScale)
        vectorScale = 200; % set vectorScale to 200
    end

    %% Extract the points from the figure and center them
    [points,~,~] = extract3DPointsFromFig(filename);
    mNMJ = mean(points); % find the middle
    cp = points - mNMJ; % center them

    %% Perform PCA
    cov_matrix = cov(cp); % compute the covariance matrix
    [eigenvectors, eigenvalues] = eig(cov_matrix); % compute eigenvectors and eigenvalues
    [~, idx] = max(diag(eigenvalues)); % get the largest eigenvector
    pv = eigenvectors(:, idx); % direction of the vector

    %% Calculate angle of the NMJs
    ZXangle = atan(pv(3)./pv(1)) % angle of the plane

    %% Plot
    fig = openfig(filename);
    hold on;
    quiver3(mNMJ(1), mNMJ(2), mNMJ(3), ...
            pv(1)*vectorScale, pv(2)*vectorScale, pv(3)*vectorScale, ...
            'k', 'LineWidth', 2, 'MaxHeadSize', 1);
    grid on;
end