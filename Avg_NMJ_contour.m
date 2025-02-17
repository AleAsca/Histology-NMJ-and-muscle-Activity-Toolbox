%% Avg_NMJ_contour.m
%
% Created by: Alessandro Ascani Orsini
%
% Date: 4/13/2024
%
% Version: 0.0.1
%
% *Description*: Takes in the rotate coordinates of NMJs and averages them
% with standard deviation using a contour method, this means it looks for
% the densest point in the plot and gets the closest 90% of the points.
% Because of how the code is structured, it doesn't take in consideration
% the points from dir_in{1} when counting and averaging
%
% *REQUIREMENTS*: to run this code you need the muscle reconstructions to
% be in a .fig file format. You will also need the vol3D program for 3D
% heatmaps.
% - Oliver Woodford (2025). vol3d v2 (https://www.mathworks.com/matlabcentral/fileexchange/22940-vol3d-v2), MATLAB Central File Exchange. Retrieved February 12, 2025.

%% Clean everything
close all;
clear all;
clc;

%% Files location
figfile{1} = "Muscle Model Directory (.fig file)"; % Model muscle
figfile{2} = "First Directory to your first reconstructed muscle (.fig file)"; % data will actually be looked up from this file and following ones
figfile{3} = "Second Directory to your second reconstructed muscle (.fig file)"; % and so on.

%% Settings
dir_model = "directory";
fname = "NMJsCoordinates.mat";
scalePoints = 10; % scale factor for size of markers
umpx = 2.1538; % pixels/um
size_vox = 20.*umpx; % size of voxels
rnp = 60; % radius of neighbouring points
minp = 4; % minimum number of points to have a cluster

%% Load
tic
[NMJs,w,l] = cellfun(@(x) extract3DPointsFromFig(x),figfile,'UniformOutput', false); % load all the coordinates
toc
%otherData = cellfun(@(x) load(x+"/"+fname),dir_in); % load length of muscle and type
load(dir_model); % load model

%% Extract the data points
allPoints = NMJs;
refMin = cellfun(@(x) floor(min(x, [], 1)),allPoints,'UniformOutput',false); % get the reference max of each
refMax = cellfun(@(x) ceil(max(x, [], 1)),allPoints,'UniformOutput',false);
refMin = vertcat(refMin{:});
refMax = vertcat(refMax{:});
refLen = l{1};
refWid = w{1};

%% Normalize the data points and exclude dir{1}
normPoints = allPoints; % preallocate array
onModel = {}; % Initialize empty, will concatenate later without dir{1}
refRange = refMax - refMin; % get the range
for i = 1:size(allPoints,2)
    normPoints{i} = allPoints{i}./refRange(i,:); % normalize the points
    tempOnModel = (normPoints{i}); 
    tempOnModel(:,1) = (tempOnModel(:,1).*refRange(i,1).*refLen./l{i}); % scale length
    tempOnModel(:,2) = (tempOnModel(:,2).*refRange(i,2).*refWid./w{i}); % scale width
    tempOnModel(:,3) = tempOnModel(:,3).*refRange(1,3); % clamp the z
    if i ~= 1 % Skip points from dir{1}
        onModel{end+1} = tempOnModel; % Add to onModel if not from dir{1}
    end
end
onModelCell = onModel; % save the cell method
onModel = vertcat(onModel{:}); % Concatenate all points excluding dir{1}
onModel(:,1) = onModel(:,1);
onModel(:,2) = onModel(:,2);

%% Define grid
x = linspace(min(onModel(:,1)), max(onModel(:,1)), ceil((max(onModel(:,1)) - min(onModel(:,1)))/size_vox) + 1);
y = linspace(min(onModel(:,2)), max(onModel(:,2)), ceil((max(onModel(:,2)) - min(onModel(:,2)))/size_vox) + 1);
z = linspace(min(onModel(:,3)), max(onModel(:,3)), ceil((max(onModel(:,3)) - min(onModel(:,3)))/size_vox) + 1);

%% check which points fall in each square
figure
[cx,~,ix] = histcounts(onModel(:,1), x);
[cy,~,iy] = histcounts(onModel(:,2), y);
[cz,~,iz] = histcounts(onModel(:,3), z);

a = zeros(length(y),length(x),length(z)); % matrix to keep track of the count
for i = 1:size(onModel, 1)
    a(iy(i), ix(i), iz(i)) = a(iy(i), ix(i), iz(i)) + 1; % Increment the voxel count by 1 for the voxel that this point falls into
end
%a(a<=1) = 0;
xlm = x; % calculate the x limits
ylm = y; % calculate the y limits
vol3d(CData=a./(length(figfile)-1),XData=[xlm(1),xlm(end)+size_vox],YData=[ylm(1),ylm(end)+size_vox],ZData=[z(1),z(end)+size_vox]);
colorbar
colormap hot;
axis equal

%% Get the average coordinate in each vox
[C,ia,ic] = unique([ix,iy,iz],"rows");
avgC = [accumarray(ic, onModel(:,1), [], @mean),accumarray(ic, onModel(:,2), [], @mean),accumarray(ic, onModel(:,3), [], @mean)];
% diffC = avgC-onModel(ia,:); % remove the ones where there is less than 2 coordinates
% avgC = avgC((diffC(:,1)~=0),:);
hold on; 
plot3(avgC(:,1),avgC(:,2),avgC(:,3),'b.',MarkerSize = 10)

%% plot the model
pxH = z(end):(z(1)-z(end))/size(edit,1):z(1);
for i = 1:size(edit,1)
    [h, w, ~] = size(edit{i});
    [X, Y] = meshgrid(linspace(-w./2, w./2, w), linspace(-h./2, h./2, h)); % centering images
    Z = pxH(i).*ones(size(X)); % Different depths to give 3d size
    CDataCopy = edit{i};
    CDataCopy(CDataCopy == 0) = NaN; % Replace zeros with NaN
    surf_h = surface('XData', X, 'YData', Y, 'ZData', Z, 'CData', CDataCopy./(length(figfile)-1), ...
        'EdgeColor', 'none');
    alpha(surf_h,0.05); % add transparency
end
axis off equal;
view(3);
hold off

%% Define contour
figure;
sgtitle("Overlap of individual clusters")
hold on;
pxH = max(onModel(:,3)):(min(onModel(:,3))-max(onModel(:,3)))/size(edit,1):min(onModel(:,3));
for i = 1:size(edit,1)
    [h, w, ~] = size(edit{i});
    [X, Y] = meshgrid(linspace(-w./2, w./2, w), linspace(-h./2, h./2, h));
    Z = pxH(i).*ones(size(X));
    CDataCopy = edit{i};
    CDataCopy(CDataCopy == 0) = NaN;
    surf_h = surface('XData', X, 'YData', Y, 'ZData', Z, 'CData', CDataCopy, ...
        'EdgeColor', 'none');
    alpha(surf_h,0.05);
end

for i = 2:size(allPoints,2)
    % Iterate through each file and each cluster and plot
    IDX = dbscan(onModelCell{i-1}, rnp, minp);
    denseCluster = mode(IDX(IDX~=-1));
    selected_points = onModelCell{i-1}(IDX == denseCluster, :);
    k=0;
    unique_clusters = unique(IDX);
    itCol = rand(1,3);
    for cluster_id = unique_clusters'
        if cluster_id ~= -1 % Ignore noise points
            cluster_points = onModelCell{i-1}(IDX == cluster_id, :);
            if size(cluster_points,1)>4
                if size(unique(cluster_points),2)==1 % if the points are coplanar
                    % Add a small perturbation to one of the points to break coplanarity
                    cluster_points(1,3) = cluster_points(1,3) + 1; % adjust z-coordinate slightly
                end
                k=k+1;
                K{k,i} = convhull(cluster_points(:,1), cluster_points(:,2), cluster_points(:,3)); % recompute convex hull and save
                trisurf(K{k,i}, cluster_points(:,1), cluster_points(:,2), cluster_points(:,3), ...
                    'FaceColor', itCol, 'FaceAlpha', 0.5); % Random color for each cluster
            end
        end
    end
end

plot3(onModel(:,1), onModel(:,2), onModel(:,3), '.b', 'MarkerSize', 20);
axis equal off;
view(2);
colorbar;
colormap hot;
clim([0,1]);
hold off;

%% Define contour overall
figure;
sgtitle("Overall clusters");
hold on;
pxH = max(onModel(:,3)):(min(onModel(:,3))-max(onModel(:,3)))/size(edit,1):min(onModel(:,3));
for i = 1:size(edit,1)
    [h, w, ~] = size(edit{i});
    [X, Y] = meshgrid(linspace(-w./2, w./2, w), linspace(-h./2, h./2, h));
    Z = pxH(i).*ones(size(X));
    CDataCopy = edit{i};
    CDataCopy(CDataCopy == 0) = NaN;
    surf_h = surface('XData', X, 'YData', Y, 'ZData', Z, 'CData', CDataCopy, ...
        'EdgeColor', 'none');
    alpha(surf_h,0.05);
end

IDX = dbscan(onModel, rnp, minp);
denseCluster = mode(IDX(IDX~=-1));
selected_points = onModel(IDX == denseCluster, :);
k=0;
unique_clusters = unique(IDX);
itCol = rand(1,3);
for cluster_id = unique_clusters'
    if cluster_id ~= -1 % Ignore noise points
        cluster_points = onModel(IDX == cluster_id, :);
        if size(cluster_points,1)>4
            if size(unique(cluster_points),2)==1 % if the points are coplanar
                % Add a small perturbation to one of the points to break coplanarity
                cluster_points(1,3) = cluster_points(1,3) + 1; % adjust z-coordinate slightly
            end
            k=k+1;
            K{k,i} = convhull(cluster_points(:,1), cluster_points(:,2), cluster_points(:,3)); % recompute convex hull and save
            trisurf(K{k,i}, cluster_points(:,1), cluster_points(:,2), cluster_points(:,3), ...
                'FaceColor', itCol, 'FaceAlpha', 0.5); % Random color for each cluster
        end
    end
end
plot3(onModel(:,1), onModel(:,2), onModel(:,3), '.b', 'MarkerSize', 20);
axis equal off;
view(2);
colorbar;
colormap hot;
clim([0,1]);
hold off;
