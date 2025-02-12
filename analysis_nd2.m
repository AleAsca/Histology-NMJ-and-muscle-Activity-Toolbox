%% analysis.m
%
% Created by: Alessandro Ascani Orsini
%
% Date: 2/15/2024
%
% Version: 1.0.0
%
% *Description*: Analyses the histology slides nd2 format and scales them
% properly, producing a 3D visualization of the slides stacked. It then
% plots and color codes the NMJ based on an excel sheet present in the same
% folder. 
%
% *REQUIREMENTS*: This code requires the nd2 file opener and Bresenham (for
% length estimation)
% - CY Y (2025). nd2read (https://github.com/joe-of-all-trades/nd2read), GitHub. Retrieved February 12, 2025.
% - Aaron Wetzler (2025). Bresenham optimized for Matlab (https://www.mathworks.com/matlabcentral/fileexchange/28190-bresenham-optimized-for-matlab), MATLAB Central File Exchange. Retrieved February 12, 2025.

%% Clean everything
close all;
clear all;
clc;

%% File location
dir_in = "Directory to your images"; % folder containing the images
excelN = "/NMJ.xlsx"; % name of the excel file. Add "/" before the name
dir_out = "Edited pics"; % if savePNG is on, this is the name of the folder that will contain the produced images

%% Colors Legends
col1 = [.25 0 0]; % color of NMJs type 1
col2 = [.5 0 0]; % color of NMJs type 2
col3 = [.75 0 0]; % color of NMJs type 3
col4 = [1 0 0]; % color of NMJs type 4
colnan = [0 1 1]; % color of NMJs type nan

%% Settings
printIm = false; % print the images once imported
savePNG = false; % save the images as png
slideS = .120; % spacing between slides in mm
scale = .10; % resolution of the 3D model from 0 to 1. 1 means that the resolution of the model is the same of the images. To avoid crash, keep low (reccomend 0.1)
umpx = 2.1538; % conversion from um to px might depend from pic to pic. Number of pixels per um

%% Import images
fprintf("IMPORTING IMAGES... \n");
files = struct2cell(dir(fullfile(dir_in, '*.nd2'))); % finds all the .nd2 files in the images folder
files = files(1, :); % keep only the names of the files
nim = size(files, 2); % number of images

% preallocate cell arrays for the images
imDouble = cell(nim, 1); % for storing double precision images
edit = cell(nim, 1); % cell array with edited images

%% Locate the NMJs based on excel file
sheetNames = sheetnames(dir_in+excelN); % get the list of sheets
dataExcel = cell(length(sheetNames), 1); % Preset cell array
scaledpx = dataExcel;
type = dataExcel;
for i = 1:length(sheetNames) % Loop through each sheet name and read the data into a table
    dataExcel{i} = readtable(dir_in+excelN, 'Sheet',sheetNames{i});
    if size(dataExcel{i},1)>0
        scaledpx{i} = [table2array(dataExcel{i}(:,2)),table2array(dataExcel{i}(:,3))].*umpx;
        type{i} = table2array(dataExcel{i}(:,1)); % extract the classifier
        if class(type{i}) == "cell" % check in case there is cell
            type{i}=str2double(type{i}); % set to nan since it's a character array, otherwise it would have converted
        end
    else
        type{i} = nan;
        scaledpx{i} =[nan,nan];
    end
end

%% Main loop
NMJs = edit;
for i = 1:nim % loop through the images
    fprintf("Histology slide: #%i \n",i);
    [a, ~] = nd2read(fullfile(dir_in, files{1, i})); % import the channels
    % (assuming nd2read function can handle the path and filename correctly).
    % If this returns Error, most likely nd2 file is corrupted (might happen when downloading)
    a = imadjust(a); % Enhance contrast
    if printIm
        figure;
        imshow(a);
    end
    imDouble = im2double(a); % Convert to double precision
    clear a; % these files are redundant and extremely heavy, removing them improves the performance of the code

    %% Save as png
    if savePNG
        if ~exist(dir_out, 'dir')
            mkdir(dir_out);
        end
        filePath = fullfile(dir_out, num2str(i)+".png");
        imwrite(imDouble,filePath);
    end

    %% Edit images - reduce their size and align
    fprintf("SELECT FEATURES... please press Enter to confirm. \n")
    edit{i} = imresize(imDouble, scale);
    clear imDouble;
    NMJs{i} = edit{i};
    NMJs{i}(:,:) = NaN;
    x = round(scaledpx{i}(:,1).*scale);
    y = round(scaledpx{i}(:,2).*scale);
    for j=1:length(x)
        if ~isnan(x(j)) && ~isnan(y(j))
            NMJs{i}(y(j),x(j)) = 1;
        end
    end
    [edit{i},~,shift] = removeBG(edit{i});
    if shift(2)<0 % Expand the image to prevent cropping
        NMJs{i} = [zeros(size(NMJs{i},1),-1.*shift(2)),NMJs{i}];
    else
        NMJs{i} = [NMJs{i},zeros(size(NMJs{i},1),shift(2))];
    end
    if shift(1)<0
        NMJs{i} = [zeros(-1.*shift(1),size(NMJs{i},2));NMJs{i}];
    else
        NMJs{i} = [NMJs{i};zeros(shift(1),size(NMJs{i},2))];
    end
    NMJs{i} = circshift(NMJs{i}, shift); % move the image
    figure;
    imshow(edit{i});
    hold on;
    %% Loop for feature selection until Enter is pressed
    j = 0;
    angle = [];
    while true
        j=j+1;
        [x, y] = ginput(2); % Select features
        if isempty(x) && j>1 % Break the loop if Enter is pressed
            theta = mean(angle);
            edit{i} = imrotate(edit{i}, theta, 'bilinear', 'loose'); % rotate the image as an average of the features
            %% rotate NMJs coordinates
            [yo,xo] = find(NMJs{i}==1); % get original coordinates of all the points
            oldC = size(NMJs{i})./2; %y and x of the original coordinates
            newC = size(edit{i})./2; %y and x of the new center coordinates
            dC = newC-oldC; % get the delta
            newXY = rotateCoordinates([xo,yo],oldC(end:-1:1),-theta); % rotate around the old origin
            NMJs{i} = [newXY(:,1)+dC(2),newXY(:,2)+dC(1)]; % translate to new origin and store
            figure; imshow(edit{i}); hold on; plot(NMJs{i}(:,1),NMJs{i}(:,2),'r.',MarkerSize = 40); % show the image
            break;
        end
        plot(x, y, 'r', 'LineWidth', 2); % Draw the diagonal
        dy = diff(y); % Change in y
        dx = diff(x); % Change in x
        angle = [angle; atan2d(dy, dx)]; % Angle in degrees
    end
    hold off
    type{i} = type{i}(1:size(NMJs{i},1));
    % Display the rotated image
    figure;
    imshow(edit{i});
    title('Rotated Image');
    close all;
end
fprintf("IMAGES SUCCESSFULLY ORIENTED. \n");
fprintf("IMPORT COMPLETED. \n");

%% Calculate volume of muscle
fprintf("ESTIMATING TOTAL VOLUME OF THE MUSCLE... \n");
tic
V = cellfun(@(x) nnz(x), edit)./(scale.^2); % number of pixel in each slide
V = V.*slideS.*1e3.*(1./umpx).^2.*1e-9; % volume of a slide in mm3
Vtot = sum(V); % total volume of the muscle in mm3
toc
fprintf("VOLUME MUSCLE: %f mm3 \n",Vtot);
save(dir_in+"/Vmuscle"+".mat","V","Vtot"); % save
fprintf("MUSCLE VOLUME SAVED TO FILE. \n");

%% Get the muscle longest distance
% To calculate the longest distance, we take the x axis and t
fprintf("DETERMINING MUSCLE SIZE... \n");
tic
d = @(a,b) sqrt(sum((a-b).^2,2)); % distance function
theta = 0:.1:180; % create rotation array
m = tand(theta)'; % get the slope
c = 0; % temporary intercept
sf = @(x) m.*x+c; % slope function
isf = @(y) (y-c)./m; % inverse slope function
refDist = 0; % preallocate reference distance
for i = 1:length(edit)
    kernel = ones(3, 3); % Define a 3x3 matrix of ones
    convResult = conv2(edit{i}>0.1, kernel, 'same'); % Apply convolution to get the border
    convResult = (convResult<9 & convResult>0); % get the contour
    [cY,cX] = size(edit{i});cX = cX./2; cY =cY./2; % calculate the center of the image
    c = cY-cX.*m; % get the intercept
    y1 = sf(1); % border 1
    y2 = sf(2.*cX); % border 2
    x1 = isf(1); % border3
    x2 = isf(2.*cY); % border4
    coo = [ones(length(m),1),y1.*((y1 >= 1)&(y1 <= 2.*cY))]; % coordinates of the borders based on the angle
    coo = [coo,ones(length(m),1).*2.*cX,y2.*((y2 >= 1)&(y2 <= 2.*cY))];
    coo = [coo,x1.*((x1 >= 1)&(x1 <= 2.*cX)),ones(length(m),1)];
    coo = [coo,x2.*((x2 >= 1)&(x2 <= 2.*cX)),ones(length(m),1).*2.*cX];
    for j = 1:size(coo,1)
        pt = zeros(1,4);
        num = 1;
        if coo(j,2)~=0
            pt(num) = coo(j,1);
            pt(num+1) = coo(j,2);
            num = num+2;
        end
        if coo(j,4)~=0
            pt(num) = coo(j,3);
            pt(num+1) = coo(j,4);
            num = num+2;
        end
        if coo(j,5)~=0
            pt(num) = coo(j,5);
            pt(num+1) = coo(j,6);
            num = num+2;
        end
        if coo(j,7)~=0
            pt(num) = coo(j,7);
            pt(num+1) = coo(j,8);
        end
        if max(isnan(pt))==0
            [xB,yB] = bresenham(pt(1),pt(2),pt(3),pt(4));
            xB = xB(xB>0 & xB<2.*cX & yB>0 & yB<2.*cY);
            yB = yB(yB>0 & yB<2.*cY);
            dx = ((min(xB)+max(xB))./2-cX);
            dy = ((min(yB)+max(yB))./2-cY);
            xB = xB-dx;
            yB = yB-dy; % adjust the y
            for k=1:length(xB)
                diag(k)=convResult(round(yB(k)),round(xB(k)));
            end
            [~,allIdx]=find(diag>0); %get the idx of the max
            maxIdx = max(allIdx); % get the
            minIdx = min(allIdx);
            dist=abs(maxIdx-minIdx); % get the distance
            if dist>refDist
                refDist = dist;
            end
        end
    end
end
toc
refDist = refDist./(umpx.*scale).*1e-3; % length of the muscle converted in mm
fprintf("MUSCLE SIZE SAVED. Muscle: %f mm \n", refDist);

%% Rebuild in 3D
fprintf("PRODUCING 3D RECONSTRUCTION... \n");
%px = max(cellfun(@(x) length(x),imDouble)); % get the maximum pixel count on the longest side of a slide
pxH = round(slideS.*1e3.*umpx).*scale; % calculate the amount of pixels to scale the image properly
fig = figure;
hold on
for i = 1:nim
    [h, w, ~] = size(edit{i});
    [X, Y] = meshgrid(linspace(-w/2, w/2, w), linspace(-h/2, h/2, h)); % centering images
    Z = -pxH.*i.*ones(size(X)); % Different depths to give 3d size
    CDataCopy = edit{i}.*(edit{i}>0); % brightness based filter - it removes parts of the muscle below a certain threshold. can be useful
    CDataCopy(CDataCopy == 0) = NaN; % Replace zeros with NaN
    surf_h = surface('XData', X, 'YData', Y, 'ZData', Z, 'CData', CDataCopy, ...
        'EdgeColor', 'none');
    alpha scaled
    alpha(0.3); % add transparency
    %% NMJs plotting
    z = (NMJs{i}(:,1)./NMJs{i}(:,1)).*Z(1);
    NMJs{i}(:,3) = z; % save z coordinates
    c1 = find(type{i}==1);
    c2 = find(type{i}==2);
    c3 = find(type{i}==3);
    c4 = find(type{i}==4);
    cnan = find(isnan(type{i}));
    plot3((NMJs{i}(c1,1)-w./2),(NMJs{i}(c1,2)-h./2),z(c1), '.', Color=col1, MarkerSize =30); % plot the nmjs
    plot3((NMJs{i}(c2,1)-w./2),(NMJs{i}(c2,2)-h./2),z(c2), '.', Color=col2, MarkerSize =30);
    plot3((NMJs{i}(c3,1)-w./2),(NMJs{i}(c3,2)-h./2),z(c3), '.', Color=col3, MarkerSize =30);
    plot3((NMJs{i}(c4,1)-w./2),(NMJs{i}(c4,2)-h./2),z(c4), '.', Color=col4, MarkerSize =30);
    plot3((NMJs{i}(cnan,1)-w./2),(NMJs{i}(cnan,2)-h./2),z(cnan), '.', Color=colnan, MarkerSize =30);
end
% Make dummy plots for the legend
h1 = plot3(NaN, NaN, NaN, '.', 'Color', col1, 'MarkerSize', 30);
h2 = plot3(NaN, NaN, NaN, '.', 'Color', col2, 'MarkerSize', 30);
h3 = plot3(NaN, NaN, NaN, '.', 'Color', col3, 'MarkerSize', 30);
h4 = plot3(NaN, NaN, NaN, '.', 'Color', col4, 'MarkerSize', 30);
h5 = plot3(NaN, NaN, NaN, '.', 'Color', colnan, 'MarkerSize', 30);
legend([h1, h2, h3, h4, h5], {'Type 1', 'Type 2', 'Type 3', 'Type 4', 'Type NaN'}); % add legend
view(3);
colormap(gray); % bone is also pretty cool
axis equal off;
grid on;
hold off;
fprintf("RECONSTRUCTION COMPLETED! \n");

%% Save the image
fprintf("SAVING MODEL...\n");
set(gcf, 'Position', get(0, 'ScreenSize'));
figName = strsplit(dir_in,"/");
savefig(fig, dir_in+'/'+figName(end)+'.fig'); % Save the figure as fig
exportgraphics(fig, dir_in+'/'+figName(end)+'.png'); % Save the figure as PNG

%% Save the coordinates of the NMJs
fprintf("SAVING COORDINATES OF NMJs... \n");
save(dir_in+"/NMJsCoordinates"+".mat","NMJs","type","refDist");
fprintf("DONE! \n");

%% rotateCoordinates function
% Description: This function rotates an array of coordinates
% [x,y;x1,y1;...] around a pivot by an angle in a clockwise manner.
function nC = rotateCoordinates(C, pivot, angle)
rM = [cosd(angle) -sind(angle); sind(angle) cosd(angle)]; % rotation matrix
Cs = C - pivot; % shifted coordinates to origin
Cr = (rM * Cs')'; % rotate and transpose back
nC = Cr + pivot; % shift coordinates back to original pivot
end

%% removeBG function
% Description: This function takes a greyscale image, detects the largest
% connected area in it, and applies a mask, returning the isolated area.
% It's ideal to be used for slice detection in histology. It also centers
% the slice
function [cleanimg, mask, shift] = removeBG(img_gray)
%img_gray = rgb2gray(img_gray);
threshold = graythresh(img_gray);
img_binary = imbinarize(img_gray, threshold); % Threshold the image to create a binary image
img_filled = imfill(img_binary, 'holes'); % Fill holes in the binary image
se = strel('disk', 10); % create a structuring element 10 pixels in radius
img_dilated = imdilate(img_filled, se); % Dilate the image slightly to connect components
cc = bwconncomp(img_dilated); % Find connected components (ideally bg)
stats = regionprops(cc, 'Area', 'PixelIdxList');
[~, largestIdx] = max([stats.Area]); % Find the largest connected component
mask = false(size(img_gray)); % Preallocate mask
mask(stats(largestIdx).PixelIdxList) = true; % Change central region to 1
% Calculate the centroid of the mask
mstat = regionprops(mask,'Centroid');
centroid = mstat.Centroid;
shift = round(size(img_gray)/2 - centroid(2:-1:1)); % Calculate the shift required to center the centroid
cleanimg = zeros(size(img_gray), 'like', img_gray);
%cleanimg(:) = NaN; % Uncomment if you prefer NaN
cleanimg(mask) = img_gray(mask);
if shift(2)<0 % Expand the image to prevent cropping
    cleanimg = [zeros(size(cleanimg,1),-1.*shift(2)),cleanimg];
else
    cleanimg = [cleanimg,zeros(size(cleanimg,1),shift(2))];
end
if shift(1)<0
    cleanimg = [zeros(-1.*shift(1),size(cleanimg,2));cleanimg];
else
    cleanimg = [cleanimg;zeros(shift(1),size(cleanimg,2))];
end
cleanimg = circshift(cleanimg, shift);
end

%% Scale images
% Disabled because the images are already properly scaled
%
% fprintf("SCALING IMAGES... \n");
% pdist = zeros(nim,1);
% profS = cell(nim,1);
% % determine spacing
% for i=1:nim
%     bE = im{i}(end-300:end-60,:); % Get the bottom profile
%     prof = sum(bE,1); % Get the profile by summing the pixel values vertically
%     prof = movmean(prof, 100); % filters through moving average
%     prof = prof-min(prof(200:end-200)); % make sure it is at 0.
%     prof = prof.*(prof>1.5e6); % adjust to better determine the peaks
%     profS{i} = prof;
%     [peak,loc] = findpeaks(prof); % extract location of the peaks
%     dist = diff(loc); % calculate the difference
%     dist = dist(3:end-3);
%     pt = findpeaks(dist.*(dist>650)); % calculate the distance between the important peaks
%     pdist(i) = mean(pt(1:3)); % take first peak
% end
% sf = pdist./max(pdist); % determine scaling factor for each
% fprintf('Scaling factor: x %f \n', sf); % show the scaling factors
% for i = 1:nim
%     edit{i} = imresize(im{i}, sf(i));
%     if printIm
%         figure;
%         imshow(edit{i});
%     end
% end
