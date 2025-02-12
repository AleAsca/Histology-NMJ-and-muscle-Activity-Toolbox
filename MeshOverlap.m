%% MeshOverlap.m
%
% Created by: Alessandro Ascani Orsini
%
% Date: 4/26/2024
%
% Version: 1.0.0
%
% *Description*: Takes in the image of the muscle and taken the manually
% calculated ratio from the image it overlaps with the NMJ locations. This
% assumes that the direction of slicing is the same along which the
% electrodes are placed

%% Clean everything
close all;
clear all;
clc;

%% Files location
figfile = "Directory of the .fig file of the muscle reconstruction"; % reference

%% Settings
umpx = .21538; % scaled pixels/um
dT = 0.3279; % distance from the top normalized
sR = 0.3443; % size ratio of electrodes compared to muscle
topLR = 1; % indicate position of the top of the muscle. Right = 1, Left = 0
damp = -.015; % damping factor for tissue
showIMG = false; % set to true to view each layer
hframe = 50; % frame of the heatmap to plot (80 was the one used to calibrate)
heatmapSample = "Directory of the heatmap as a 2D .mat file"; % directory to sample heatmap

%% Load
tic
[NMJs,w,l,muscle] = extract3DPointsFromFig(figfile); % load all the coordinates
toc
load(heatmapSample);

%% Recreate hull in the electrode region
% the mask needs to be adjusted so that it's considering a subtraction from
% the outer region of the electrodes, rather than adding it full always
% (causes a shift from the center)
lmax = max(l);
LdT = lmax.*dT; % scaled distance from top
LsR = lmax.*sR; % scaled length of the electrodes
cL = lmax./2; % center of the length to calibrate
lastR = round((1-dT-sR).*lmax); % remaining of the length
% define distances from center
distC1 = cL-LdT;
distC2 = cL-lastR;
for j = size(muscle,1):-1:1
    mask{j} = NaN(size(muscle{j,1})); % preallocate
    cM{j} = size(mask{j})./2;% center of the image
    cLocal = cM{j}(2); % get the x axis
    if topLR==1
        ll = cLocal-distC1; % lower limit
        ul = cLocal+distC2; % upper limit
    else
        ll = cLocal-distC2;
        ul = cLocal+distC1;
    end
    mask{j}(:,round(ll.*(ll>0)+(1.*(ll<1)):(ul.*(ul<=cLocal.*2)+cLocal.*2.*(1-(ul<=cLocal.*2))))) = 1; % highlight the region between upper and lower limit fixing it between 0 and the maximum center
    mask{j} = (~isnan(muscle{j,1}.*mask{j}(1:size(muscle{j,1},1),1:size(muscle{j,1},2)))>0); % fix the mask
end

%% Overlay electrodes
clear testH;
testH(:,:) = heatmap(hframe,:,:);
imagesc(testH);
testH = rot90(imresize(testH,umpx./0.1),3); % sample heatmap timeframe scaled (100 px = 1mm heatmap --> 1 px = 10µm ) and rotated
for j = 1:size(muscle,1)
    eVal{j} = zeros(size(mask{j}));
    %eVal{i,j} = NaN(size(mask{i,j}));
    rX = round(cM{j}(2)-size(testH,2)./2):round(cM{j}(2)+size(testH,2)./2)-1; % range X
    rX = rX(rX>0 & rX<=size(eVal{j},2));
    rY = round(cM{j}(1)-size(testH,1)./2):round(cM{j}(1)+size(testH,1)./2)-1; % range Y
    rY = rY(rY>0 & rY<=size(eVal{j},1));
    eVal{j}(rY,rX) = testH(1:length(rY),1:length(rX));
    EMGslide{j} = eVal{j}.*mask{j};
    %EMGslide{i,j}(EMGslide{i,j} == 0) = NaN;
    if showIMG
        figure;
        imagesc(EMGslide{j});
    end
end

%% Now decrease the intensity of already drawn areas
for j = 1:size(muscle,1)
    if j>1
        j = size(muscle,1)+1-j; % adjust so that it's in inverse order
        EMGslidediff{j} = centeredSubtraction(EMGslide{j},cumEMGdiff{j+1}.*damp); % subtract the difference
        cumEMGdiff{j} = centeredSubtraction(cumEMGdiff{j+1},-EMGslidediff{j}); % make the difference cumulative
        if showIMG
            figure;
            imagesc(EMGslidediff{j});
        end
    else
        j = size(muscle,1)+1-j;
        EMGslidediff{j} = EMGslide{j};
        cumEMGdiff{j} = EMGslide{j};
    end
end

%% Plot 3D
figure
hold on
for j = 1:size(muscle,1)
    EMGslidediff{j}(EMGslidediff{j}==0) = nan;
    [h, w, ~] = size(EMGslidediff{j});
    [Xd, Yd] = meshgrid(linspace(-w./2, w./2, w), linspace(-h./2, h./2, h));
    Zd = muscle{j,4}(1).*ones(size(EMGslidediff{j}));
    surf_h = surf('XData', Xd,'YData', Yd, 'ZData', Zd, ...
        'CData', EMGslidediff{j},'EdgeColor', 'none');
    %alpha(surf_h,0.2);
end
plot3(NMJs(:,1),NMJs(:,2),NMJs(:,3),'.r',MarkerSize=20);
hold off
axis equal off
colorbar
view(3)
clim([0,100]);

%% Subtract two matrices from the center
function C = centeredSubtraction(A, B)
[sizeA1, sizeA2] = size(A); % Get the sizes of A and B
[sizeB1, sizeB2] = size(B);
sizeC1 = max(sizeA1, sizeB1); % Determine the dimensions of the resulting matrix
sizeC2 = max(sizeA2, sizeB2);
C_A = zeros(sizeC1, sizeC2); % Initialize matrices C_A and C_B with zeros of dimensions sizeC1 x sizeC2
C_B = zeros(sizeC1, sizeC2);
startA1 = floor((sizeC1 - sizeA1) / 2) + 1; % Calculate the starting indices for A in C_A
startA2 = floor((sizeC2 - sizeA2) / 2) + 1;
C_A(startA1:startA1+sizeA1-1, startA2:startA2+sizeA2-1) = A;% Place A in the center of C_A
startB1 = floor((sizeC1 - sizeB1) / 2) + 1;  % Calculate the starting indices for B in C_B
startB2 = floor((sizeC2 - sizeB2) / 2) + 1;
C_B(startB1:startB1+sizeB1-1, startB2:startB2+sizeB2-1) = B; % Place B in the center of C_B
C = C_A - C_B; % Subtract C_B from C_A and clamp values at 0
C(C<1) = 0;
end
