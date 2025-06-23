% =========================================================================
% Script Name: Setregions.m
%
% Description:
%   This program allows the user to manually select four square regions 
%   corresponding to the polarized projections of bead images.
%   These selected regions are used for correcting optical distortions 
%   in the 4polar3D analysis workflow.
%
%   The user is guided step-by-step to select:
%   - The first rectangle (left-bottom)
%   - The second rectangle (right-bottom)
%   - The third rectangle (left-top)
%   - The fourth rectangle (right-top)
%
%   The selected regions are cropped, optionally padded, and filled with 
%   random noise in the background to prevent zero-valued pixels.
%   Finally, the selected regions are saved into a .mat file.
%
% Authors:
%   Cesar Valades-Cruz        - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

%% Initialize
clear all
close all
clc

%% Select TIFF file
[FileName, PathName, FilterIndex] = uigetfile('D:\4POLARSTORM\_RAW DATA\*.tif');
tiff_file = [PathName '\' FileName];

%% Load the first image
data1 = double(imread(tiff_file, 1));

%% Display image for manual region selection
h1 = figure(1);
imagesc(data1)
colormap('gray')
axis('image')
contrasth = imcontrast(h1); % Interactive contrast adjustment

%% Select the first rectangle (left-bottom)
h = imrect;
verticesleft = wait(h);
BW2l = createMask(h);
rectveL = round(getPosition(h));

%% Select the second rectangle (right-bottom)
msgbox('Select the second rectangle (right-bottom)', 'Hi', 'modal')
verticesright = wait(h);
BW2r = createMask(h);
rectveR = round(getPosition(h));

%% Select the third rectangle (left-top)
msgbox('Select the third rectangle (Left-top)', 'Hi', 'modal')
verticesleftup = wait(h);
BW2lu = createMask(h);
rectveLU = round(getPosition(h));

%% Select the fourth rectangle (right-top)
msgbox('Select the fourth rectangle (Right-top)', 'Hi', 'modal')
verticesrightup = wait(h);
BW2ru = createMask(h);
rectveRU = round(getPosition(h));

%% Apply masks to isolate each selected region
data1sel = data1 .* BW2l;
data2sel = data1 .* BW2r;
data3sel = data1 .* BW2lu;
data4sel = data1 .* BW2ru;

%% Crop selected regions
data1selcrop = data1sel(rectveL(2):(rectveL(2)+rectveL(4)), rectveL(1):(rectveL(1)+rectveL(3)), :);
data2selcrop = data2sel(rectveR(2):(rectveR(2)+rectveR(4)), rectveR(1):(rectveR(1)+rectveR(3)), :);
data3selcrop = data3sel(rectveLU(2):(rectveLU(2)+rectveLU(4)), rectveLU(1):(rectveLU(1)+rectveLU(3)), :);
data4selcrop = data4sel(rectveRU(2):(rectveRU(2)+rectveRU(4)), rectveRU(1):(rectveRU(1)+rectveRU(3)), :);

%% Determine maximum region size to standardize all crops
maxM = max([size(data1selcrop,1), size(data2selcrop,1), size(data3selcrop,1), size(data4selcrop,1)]);
maxN = max([size(data1selcrop,2), size(data2selcrop,2), size(data3selcrop,2), size(data4selcrop,2)]);

%% Optional contour (set to zero by default)
contourIm = 0;

%% Initialize padded arrays for each region
data1selcropnew = zeros(maxM + 2*contourIm, maxN + 2*contourIm, size(data1selcrop,3));
data2selcropnew = zeros(maxM + 2*contourIm, maxN + 2*contourIm, size(data1selcrop,3));
data3selcropnew = zeros(maxM + 2*contourIm, maxN + 2*contourIm, size(data1selcrop,3));
data4selcropnew = zeros(maxM + 2*contourIm, maxN + 2*contourIm, size(data1selcrop,3));

%% Shift cropped regions into new padded arrays
data1selcrop = imtranslate(data1selcrop, [contourIm, contourIm], 'OutputView', 'full');
data1selcropnew(1:size(data1selcrop,1), 1:size(data1selcrop,2), :) = data1selcrop;

data2selcrop = imtranslate(data2selcrop, [contourIm, contourIm], 'OutputView', 'full');
data2selcropnew(1:size(data2selcrop,1), 1:size(data2selcrop,2), :) = data2selcrop;

data3selcrop = imtranslate(data3selcrop, [contourIm, contourIm], 'OutputView', 'full');
data3selcropnew(1:size(data3selcrop,1), 1:size(data3selcrop,2), :) = data3selcrop;

data4selcrop = imtranslate(data4selcrop, [contourIm, contourIm], 'OutputView', 'full');
data4selcropnew(1:size(data4selcrop,1), 1:size(data4selcrop,2), :) = data4selcrop;

%% Replace zero-valued pixels with random noise to avoid artifacts
randmatrix = rand(size(data1selcropnew));
data1selcropnew(data1selcropnew == 0) = randmatrix(data1selcropnew == 0);

randmatrix = rand(size(data2selcropnew));
data2selcropnew(data2selcropnew == 0) = randmatrix(data2selcropnew == 0);

randmatrix = rand(size(data3selcropnew));
data3selcropnew(data3selcropnew == 0) = randmatrix(data3selcropnew == 0);

randmatrix = rand(size(data4selcropnew));
data4selcropnew(data4selcropnew == 0) = randmatrix(data4selcropnew == 0);

%% Save all selected and processed regions
close all
uisave(who, ['Regions_beads_' FileName(1:end-4) '.mat'])
