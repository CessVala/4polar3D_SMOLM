% =========================================================================
% Script Name: SECOND_4polar3D_AnalyzePair_0_90_MultiCells.m
%
% Description:
%   This script performs the analysis of the first polarized projection pair: 0° and 90°.
%   It is intended for processing multiple cells (multi-selected stacks) in the 4polar3D pipeline.
%
%   The script:
%     - Allows the user to multi-select the first stack(s) of the cells to be analyzed.
%     - Requires the user to draw a single ROI around the top right (0°) polarized projection.
%     - The user must drag the same ROI across to the top left (90°) polarized projection to ensure 
%       the selected regions have exactly the same number of pixels.
%     - Processes each selected stack sequentially (multi-cell analysis).
%
% Instructions:
%   1. Run this script.
%   2. Multi-select the first TIFF stack(s) of the cells to be analyzed.
%   3. Draw a single ROI around the top right (0°) polarized projection.
%   4. Drag this same ROI across to the top left (90°) polarized projection to guarantee 
%      identical pixel coverage.
%   5. The script will automatically apply the alignment and 4polar3D analysis to each selected stack.
%
% Notes:
%   - This script must be run after FIRST_4polar3D_Preprocess.m.
%   - The STHERM_param script must be properly configured before starting.
%   - This script is designed for multi-cell analysis.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel  
%   Miguel Sison              - Institut Fresnel  
%   Cesar Valades-Cruz        - Institute of Hydrobiology (IHB), CAS
%   Arturo Vesga
%
% Date: June 2025
% =========================================================================

%% Initialize
clear all
close all
clc

%% Define global variables for ROI and alignment
global name_imagecor
global contourIm
global verticesleft
global rectveL 
global verticesright
global rectveR

%% Load STHERM parameters
STHERM_param

%% Select the image stacks to be analyzed
[FileName_dum, PathName, FilterIndex] = uigetfile('D:\4POLARSTORM\_RAW DATA\*.tif', 'MultiSelect', 'on');

%% Prepare file list
if iscell(FileName_dum)
    num_selected = size(FileName_dum, 2); % Number of selected files
else
    num_selected = 1;
    dummy = FileName_dum; % Ensure consistent format
    clear FileName_dum
    FileName_dum{1} = dummy;
end

%% Process each selected stack
for k_index = 1:num_selected
    
    FileName = FileName_dum{k_index}; % Get the file name

    tiff_file = [PathName FileName];
    tiffInfo = imfinfo(tiff_file);
    sizestack = numel(tiffInfo); % Number of frames in the stack

    % Path to the corrected images
    name_imagecor = [PathName '\' FileName(1:end-4) '\images1Corrected\Image_'];

    %% Calculate frame to align
    if NumStart > 100
        Imagetoalign = (floor(NumStart / 100) - 1) * 100 + 1;
    else
        Imagetoalign = 1;
    end

    contourIm = 10; % Optional contour padding
    filenamecor = strcat(name_imagecor, num2str(Imagetoalign));
    m = matfile(filenamecor);
    data1 = m.data;

    %% Align and select ROIs
    if k_index == 1
        % First stack: manually select ROIs
        [data1, verticesleft, rectveL, verticesright, rectveR] = alignIm(data1, contourIm);
        assignin('base', ['data' num2str(Imagetoalign)], data1);
    else
        % Subsequent stacks: use the same ROIs
        [data1, verticesleft, rectveL, verticesright, rectveR] = alignIm(data1, contourIm, verticesleft, rectveL, verticesright, rectveR);
        assignin('base', ['data' num2str(Imagetoalign)], data1);
    end

    sep = round(size(data1(:,:,1), 2) / 2); % Image separation for ROI division (optional)

    %% Run analysis for the 0° and 90° projections
    try
        tic;
        [tab_traj, tab_traj_all] = STHERMver2([NumStart, NumImgVector], [NumStart, NumImg], ...
            [PathName '\' FileName(1:end-4) '\images1Corrected\'], 'Image_', ...
            [PathName '\' FileName(1:end-4) '\'], ...
            [FileName(1:end-4) '_' num2str(pfa) '_r_' num2str(r0) '_W_' num2str(wn) '_incl0Newa.mat']);
        toc;
    catch exception
        rethrow(exception); % Display error details
    end

    %% Display processed file info
    PathName
    FileName
    
end
