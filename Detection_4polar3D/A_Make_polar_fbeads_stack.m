% =========================================================================
% Script Name: Makepolar_fbeads_stack.m
%
% Description:
%   This script creates a stack of polar and isobead calibration images.
%   It reads all TIFF files in the specified folder and builds an extended 
%   image stack to comply with the processing requirements of the 4polar3D 
%   MATLAB pipeline, which expects stacks in multiples of 100 frames.
%
%   For each image, it replicates the frame N times (default N = 20) to 
%   generate a sufficiently large stack. It also inserts replicated 
%   "background" frames at the beginning and end of the stack to pad the 
%   stack correctly.
%
%   Use this when preparing calibration datasets for polar and isobeads.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel  
%   Miguel Sison              - Institut Fresnel  
%   Cesar Valades-Cruz        - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

% --- Define input folder containing calibration TIFF files ---
fldr = ['C:\Users\MOSAIC STORM 5\Documents\Temporary Directory for Processing\Charitra\20231010_calibration\iso_multicube'];
img_list = dir(fullfile(fldr, '*.tif'));

% --- Read first image to get image dimensions ---
temp = imfinfo(fullfile(fldr, img_list(1).name));

% --- Number of times to replicate each frame ---
N = 20;

% --- Pre-allocate image stack with appropriate size ---
% The total stack size is rounded to the next multiple of 100
img_stck = NaN(temp.Height, temp.Width, ceil(numel(img_list) * 2 / 10) * 100);

% --- Load each image and replicate N times in the stack ---
for k = 1:numel(img_list)
    temp = loadtiff(fullfile(fldr, img_list(k).name));
    img_stck(:, :, (k - 0.5) * N + (1:N)) = repmat(temp, [1, 1, N]);
    disp(k)  % Display progress
end

% --- Pad the stack with background frames at the start and end ---
% The padding uses the mean image calculated while ignoring NaN values
img_stck(:, :, [1:N/2 end - N/2 + 1:end]) = repmat(mean(img_stck, 3, "omitnan"), [1, 1, N]);

% --- Replace NaN values with zero ---
img_stck(isnan(img_stck)) = 0;

% --- Save the final stack ---
savetiff([fldr '_img_stack.tif'], img_stck);
