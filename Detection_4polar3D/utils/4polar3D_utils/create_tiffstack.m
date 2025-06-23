% =========================================================================
% Script Name: create_tiffstack.m
%
% Description:
%   This script loads individual TIFF image files (non-MPTIFF format) and 
%   saves them as stacked TIFFs for easier processing. Files are grouped 
%   in blocks of 'n' images and saved with incremental suffixes.
%
% Usage:
%   Use this script when TIFF images are saved as separate files and 
%   need to be grouped into stacks for further processing.
%
% Authors:
%   Miguel Sison - Institut Fresnel  
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

clear all
clc

% --- Define input and output paths ---
img_fldr = 'E:\Miguel\Raw Data\20221103_B16\n2';  % Folder containing TIFF images
img_fname = 'Image1_';                            % Common prefix of the TIFF files
dest_fldr = 'D:\Miguel\Processing Data\20221103_B16'; % Destination folder for saving stacked TIFFs
nimg_fname = 'B16_AF568_UQ_180ms_NoBin_n2';       % Base name for output files

% --- Get list of TIFF files ---
img_list = dir(fullfile(img_fldr, [img_fname '*.tif']));
N = numel(img_list);  % Total number of TIFF files

% --- Set number of images per stack ---
n = 100;

% --- Process full stacks of 'n' images ---
for k0 = 0:(floor(N/n)-1)
    
    for k1 = 1:n
        % Get full filename
        fname = fullfile(img_list(k1 + k0 * n).folder, img_list(k1 + k0 * n).name);
        disp(fname)
        
        % Load TIFF image and add to stack
        temp_img = loadtiff(fname);
        temp_stk(:,:,k1) = temp_img;
    end

    % Save stack to new TIFF file
    if k0 == 0
        % First stack saved with base name
        savetiff(fullfile(dest_fldr, [nimg_fname '.tif']), temp_stk);
    else
        % Subsequent stacks get numbered suffix
        savetiff(fullfile(dest_fldr, [nimg_fname '_X' num2str(k0+1) '.tif']), temp_stk);
    end

    clear temp_img temp_stk
    last_file_num = k1 + k0 * n;  % Track index of last processed file
end

% --- Process remaining images (less than 'n') if any ---
inds = (last_file_num + 1):N;

if numel(inds) ~= 0
    for k1 = 1:numel(inds)
        fname = fullfile(img_list(inds(k1)).folder, img_list(inds(k1)).name);
        disp(fname)
        
        temp_img = loadtiff(fname);
        temp_stk(:,:,k1) = temp_img;
    end

    % Save remaining images to new TIFF file
    savetiff(fullfile(dest_fldr, [nimg_fname '_X' num2str(k0+2) '.tif']), temp_stk);
    clear temp_img temp_stk
end

% --- Done ---
disp('...done')
