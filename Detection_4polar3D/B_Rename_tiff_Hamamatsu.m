% =========================================================================
% Script Name: Rename_tiff_Hamamatsu.m
%
% Description:
%   This script processes TIFF files acquired with Hamamatsu camera, 
%   where image filenames include an index inside parentheses 
%   (e.g., filename(0).tif). It extracts the index, increments it by 1 
%   (to shift from 0-based to 1-based numbering), and renames the files.
%
%   If the new index is 1, the filename is simplified to omit the suffix.
%   Otherwise, the suffix "_X<num>" is appended. All renamed files are 
%   copied into a subfolder named "renamed".
%
% Usage:
%   Place this script in your MATLAB workspace and update the 
%   `source_folder` path. The renamed files will be stored in 
%   "source_folder/renamed".
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel  
%   Miguel Sison              - Institut Fresnel  
%   Cesar Valades-Cruz        - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

clear all

% --- Define source and destination folders ---
source_folder = 'C:\Users\MOSAIC STORM 5\Documents\Temporary Directory for Processing\Charitra\20231009_lipidbeads_3um\multi_cube';
destination_folder = fullfile(source_folder, 'renamed');

% --- Get list of TIFF files with parentheses pattern ---
i_list = dir(fullfile(source_folder, '*(*).tif'));

% --- Create output directory if it doesn't exist ---
if ~isfolder(destination_folder)
    mkdir(destination_folder)
end

% --- Loop through each file for renaming and copying ---
for k = 1:numel(i_list)

    % Find parentheses positions
    ind1 = find(i_list(k).name == '(');
    ind2 = find(i_list(k).name == ')');

    disp(['Processing file ' num2str(k) ' of ' num2str(numel(i_list))])

    % Ensure proper index pattern
    if ind2 > ind1
        % Extract number between parentheses and increment
        num = str2double(i_list(k).name(ind1+1 : ind2-1));
        num = num + 1;  % Convert from 0-based to 1-based indexing

        % Define new filename
        if num == 1
            new_name = [i_list(k).name(1 : ind1-1) '.tif'];
        else
            new_name = [i_list(k).name(1 : ind1-1) '_X' num2str(num) '.tif'];
        end

        % Full path for old and new files
        old_file = fullfile(i_list(k).folder, i_list(k).name);
        new_file = fullfile(destination_folder, new_name);

        % Copy file with new name
        copyfile(old_file, new_file)
    end
end
