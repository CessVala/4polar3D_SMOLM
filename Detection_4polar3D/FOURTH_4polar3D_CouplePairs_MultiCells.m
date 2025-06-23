% =========================================================================
% Script Name: FOURTH_4polar3D_CouplePairs_MultiCells.m
%
% Description:
%   This script performs the coupling of the two previously analyzed polarized projection pairs:
%   the 0°/90° pair and the 45°/135° pair, for multi-cell datasets in the 4polar3D pipeline.
%
%   The script:
%     - Asks the user to select the parent folder that contains all the folders generated 
%       by the previous three 4polar3D scripts.
%     - For each subfolder, automatically loads the results from the SECOND and THIRD 
%       4polar3D analysis steps.
%     - Couples the trajectories from the two polarized pairs (P0 and P45).
%     - Computes and saves the merged trajectory table (tab_traj_4x4) that contains pairing information.
%
% Instructions:
%   1. Run this script.
%   2. Select the parent folder that contains all the individual folders processed 
%      by the previous three 4polar3D steps.
%   3. The script will automatically process all folders inside.
%   4. Each coupled result will be saved as *_4x4_New.mat inside the corresponding folder.
%
% Notes:
%   - This script must be run **after**:
%       - FIRST_4polar3D_Preprocess.m
%       - SECOND_4polar3D_AnalyzePair_0_90_MultiCells.m
%       - THIRD_4polar3D_AnalyzePair_45_135_MultiCells.m
%   - The SECOND and THIRD scripts can be run **in parallel in two MATLAB windows.**
%   - This coupling step matches the vectors between the corresponding polarized projections.
%   - The STHERM_param script must be properly configured before running this script.
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

%% Load STHERM parameters
STHERM_param

%% Select the parent folder containing all processed cells
PathName_dum = uigetdir('D:\4POLARSTORM\_RAW DATA\');

%% List all files and folders in the selected directory
files = dir(PathName_dum);
dirFlags = [files.isdir]; % Identify folders
subFolders = files(dirFlags);
Num_folders = length(subFolders) - 2; % Exclude '.' and '..'

%% Process each subfolder
for subf = 1:Num_folders

    %% Get the subfolder path
    PathName = [PathName_dum '\' subFolders(subf+2).name '\'];
    
    %% Load 0°/90° result
    filename_dum = dir([PathName '*incl0Newa*']);
    filename1 = [[PathName filename_dum.name]];
    
    if isempty(filename_dum)
        continue % Skip if no matching file found
    end

    load(filename1, 'tab_traj')

    %% Prepare P0 pair trajectories
    tab_trajP0init = tab_traj;
    tab_trajP0init = [tab_trajP0init; 1:size(tab_traj, 2)];
    tab_trajP0 = tab_traj(:, 2:2:end);
    tab_trajP0 = [tab_trajP0; 2:2:size(tab_traj, 2)];

    %% Load 45°/135° result
    filename_dum = dir([PathName '*incl45New*']);
    FileName = filename_dum.name;
    filename2 = [[PathName filename_dum.name]];
    load(filename2, 'tab_traj')

    %% Prepare P45 pair trajectories
    tab_trajP45init = tab_traj;
    tab_trajP45init = [tab_trajP45init; 1:size(tab_traj, 2)];
    tab_trajP45 = tab_traj(:, 2:2:end);
    tab_trajP45 = [tab_trajP45; 2:2:size(tab_traj, 2)];

    %% Reposition trajectories for proper coupling
    displacementjmin1 = min(tab_trajP45(3, :));
    tab_trajP45(3, :) = tab_trajP45(3, :) - displacementjmin1 + 1;

    displacementjmin = min(tab_trajP0(3, :));
    tab_trajP0(3, :) = tab_trajP0(3, :) - displacementjmin + 1;
    displacementj = max(tab_trajP0(3, :));

    tab_trajP45(3, :) = tab_trajP45(3, :) + displacementj + 100;
    tab_trajsel = [tab_trajP0 tab_trajP45];

    %% Plot coupling result
    figure(1)
    set(gcf, 'Color', 'w')
    scatter(tab_trajP45(3, :), tab_trajP45(2, :), 10, 'Filled')
    hold on
    scatter(tab_trajP0(3, :), tab_trajP0(2, :), 10, 'Filled')
    set(gca, 'YDir', 'reverse')
    axis('image')

    sep = round((max(tab_trajsel(3, :))) / 2);

    %% Run coupling
    try
        tic;
        [tab_trajcoupled, tab_traj_allcoupled] = STHERMver3([NumStart, NumImgVector], [NumStart, NumImg], ...
            PathName, 'Image_', PathName, [FileName(1:end-11) '_New0-45.mat'], tab_trajsel);
        toc;
    catch exception
        rethrow(exception);
    end

    %% Prepare final trajectory structure
    newtrab_traj = zeros(size(tab_trajcoupled, 1), 2 * size(tab_trajcoupled, 2));
    for i = 1:2:size(tab_trajcoupled, 2)

        if tab_trajcoupled(9, i+1) == 0
            add1 = tab_trajcoupled(:, i+1);
        else
            add1 = tab_trajP0init(:, tab_trajcoupled(9, i+1) - 1);
        end

        if tab_trajcoupled(9, i) == 0
            add2 = tab_trajcoupled(:, i);
        else
            add2 = tab_trajP45init(:, tab_trajcoupled(9, i) - 1);
        end

        newtrab_traj(:, (2*i-1):(2*i-1+3)) = [tab_trajcoupled(:, i+1) tab_trajcoupled(:, i) add1 add2];
    end

    %% Final merged trajectory table
    tab_traj_4x4 = [newtrab_traj; sort(repmat(1:size(tab_trajcoupled, 2) / 2, [1, 4]))];

    %% Save the result
    name_savefile = [PathName FileName(1:end-11) '_4x4_New.mat'];
    save(name_savefile);

    %% Display processed folder info
    PathName
    FileName

end
