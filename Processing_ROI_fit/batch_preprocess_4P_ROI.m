%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: batch_preprocess_4P_ROI.m
%
% Purpose:
%   Batch processing of 4polar results using a precomputed K matrix and fitting approach.
%   Steps:
%     1. Fit raw data using 4polar model
%     2. Generate polarization-resolved images
%     3. Export per-ROI data for further analysis
%
% Authors: Charitra S. Senthil Kumar, Cesar Valades Cruz
% Date: March 17, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 1: FIT RAW DATA WITH 4POLAR MODEL
clear all

% Load matrix (K)
load('C:\Users\MOSAIC STORM 5\Documents\Temporary Directory for Processing\demo procesing\K_matrix.mat');
K = K_final;

% Define input directory containing 4x4 files
hdir = 'C:\Users\MOSAIC STORM 5\Documents\Temporary Directory for Processing\demo procesing\New folder';

% List all *_4x4_New.mat files
flist = dir([hdir '\*_4x4_New.mat']);
sdir = [hdir];
% Define output directory for fitted results
% sdir = [hdir '\fit_results'];
% if ~isfolder(sdir)
%      mkdir(sdir)
% end

str_4 = '_4x4_New'; % Identifier string for raw input files

% Loop over all files for fitting
for k = 1:numel(flist)
    fname = fullfile(flist(k).folder, flist(k).name); % Full path to input
    sname = fullfile(sdir, strrep(flist(k).name, str_4, '_fit_results')); % Replace suffix
    [fit_results, old_data] = fit_4P_ISOmodel(K, fname); % Fit model
    save(sname, 'fit_results', 'old_data', '-v7.3'); % Save results
end

%% STEP 2: GENERATE POLARIZATION-RESOLVED IMAGES
str_f = 'fit_results';

% Reload file lists
flist = dir([hdir '\*_4x4_New.mat']);
list_fit_res = dir([sdir '\*' str_f '.mat']);


% Set constants
g_sig   = 0.020;
px_size = 0.130;          % in microns
img_px  = 2^11;           % Output image size in pixels


% Process each file again
for k = 1:numel(flist)
    disp(k); % Display progress
    clear data fit_results img_struct

    % Input file paths
    fname    = fullfile(flist(k).folder, flist(k).name);
    fit_file = fullfile(list_fit_res(k).folder, list_fit_res(k).name);

    % Recompute raw 4polar quantities
    data = calc_4P_Polar_params_ISO(fname, K);
    load(fit_file)

    % Ensure number of fitted values matches detections
    if numel(data.It) ~= size(fit_results.rho,1)
        disp('Number of detections are not equal. Check file names!')
        break
    else
        disp('Files compatible. Proceeding...')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Drift correction block (optional, currently disabled) from
    % ThunderSTORM drift correction.
%     drift_temp = readmatrix( [flist(k).folder '\' flist(k).name(1:2) '_allDriftTable.csv']);    
%     y_loc = (data.Y(:, 1) + drift_temp(data.T,3)) * px_size;
%     x_loc = (data.X(:, 1) + drift_temp(data.T,2)) * px_size;
%     x_loc = x_loc - min(x_loc(:)) + 1;
%     y_loc = y_loc - min(y_loc(:)) + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Use raw localization positions (no drift correction)
    y_loc = data.Y(:, 1) * px_size;
    x_loc = data.X(:, 1) * px_size;

    % Convert total intensity to photons
    It = data.It * 0.24 / 0.775;
    It_max = max(It);

    % Extract polarization parameters
    rho   = fit_results.rho(:,2);
    eta   = fit_results.eta(:,2);
    delta = fit_results.delta(:,2);

    
    maxY = ceil(max(y_loc)/5) * 5; % get max position Y-direction rounded up to nearest 10um
    maxX = ceil(max(x_loc)/5) * 5;
    lim_img = max(maxX, maxY);
    x_edges = linspace(0, lim_img, img_px+1);
    y_edges = linspace(0, lim_img, img_px+1);

    % Calculate bin centers
    bn_size = mean(diff([x_edges; y_edges], 1, 2), 2);
    x_bin = x_edges(1:end-1) + bn_size(1)/2;
    y_bin = y_edges(1:end-1) + bn_size(2)/2;

    % Create polarization-resolved images
    img_struct = create_4P_images(x_loc, y_loc, It, rho, eta, delta, x_edges, y_edges);

    % Store bin metadata
    img_struct.x_bin   = x_bin;
    img_struct.y_bin   = y_bin;
    img_struct.bn_size = bn_size;

    % Save image structure
    sname = fullfile(list_fit_res(k).folder, ...
		strrep(list_fit_res(k).name, str_f, 'images'));
    save(sname, 'img_struct', '-v7.3')
end

%% STEP 3: EXPORT ROI-LEVEL DATA STRUCTURES
clear all

% Reload K matrix
load('C:\Users\MOSAIC STORM 5\Documents\Temporary Directory for Processing\demo procesing\K_matrix.mat');
K = K_final;

% Define paths
data_dir = 'C:\Users\MOSAIC STORM 5\Documents\Temporary Directory for Processing\demo procesing\New folder';
str_4 = '4x4_New';
str_f = 'fit_results';
str_i = 'images';

% List input files again
files_list = dir(fullfile(data_dir, ['*' str_4 '*.mat']));

% Constants
px_size = 0.130; 
count2ph = 0.24 / 0.775;

% Define output directory
save_dir = [data_dir '\ROI_Processing_directory\'];
if ~isfolder(save_dir)
    mkdir(save_dir)
end

% Loop over each dataset
for k = 1:numel(files_list)
    fname = fullfile(files_list(k).folder, files_list(k).name)
    data = calc_4P_Polar_params_ISO(fname, K);

    % Load fit results from fit_4P_ISOmodel
    fname = fullfile(files_list(k).folder, strrep(files_list(k).name, str_4, str_f));
    load(fname);

    % Load composite 4P orientation images
    fname = fullfile(files_list(k).folder, strrep(files_list(k).name, str_4, str_i));
    load(fname);

    % Create full-frame mask
    img_struct.img_mask = ones(size(img_struct.img_strm(:,:,1)));

    % Construct output structure
    data_struct.rho           = fit_results.rho(:,2);
    data_struct.eta           = fit_results.eta(:,2);
    data_struct.delta         = fit_results.delta(:,2);
    data_struct.rmse          = fit_results.rmse(:,2);
    data_struct.radius_mean     = old_data.radius_mean(:,1);
    data_struct.radius_max      = old_data.radius_max(:,1);
    data_struct.radius          = old_data.radius(:,:);
    data_struct.y_loc         = img_struct.y_loc;  % x and y coordinates in real world measurements (micrometer)
    data_struct.x_loc         = img_struct.x_loc;
    data_struct.frame         = data.T;
    data_struct.I_tot         = data.It(:,1) * count2ph;
    data_struct.I             = data.I(:,:);
    data_struct.loc_precision = data.sigma(:,:);
    data_struct.noise         = data.noise(:,:);

    % Define output filename
    save_name = strrep(files_list(k).name, str_4, 'Processing_data');
    save_name = fullfile(save_dir, save_name);

    % Save result
    save(save_name, 'img_struct', 'data_struct', 'K', 'px_size', 'count2ph', '-v7.3')
end
