% =========================================================================
% Script Name: sample_roi_4p_v2.m
%
% Description:
%   This script is used to:
%     - Load a 4polar3D dataset or a previously saved ROI processing file.
%     - Initialize or update ROI selection using an interactive interface.
%     - Allow the user to draw and save freehand ROIs.
%     - Display all selected ROIs.
%     - Generate 4polar3D stick plots and histograms for selected ROIs.
%     - Save the ROI data as Excel files.
%
% Instructions:
%   - Run the script and follow the prompts to load data, select save location, and draw ROIs.
%   - The script supports working with both new and previously saved ROI sessions.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

%%
clear all
% clc
fs = 14;

% Load the .mat data file
[fname, fpath] = uigetfile('*.mat');
load(fullfile(fpath, fname));

% Create a new save file name based on current time
time_temp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));

% Detect if the loaded file is an ROI save file
try datetime(fname(1:19), 'InputFormat', 'yyyy-MM-dd_HH-mm-ss');
    opts.Interpreter = 'tex';
    opts.Default = 'Yes';
    answer = questdlg('\bf\fontsize{14}You have loaded an "ROI processing SAVE file". Do you want to CREATE NEW save file?',...
        '', 'Yes', 'Still Yes', 'No', opts);
    switch answer
        case 'Yes'
            save_file_name = fullfile(fpath, [time_temp '_' fname(21:end)]);
        case 'No'
            save_file_name = fullfile(fpath, fname);
        case 'Still Yes'
            save_file_name = fullfile(fpath, [time_temp '_' fname(21:end)]);
    end
    clear opts
catch
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'non-modal';
    uiwait(msgbox('\bf\fontsize{12}You have loaded an "Original data file". A separate SAVE file will be created.', '', CreateStruct))
    save_file_name = fullfile(fpath, [time_temp '_' fname]);
    clear CreateStruct
end

% Create figure for ROI selection
try any([isa(roi_struct, "struct"), isfield(roi_struct, 'list_form'), isfield(roi_struct, 'mask_form')]);
    [~, ax_i, ax_r, ax_e, ax_d, ax_m, t_layout] = gen_fig_4P_roi_select(img_struct, roi_struct);
catch
    [roi_struct, ax_i, ax_r, ax_e, ax_d, ax_m, t_layout] = gen_fig_4P_roi_select(img_struct);
    roi_data = roi_get_4Pdata(0, data_struct, roi_struct.list_form);
end

[old_sdir, current_save_name] = fileparts(save_file_name);

CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';
uiwait(msgbox('\bf\fontsize{12}Select SAVE folder.', '', CreateStruct))

% Ask user to select save directory
fpath_new = uigetdir;
save_file_name = fullfile(fpath_new, [current_save_name '.mat']);

% Save current data
save(save_file_name, "img_struct", "roi_struct", "data_struct", "roi_data", "-v7.3");

%% Ask user to select ROIs
roi_struct = roi_select_4P_freehand(roi_struct, ax_i, ax_r, ax_e, ax_d, ax_m);
temp = size(roi_struct.list_form, 2);

roi_data(temp) = roi_get_4Pdata(temp - 1, data_struct, roi_struct.list_form);

% Save updated ROI data
save(save_file_name, "roi_struct", "roi_data", "-append");

disp(['You currently have ' num2str(size(roi_struct.list_form, 2) - 1) ' user defined ROI...'])

%% Display all ROIs

figure;
tiledlayout('flow')

% Display STORM image
nexttile
imagesc(img_struct.img_strm(:, :, 2));
axis image; set(gca, 'YDir', 'normal')
set(gca, 'Colormap', gray)
clim(prctile(img_struct.img_strm(:, :, 2), [0 99.975], 'all'))
title('STORM image')

% Display all ROI masks
nexttile
temp_im = imagesc(sum(roi_struct.mask_form(:, :, 2:end) .* reshape(1:(size(roi_struct.mask_form, 3) - 1), 1, 1, []), 3));
axis image
set(gca, 'YDir', 'normal')
for k = 2:size(roi_struct.mask_form, 3)
    [row, col] = find(roi_struct.mask_form(:, :, k));
    text(max(col) + 2, mean(unique(row)), num2str(k - 1), 'FontSize', 14, 'Color', [1 1 0] * 0.9)
end
set(gca, 'Colormap', jet)

% Display individual ROI overlays
for k = 2:size(roi_struct.mask_form, 3)
    nexttile
    temp_im = imagesc(min(img_struct.img_strm(:, :, 2), prctile(img_struct.img_strm(:, :, 2), 99.99, 'all'))); % STORM image as background
    %     temp_im = imagesc(min(imgFluo, prctile(imgFluo(:),99.99))); % use fluoresence image as background
    temp_im.AlphaData = (~roi_struct.mask_form(:, :, k)) + roi_struct.mask_form(:, :, k) * 0.55;
    set(gca, 'Color', [1 1 0] .* 0.9);
    axis image; set(gca, 'YDir', 'normal');
    title(['ROI ' num2str(k - 1)]);
    set(gca, 'Colormap', gray)
end

%% Save ROI data to Excel file
save_roi_data(roi_data)

%% Apply filtering and generate stick plots
rho = data_struct.rho;
eta = data_struct.eta;
delta = data_struct.delta;
tot_int = data_struct.I_tot;
filter = eta < deg2rad(80) & delta > deg2rad(10) & tot_int > 1000;

imgStrm = img_struct.img_strm(:, :, 2);
x_bin = img_struct.x_bin;
y_bin = img_struct.y_bin;
x_loc = data_struct.x_loc;
y_loc = data_struct.y_loc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edited on 14/03/2023 to show one or more than one ROI simultaneously
% figure;
% imagesc(imgStrm)
% %,'XData',[0 50],'YData',[0 50]);
% colormap gray
% clim([0 1000])
% hold on
% scatter(y_cent(1:end),x_cent(1:end),100,'o','filled','white')
% set(gca,'YDir','normal')

roi_num = [1:numel(roi_data)];
roi_list = zeros(numel(x_loc), 1);
sl = 0.1;      % Stick length (microns)
sw = sl * 0.1; % Stick width (microns)
sa = 0.95;     % Stick opacity 0-1

flag = numel(roi_num);

if flag == 1

    % Single ROI
    roi_list = roi_struct.list_form(:, roi_num);
    roi_list = roi_list & filter;

    % Make sticks figures

    %%%% Sticks image with Rho as colormap

    % time
    c_line = data_struct.frame; 	% assign colors to be based on frame number

    fig_t = figure('Position', [50 50 940 700], 'Resize', 'off');
    [ax_img_r, ax_stk_r, fig_t] = create_single_4Psticks_image(fig_t, imgStrm, x_bin, y_bin, ...
        x_loc(roi_list), y_loc(roi_list), rho(roi_list), c_line(roi_list), sl, sw, sa);

    colormap(ax_stk_r, "hsv");
    ax_stk_r.CLim = [0 2000];				% max number of frames
    title('frames')

    % Stick plot: rho
    c_line = rad2deg(rho);		% assign colors to be based on Rho

    fig_r = figure('Position', [50 50 940 700], 'Resize', 'off');
    [ax_img_r, ax_stk_r, fig_r] = create_single_4Psticks_image(fig_r, imgStrm, x_bin, y_bin, ...
        x_loc(roi_list), y_loc(roi_list), rho(roi_list), c_line(roi_list), sl, sw, sa);

    colormap(ax_stk_r, "hsv");
    ax_stk_r.CLim = [0 180];		 % set range of colors (here it is 0-180) because Rho goes from 0 to 180 degrees
    title('\rho')

    % Stick plot: eta
    c_line = rad2deg(eta);		% assign colors to be based on eta

    fig_e = figure('Position', [50 50 940 700], 'Resize', 'off');
    [ax_img_r, ax_stk_r, fig_e] = create_single_4Psticks_image(fig_e, imgStrm, x_bin, y_bin, ...
        x_loc(roi_list), y_loc(roi_list), rho(roi_list), c_line(roi_list), sl, sw, sa);

    colormap(ax_stk_r, "plasma");
    ax_stk_r.CLim = [0 90];		% set range of colors (here it is 0-90) because eta goes from 0 to 90 degrees
    title('\eta')

    % %binned eta
    %
    % [ax_img_e0, ax_stk_e0, fig_e0] = bin_eta(eta, roi_list, imgStrm,  x_bin, y_bin, x_loc, y_loc, rho, c_line, sl, sw, sa);
    % Stick plot: delta
    c_line = rad2deg(delta); 	% assign colors to be based on delta

    fig_d = figure('Position', [50 50 940 700], 'Resize', 'off');
    [ax_img_r, ax_stk_r, fig_d] = create_single_4Psticks_image(fig_d, imgStrm, x_bin, y_bin, ...
        x_loc(roi_list), y_loc(roi_list), rho(roi_list), c_line(roi_list), sl, sw, sa);

    colormap(ax_stk_r, "jet");
    ax_stk_r.CLim = [0 180];	% set range of colors (here it is 0-180) because delta goes from 0 to 180 degrees
    title('\delta')

    % Histograms for ROI
    figure;
    t_layout = tiledlayout(2, 2);
    nexttile;
    histogram(rad2deg(rho(roi_list)), 'BinLimits', [0, 180]); title('\rho')
    nexttile;
    histogram(rad2deg(eta(roi_list)), 'BinLimits', [0, 90]); title('\eta')
    nexttile;
    histogram(rad2deg(delta(roi_list)), 'BinLimits', [0, 180]); title('\delta')
    nexttile;
    histogram(tot_int, 'BinLimits', [0, 15000]); title('Total intensity');

else

    % Multiple ROIs
    for ij = 1:numel(roi_num) - 1
        temp_roi_list = roi_struct.list_form(:, roi_num(ij) + 1);
        roi_list = or(roi_list, temp_roi_list);
    end
    roi_list = roi_list & filter;

    % Stick plot: time (frame number)

    c_line = data_struct.frame;			% assign colors to be based on frame number

    fig_t = figure('Position', [50 50 940 700], 'Resize', 'off');
    [ax_img_t, ax_stk_t, fig_t] = create_single_4Psticks_image(fig_t, imgStrm, x_bin, y_bin, ...
        x_loc(roi_list), y_loc(roi_list), rho(roi_list), c_line(roi_list), sl, sw, sa);

    colormap(ax_stk_t, "turbo");
    ax_stk_t.CLim = [0 40000];			% max number of frames
    title('frames')

    % Stick plot: rho
    c_line = rad2deg(rho);			% assign colors to be based on Rho

    fig_r = figure('Position', [50 50 940 700], 'Resize', 'off');
    [ax_img_r, ax_stk_r, fig_r] = create_single_4Psticks_image(fig_r, imgStrm, x_bin, y_bin, ...
        x_loc(roi_list), y_loc(roi_list), rho(roi_list), c_line(roi_list), sl, sw, sa);

    colormap(ax_stk_r, "hsv");
    ax_stk_r.CLim = [0 180];	% set range of colors (here it is 0-180) because Rho goes from 0 to 180 degrees
    title('\rho')

    % Stick plot: eta
    c_line = rad2deg(eta);		% assign colors to be based on eta

    fig_e = figure('Position', [50 50 940 700], 'Resize', 'off');
    [ax_img_e, ax_stk_e, fig_e] = create_single_4Psticks_image(fig_e, imgStrm, x_bin, y_bin, ...
        x_loc(roi_list), y_loc(roi_list), rho(roi_list), c_line(roi_list), sl, sw, sa);

    colormap(ax_stk_e, "plasma");
    ax_stk_e.CLim = [0 80];	 % set the color range: eta goes from 0 to 90 degrees. After filtering, it is adjusted to 80 (see the 'filter' variable above)
    title('\eta')

    % %binned eta
    % [ax_img_e0, ax_stk_e0, fig_e0] = bin_eta(eta, roi_list, imgStrm,  x_bin, y_bin, x_loc, y_loc, rho, c_line, sl, sw, sa);
    % Stick plot: delta
    c_line = rad2deg(delta);	% assign colors to be based on delta

    fig_d = figure('Position', [50 50 940 700], 'Resize', 'off');
    [ax_img_d, ax_stk_d, fig_d] = create_single_4Psticks_image(fig_d, imgStrm, x_bin, y_bin, ...
        x_loc(roi_list), y_loc(roi_list), rho(roi_list), c_line(roi_list), sl, sw, sa);

    colormap(ax_stk_d, "jet");
    ax_stk_d.CLim = [0 180];	% set range of colors (here it is 0-180) because delta goes from 0 to 180 degrees
    title('\delta')

    % Link axes
    linkaxes([ax_stk_r ax_img_r ax_stk_e ax_img_e ax_stk_d ax_img_d], 'xy');

    % Histograms for ROI
    figure;
    t_layout = tiledlayout(2, 2);
    nexttile;
    histogram(rad2deg(rho(roi_list)), 'BinLimits', [0, 180]); title('\rho')
    nexttile;
    histogram(rad2deg(eta(roi_list)), 'BinLimits', [0, 90]); title('\eta')
    nexttile;
    histogram(rad2deg(delta(roi_list)), 'BinLimits', [0, 180]); title('\delta')
    nexttile;
    histogram(tot_int, 'BinLimits', [0, 15000]); title('Total intensity');
end