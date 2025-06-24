% =========================================================================
% Function Name: gen_fig_4P_roi_select.m
%
% Description:
%   This function generates an interactive figure for ROI selection in the
%   4polar3D visualization pipeline. It displays the STORM image, the ROI mask,
%   and the 4polar3D parameter maps (rho, eta, delta) in a tiled layout.
%
%   The function:
%     - Prepares the image layers and parameter maps.
%     - Handles optional input arguments for ROI structures and display settings.
%     - Creates a multi-panel figure with linked axes for ROI selection.
%
% Instructions:
%   - Call this function to open the ROI selection interface.
%   - Optionally provide existing ROI data and display thresholds.
%
% Notes:
%   - Gaussian bluring is applied to the STORM image for display.
%   - Pixel binning sizes are considered for filtering and scaling.
%
% Inputs:
%   - img_struct: Structure containing the STORM image, parameter maps, and binning info.
%   - varargin (optional):
%       1st argument: Existing ROI structure.
%       2nd argument: Gaussian sigma for bluring.
%       3rd argument: Intensity percentile cutoff for display.
%
% Outputs:
%   - roi_struct: Updated or initialized ROI structure.
%   - ax_i: Axis handle for the STORM image.
%   - ax_r: Axis handle for the rho map.
%   - ax_e: Axis handle for the eta map.
%   - ax_d: Axis handle for the delta map.
%   - ax_m: Axis handle for the ROI mask.
%   - t_layout: Tiled layout object.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function [roi_struct, ax_i, ax_r, ax_e, ax_d, ax_m, t_layout] = gen_fig_4P_roi_select(img_struct, varargin)

% Input handling for optional ROI and display parameters
if nargin == 1
    g_sig = 0.020; % Gaussian sigma for bluring (um)
    a_lim = 99.975; % Percentile cutoff for display
elseif nargin == 2
    roi_struct = varargin{1};
    g_sig = 0.020; % Default bluring
    a_lim = 99.975;
elseif nargin == 3
    g_sig = varargin{2};
    a_lim = 99.975;
elseif nargin == 4
    g_sig = varargin{2};
    a_lim = varargin{3};
else
    errordlg('Wrong Inputs!!!')
    return
end

% Set figure aspect ratio
ar = 0.45;

% Extract image data and binning
bn_size = img_struct.bn_size;
x_bin = img_struct.x_bin;
y_bin = img_struct.y_bin;

img_strm = img_struct.img_strm(:,:,1);
img_rho = img_struct.img_rho(:,:,2);
img_eta = img_struct.img_eta(:,:,2);
img_delta = img_struct.img_delta(:,:,2);

% Initialize or load ROI structure
roi_struct.x_loc = img_struct.x_loc;
roi_struct.y_loc = img_struct.y_loc;
roi_struct.Xbinned = img_struct.Xbinned;
roi_struct.Ybinned = img_struct.Ybinned;

if ~(isfield(roi_struct, 'list_form') && isfield(roi_struct, 'mask_form'))
    roi_struct.list_form = true(numel(roi_struct.x_loc), 1);
    roi_struct.mask_form = ones(size(img_strm));
end

% Generate combined mask display
img_mask = sum(roi_struct.mask_form(:,:,1:end).*(cumsum(ones(1, 1, size(roi_struct.list_form, 2)), 3)-1), 3);

% Apply Gaussian bluring to the STORM image
img_strmg = imgaussfilt(img_strm, g_sig ./ bn_size(1));
a_data = min(img_strmg, prctile(img_strmg(:), a_lim));
a_data = a_data ./ max(a_data);

% Set figure size based on screen resolution
scr_size = get(0, 'MonitorPositions');
[~, a2] = max(scr_size(:, 3) .* scr_size(:, 4));
temp = [scr_size(a2, 3) scr_size(a2, 3) * ar;
    scr_size(a2, 4) / ar scr_size(a2, 4)];


figure('Position', [scr_size(a2,1:2) + (scr_size(a2,3:4)-temp(prod(temp,2) < prod(scr_size(a2,3:4)),:)*0.95)/2 ...
    temp(prod(temp,2) < prod(scr_size(a2,3:4)),:)*0.95] )

% Create tiled layout
t_layout = tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

% Plot STORM image
ax_i = nexttile(t_layout, 1, [2, 2]);
imagesc(x_bin, y_bin, img_strmg);
axis image
set(ax_i, 'YDir', 'normal')
set(ax_i, 'Color', 'Yellow')
colormap(ax_i, 'gray')
set(ax_i, 'CLim', prctile(img_strmg(:), [0 a_lim], 'all'))
box off
title('STORM image')

% Plot ROI mask
ax_m = nexttile(t_layout, 3);
imagesc(x_bin, y_bin, img_mask);
axis image
set(ax_m, 'YDir', 'normal')
colormap(gca, 'jet')
box off
title('ROIs')

% Plot rho map
ax_r = nexttile(t_layout, 4);
imagesc(x_bin, y_bin, rad2deg(img_rho), 'AlphaData', ~isnan(img_rho) .* a_data);
axis image
set(ax_r, 'YDir', 'normal')
set(ax_r, 'Color', 'Black')
colormap(gca, 'hsv')
box off
title('RHO')

% Plot eta map
ax_e = nexttile(t_layout, 7);
imagesc(x_bin, y_bin, rad2deg(img_eta), 'AlphaData', ~isnan(img_eta) .* a_data);
axis image
set(ax_e, 'YDir', 'normal')
set(ax_e, 'Color', 'Black')
colormap(gca, 'parula')
box off
title('ETA')

% Plot delta map
ax_d = nexttile(t_layout, 8);
imagesc(x_bin, y_bin, rad2deg(img_delta), 'AlphaData', ~isnan(img_delta) .* a_data);
axis image
set(ax_d, 'YDir', 'normal')
set(ax_d, 'Color', 'Black')
colormap(gca, 'jet')
box off
title('DELTA')

% Link all axes for synchronized zooming and panning
linkaxes(t_layout.Children, 'xy')

% Add colorbars to all parameter plots
colorbar(ax_r);
colorbar(ax_e);
colorbar(ax_d);
colorbar(ax_m);

end
