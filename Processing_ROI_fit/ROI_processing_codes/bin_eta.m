% =========================================================================
% Function Name: bin_eta.m
%
% Description:
%   This function bins the eta values into quantile-based categories and
%   visualizes them using the 4polar3D stick plot.
%
%   The function:
%     - Selects eta values corresponding to the provided ROI list.
%     - Bins eta values into four quantile-based groups.
%     - Assigns these bins as color indices for the sticks.
%     - Visualizes the result using the create_single_4Psticks_image function.
%
% Instructions:
%   - Call this function with the eta array, ROI indices, the input image,
%     coordinate and display parameters.
%   - The function automatically bins eta and updates the stick colors accordingly.
%
% Notes:
%   - This function is part of the 4polar3D visualization pipeline.
%   - The colorbar is set to display four bins.
%   - The 'plasma' colormap is used by default.
%   - Example figure export commands are provided at the end (commented).
%
% Inputs:
%   - eta: Array of eta values.
%   - roi_list: Indices of the regions of interest.
%   - img_strm: Background image for plotting.
%   - xd, yd: Coordinates for the stick representation.
%   - x_loc, y_loc: X and Y positions of the features.
%   - rho: Orientation or intensity parameter.
%   - c_line: Color assignment vector (updated inside the function).
%   - sl: Stick length parameter.
%   - sw: Stick width parameter.
%   - sa: Stick alpha (transparency) parameter.
%
% Outputs:
%   - ax_img_e0: Axis handle for the background image.
%   - ax_stk_e0: Axis handle for the stick plot.
%   - fig_e0: Figure handle of the generated plot.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function [ax_img_e0, ax_stk_e0, fig_e0] = bin_eta(eta, roi_list, img_strm, xd, yd, x_loc, y_loc, rho, c_line, sl, sw, sa)

% Number of bins to split eta into
nn = 4;

% Bin eta values for selected ROIs into nn bins using percentiles
[~,~,binned_eta] = histcounts(eta(roi_list), prctile(eta(roi_list), linspace(0,1,nn+1)*100));

% Assign binned eta values as colors for the sticks
c_line = binned_eta;

% Create figure for plotting
fig_e0 = figure('Position', [50 50 940 700], 'Resize', 'off');

% Generate the 4polar3D stick plot with binned colors
[ax_img_e0, ax_stk_e0, fig_e0] = create_single_4Psticks_image(fig_e0, img_strm, xd, yd, ...
    x_loc(roi_list), ...
    y_loc(roi_list), ...
    rho(roi_list), ...
    c_line, sl, sw, sa);

% Set colormap to plasma for visualization
colormap(ax_stk_e0, "plasma");
% Alternative colormap (commented)
% colormap(ax_stk_e0, "parula")

% Set plot title
title('\eta binned')

% Customize tick marks of the color axis (assumes first child is colorbar)
temp = get(gcf, 'Children');
temp(1).Ticks = 1:nn;

% Set figure background color to white
set(fig_e0, 'Color', 'w')
end

% % Example figure export commands:
% print([fld 'C1_eta_binned'],'-vector','-dpdf','-r1500')
% print([fld 'C1_eta_binned3'],'-dpng','-r1500')
