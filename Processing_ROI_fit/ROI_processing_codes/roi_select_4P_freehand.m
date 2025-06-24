% =========================================================================
% Function Name: roi_select_4P_freehand.m
%
% Description:
%   This function allows freehand ROI (Region of Interest) selection
%   on a 4polar3D dataset and updates the ROI structure accordingly.
%
%   The function:
%     - Lets the user draw a freehand ROI on the intensity image.
%     - Updates both list and mask representations of the selected ROI.
%     - Adjusts the transparency (AlphaData) of multiple parameter maps
%       to highlight the selected region.
%     - Updates the summary mask displayed in the multi-ROI axis.
%
% Instructions:
%   - Call this function during the interactive exploration step of the 4polar3D pipeline.
%   - Draw a freehand shape on the image and double-click to finish.
%
% Notes:
%   - This function uses MATLAB's drawfreehand tool.
%   - The input axis handles must be linked to the corresponding parameter maps.
%
% Inputs:
%   - roi_struct: Structure containing ROI information (binning, masks, lists).
%   - ax_i: Axis handle for the intensity image.
%   - ax_r: Axis handle for the rho image.
%   - ax_e: Axis handle for the eta image.
%   - ax_d: Axis handle for the delta image.
%   - ax_m: Axis handle for the multi-ROI mask display.
%
% Outputs:
%   - roi_struct: Updated structure with new ROI list and mask.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function roi_struct = roi_select_4P_freehand(roi_struct, ax_i, ax_r, ax_e, ax_d, ax_m)

% Activate the ax_i image axis for ROI selection
axes(ax_i)

% Let the user draw a freehand ROI
%     temp_roi = roipoly;
%     temp = imfreehand;
temp = drawfreehand; % Use modern freehand drawing tool
temp.wait; % Wait for user to complete the selection

temp_roi = temp.createMask; % Generate binary mask from selection
temp.delete % Remove the ROI object after selection

% Initialize selection list for the current ROI
temp_list = zeros(size(roi_struct.Xbinned));

% Check each binned point to see if it falls inside the selected ROI
for k = 1:numel(roi_struct.Xbinned)
    temp_list(k) = temp_roi(roi_struct.Ybinned(k), roi_struct.Xbinned(k));
end

% Add the new selection to the list and mask forms in the ROI structure
roi_struct.list_form(:, end + 1) = temp_list;
roi_struct.mask_form(:, :, end + 1) = temp_roi;

% Adjust display transparency for the intensity image
temp = ax_i.Children.CData;
a_data = min(temp, prctile(temp(:), 99.975)); % Limit extreme intensity values
a_data = a_data ./ max(a_data);               % Normalize

% Update transparency (AlphaData) to highlight the selected ROI
ax_i.Children.AlphaData = ones(size(temp)) - temp_roi .* 0.5;

% Update transparency for rho, eta, and delta images
ax_r.Children.AlphaData = (~isnan(ax_r.Children.CData)) .* (a_data .* 0.4 + temp_roi .* 0.6);
ax_e.Children.AlphaData = (~isnan(ax_e.Children.CData)) .* (a_data .* 0.4 + temp_roi .* 0.6);
ax_d.Children.AlphaData = (~isnan(ax_d.Children.CData)) .* (a_data .* 0.4 + temp_roi .* 0.6);

% Update the combined mask display with all selected ROIs
all_masks = sum(roi_struct.mask_form(:, :, 2:end) .* cumsum(ones(1, 1, size(roi_struct.list_form, 2) - 1), 3), 3);
ax_m.Children.CData = all_masks;

%     Optional: Alpha blending control for ROI mask
%     img_obj.AlphaData = msk_roi + (msk_roi == 0) .* 0.25;

end
