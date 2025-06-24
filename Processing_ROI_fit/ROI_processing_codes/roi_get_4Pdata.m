% =========================================================================
% Function Name: roi_get_4Pdata.m
%
% Description:
%   This function extracts the localization and parameter data corresponding
%   to a selected ROI from the 4polar3D dataset.
%
%   The function:
%     - Selects the subset of data points belonging to the specified ROI.
%     - Retrieves the parameters for those points.
%
% Instructions:
%   - Call this function to extract the data for a specific ROI index.
%   - ROI indices start from 0 in the mask_list structure.
%
% Notes:
%   - This function is part of the 4polar3D MATLAB analysis pipeline.
%   - The ROI selection is based on the binary mask list.
%
% Inputs:
%   - roi_num: Index of the ROI to extract (0-based index).
%   - data: Structure containing the full localization and parameter dataset.
%   - mask_list: List-form matrix indicating ROI membership for each localization.
%
% Outputs:
%   - roi_data: Structure containing the extracted data for the selected ROI.
%       .rho:
%       .eta:
%       .delta:
%       .rmse:   Localization fitting error.
%       .frame:  Frame index of each localization.
%       .X, .Y:  Localization positions.
%       .I_tot:  Total intensity.
%       .I:      Per-channel intensity data.
%       .noise:  Per-channel noise estimation.
%       .loc_precision: Per-channel localization precision.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function roi_data = roi_get_4Pdata(roi_num, data, mask_list)


m_lst = mask_list(:, roi_num + 1) == 1;

% Extract rho values for the selected ROI
roi_data.rho = data.rho(m_lst);

% Extract eta values for the selected ROI
roi_data.eta = data.eta(m_lst);

% Extract delta values for the selected ROI
roi_data.delta = data.delta(m_lst);

% Extract RMSE for the selected ROI
roi_data.rmse = data.rmse(m_lst);

% Extract frame indices for the selected ROI
roi_data.frame = data.frame(m_lst);

% Extract localization X and Y positions for the selected ROI
roi_data.X = data.x_loc(m_lst);
roi_data.Y = data.y_loc(m_lst);

% Extract total intensity for each localization in the selected ROI
roi_data.I_tot = data.I_tot(m_lst);

% Extract per-channel intensity data
roi_data.I = data.I(m_lst, :);

% Extract per-channel noise estimation
roi_data.noise = data.noise(m_lst, :);

% Extract per-channel localization precision
roi_data.loc_precision = data.loc_precision(m_lst, :);

end
