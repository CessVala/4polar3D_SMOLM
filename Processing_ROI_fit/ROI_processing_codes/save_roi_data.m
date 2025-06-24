% =========================================================================
% Function Name: save_roi_data.m
%
% Description:
%   This function saves the extracted 4polar3D ROI data into an Excel table.
%
%   The function:
%     - Allows selecting specific ROIs to export, or exports all by default.
%     - Combines data from multiple ROIs into a single table.
%     - Saves the table as an Excel file via a user prompt.
%
% Instructions:
%   - Call this function to export ROI data to Excel.
%   - Optional: Provide a list of ROI indices to export (zero-based).
%
% Notes:
%   - The ROI indices provided as input should not exceed the number of ROIs.
%   - Pixel size is set to 1 (no scaling).
%
% Inputs:
%   - roi_data: Structure array containing the extracted ROI data.
%   - varargin: (Optional) List of ROI indices to export (zero-based).
%
% Outputs:
%   - None. The function saves the file and opens it automatically.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function save_roi_data(roi_data, varargin)

% Check if specific ROIs are provided, otherwise export all
if nargin > 1
    roi_num = varargin{1};
else
    roi_num = 1:(numel(roi_data) - 1);
end

disp(['Saving ' num2str(numel(roi_num)) ' out of ' num2str(numel(roi_data) - 1) ' ROIs...'])

% Validate that the input is a structure
if ~isstruct(roi_data)
    errordlg('ROI_data doesnt exist')
    return
    % Validate that the ROI indices are within range
elseif max(roi_num) > (numel(roi_data) - 1)
    errordlg('ROI number is larger than total number of ROIs')
    return
end

% Pixel size (set to 1 if no scaling is required)
% px_size = 0.130; % Example of real scaling
px_size = 1;

% Initialize arrays to collect data
Rho = [];
Eta = [];
Delta = [];
X = [];
Y = [];
I_total = [];
Frame = [];
ROI_number = [];

% Loop over selected ROIs and concatenate their data
for k = 1:numel(roi_num)
    temp_num = roi_num(k) + 1;

    Rho = cat(1, Rho, rad2deg(roi_data(temp_num).rho));
    Eta = cat(1, Eta, rad2deg(roi_data(temp_num).eta));
    Delta = cat(1, Delta, rad2deg(roi_data(temp_num).delta));
    I_total = cat(1, I_total, roi_data(temp_num).I_tot);
    Frame = cat(1, Frame, roi_data(temp_num).frame);
    X = cat(1, X, roi_data(temp_num).X .* px_size);
    Y = cat(1, Y, roi_data(temp_num).Y .* px_size);
    ROI_number = cat(1, ROI_number, roi_num(k) * ones(size(roi_data(temp_num).rho)));
end

% Create table with all collected data
data_table = table(ROI_number, Rho, Eta, Delta, X, Y, I_total, Frame);

% Prompt user to select a file name and save location
[fname, fpath] = uiputfile({'*.xlsx'; '*.xls'});

% Write table to Excel file
writetable(data_table, fullfile(fpath, fname));

% Automatically open the saved file
winopen(fullfile(fpath, fname))
end
