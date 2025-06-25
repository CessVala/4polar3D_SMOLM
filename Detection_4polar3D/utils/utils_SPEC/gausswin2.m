% =========================================================================
% Function Name: gausswin2.m
%
% Description:
%   This function generates a 2D Gaussian with a specified width (sigma)
%   of power 1, optionally offset in X and Y directions.
%
%   The function:
%     - Creates a centered or shifted 2D Gaussian.
%     - Supports both square and rectangular windows.
%     - Provides optional X and Y pixel offsets.
%
% Instructions:
%   - Call this function to generate a Gaussian of power 1.
%
% Examples:
%   - g = gausswin2(1.3, 8);                  % Square window, no offset
%   - g = gausswin2(1.3, 8, 0.3, -0.8);       % Square window with offset
%   - g = gausswin2(1.3, 8, 9, 0.3, -0.8);    % Rectangular window with offset
%

% Inputs:
%   - sig: Gaussian width (sigma).
%   - wn_i: Window size along the i axis.
%   - wn_j: Window size along the j axis (optional, defaults to wn_i).
%   - offset_i: Optional offset in pixels along i axis.
%   - offset_j: Optional offset in pixels along j axis.
%
% Outputs:
%   - g: Generated 2D Gaussian window.
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: December 2017
% =========================================================================

function g = gausswin2(sig, wn_i, wn_j, offset_i, offset_j)

if (nargin < 5)
    offset_i = 0.0;
    offset_j = 0.0;
end

if (nargin < 3)
    wn_j = wn_i;
end

% Create window axes
refi = 0.5 + (0:(wn_i - 1)) - wn_i / 2;
i = refi - offset_i;
refj = 0.5 + (0:(wn_j - 1)) - wn_j / 2;
j = refj - offset_j;

% Create grid
ii = i' * ones(1, wn_j);
jj = ones(wn_i, 1) * j;

% Generate Gaussian with power 1
g = (1 / (sqrt(pi) * sig)) * exp(-(1 / (2 * sig^2)) * (ii .* ii + jj .* jj));

end %function
