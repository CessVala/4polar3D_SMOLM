% =========================================================================
% Function Name: all_max_2d.m
%
% Description:
%   This function generates a binary map of all local maxima in a 2D matrix.
%
%   The function:
%     - Compares each pixel to its 8 neighbors.
%     - Identifies positions where the center pixel is larger than its
%       neighbors in horizontal, vertical, and diagonal directions.
%     - Returns a binary map where local maxima retain their original value
%       and all other positions are zero.
%
% Instructions:
%   - Call this function with a 2D input matrix.
%
% Notes:
%   - This is a simple local maxima detector without border processing.
%   - The output map contains the local maxima intensities.
%
% Inputs:
%   - input: 2D numerical array.
%
% Outputs:
%   - carte_max: 2D array with local maxima positions retaining original values.
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function carte_max = all_max_2d(input)

% Get input size
[N, M] = size(input);

% Extract the inner matrix (excluding borders)
ref = input(2:N-1, 2:M-1);

% Check local maxima: compare with top and bottom neighbors
pos_max_h = input(1:N-2, 2:M-1) < ref & ...
    input(3:N  , 2:M-1) < ref;

% Check local maxima: compare with left and right neighbors
pos_max_v = input(2:N-1, 1:M-2) < ref & ...
    input(2:N-1, 3:M  ) < ref;

% Check local maxima: compare with 135-degree diagonal neighbors
pos_max_135 = input(1:N-2, 1:M-2) < ref & ...
    input(3:N  , 3:M  ) < ref;

% Check local maxima: compare with 45-degree diagonal neighbors
pos_max_45  = input(3:N  , 1:M-2) < ref & ...
    input(1:N-2, 3:M  ) < ref;

% Initialize the output map
carte_max = zeros(N, M);

% Set the local maxima positions to their original intensity values
carte_max(2:N-1, 2:M-1) = pos_max_h & pos_max_v & pos_max_135 & pos_max_45;
carte_max = carte_max .* input;

end %function
