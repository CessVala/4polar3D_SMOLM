% =========================================================================
% Function Name: BCR_1G_ij_rconnu_liste_r.m
%
% Description:
%   This function calculates the Cramer-Rao lower bounds for the estimation
%   of particle positions in the case of a Gaussian (1G) with known radius.
%
%   The function:
%     - Computes the standard deviation (std_ij) using the Cramer-Rao lower bound.
%     - Assumes a symmetric PSF: std_i == std_j == std_ij.
%     - Processes a list of particles (input vectors for each parameter).
%
% Inputs:
%   - wn: Window size used for parameter estimation.
%   - list_sig2: Vector of noise variances for each particle.
%   - list_r0: Vector of known Gaussian radii.
%   - list_alpha: Vector of estimated Gaussian amplitudes.
%
% Outputs:
%   - var_ij: Vector of variances (standard deviations squared) for each particle position.
%
% Notes:
%   - The PSF is assumed to be symmetric.
%   - All input lists must be row vectors.
%
% Author:
%   Nicolas Bertaux - Institut Fresnel
%
% Date: December 2012
% =========================================================================

function var_ij = BCR_1G_ij_rconnu_liste_r(wn, list_sig2, list_r0, list_alpha)

% Define window coordinates centered at zero
i = 0.5 + (0:wn-1) - wn/2;
ii = i' * ones(1, wn);

% Initialize output vector
N = length(list_alpha);
var_ij = zeros(1, N);

% Loop through each particle
for n = 1:N

    r0 = list_r0(n);         % Known Gaussian radius
    sig2 = list_sig2(n);     % Noise variance
    alpha = list_alpha(n);   % Estimated Gaussian amplitude

    % Generate the Gaussian
    G = gausswin2(r0, wn, wn, 0, 0);

    % Compute partial derivatives of G
    derG_i0 = ii .* G / (r0^2);


    CTE = sum(sum(derG_i0.^2));
    F_ii = ((alpha^2) / sig2) * CTE;

    % Cramer-Rao variance
    var_ij(n) = 1 / (F_ii);

end

end %function
