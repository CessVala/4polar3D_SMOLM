% =========================================================================
% Function Name: mesure_intensite_parametres_connus.m
%
% Description:
%   This function measures the intensity of a particle at a known position (i0, j0)
%   and with a known PSF width (r0).
%
%   The function:
%     - Extracts a local window centered on the provided position.
%     - Computes the particle intensity using maximum likelihood estimation.
%     - Estimates the background noise variance.
%     - Ignores particles that are too close to the image boundaries.
%
% Inputs:
%   - r0: Known PSF radius.
%   - i0: Particle position.
%   - j0: Particle position.
%   - wn: Window size for calculation.
%   - im_input: Input fluorescence image.
%
% Outputs:
%   - alpha: Measured particle intensity at (i0, j0) for PSF width r0.
%   - sig2: Estimated background noise variance.
%   - reste: Residual image after particle subtraction.
%
% Notes:
%   - This version includes management of undetected particles in horizontal or vertical directions.
%   - If the particle is near the image border, the function returns zero intensity.
%   - If the estimated SNR is less than 0 dB, alpha is set to zero.
%
% Author:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2012
% =========================================================================

function [alpha, sig2, reste] = mesure_intensite_parametres_connus(r0, i0, j0, wn, im_input)

[idim, jdim] = size(im_input);

% Check if the particle is too close to the border
wn2 = floor(wn / 2);
if ((i0 < wn2) || (j0 < wn2) || (i0 > (idim - wn2)) || (j0 > (jdim - wn2)))
    alpha = 0;
    sig2 = 0;
    reste = 0;
    return;
end

% Extract the calculation window
round_i0 = round(i0);
delta_i0 = i0 - round_i0;
round_j0 = round(j0);
delta_j0 = j0 - round_j0;

di = (1:wn) + round_i0 - wn2;
dj = (1:wn) + round_j0 - wn2;
x = im_input(di, dj);

% Create the reference Gaussian
N = wn * wn;
refi = (0:(wn - 1)) - wn / 2;
refj = (0:(wn - 1)) - wn / 2;

i = 0.5 + refi - delta_i0;
j = 0.5 + refj - delta_j0;
ii = i' * ones(1, wn);
jj = ones(wn, 1) * j;

% Gaussian normalized to a power of one
g = (1 / (sqrt(pi) * r0)) * exp(-(1 / (2 * r0^2)) * (ii .* ii + jj .* jj));

% Maximum likelihood estimate of alpha
gc = g - sum(g(:)) / N;
Sgc2 = sum(gc(:).^2);
alpha = sum(sum(x .* gc)) / Sgc2;

% Estimate the background mean
x_alphag = x - alpha .* g;
m = sum(sum(x_alphag)) / N;

% Compute residual and noise variance
reste = x_alphag - m;
sig2 = var(reste(:));

% If SNR < 0 dB, discard the measurement
if (alpha < sqrt(sig2))
    fprintf(display_out, 'SNR < 0 dB, intensity not measured! (=> alpha = 0)\n');
    alpha = 0;
end

end %function
