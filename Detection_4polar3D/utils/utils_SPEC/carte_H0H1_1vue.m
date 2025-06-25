% =========================================================================
% Function Name: carte_H0H1_1vue.m
%
% Description:
%   This function computes a hypothesis test map (H0: no particle / H1: particle present)
%   assuming a Gaussian background noise (iid) and a particle signal modeled as a Gaussian.
%   
%   The function:
%     - Assumes the particle amplitude and noise variance are unknown.
%     - Computes the likelihood ratio map for detecting particles at the center of a scanning window.
%     - Uses a user-specified false alarm rate (PFA) threshold.
%     - Detects local maxima and returns the detected particle list.
%
% Instructions:
%   - The window size should preferably be less than 12 pixels.
%   - A smaller window is preferable as the background is assumed to be constant.
%
% Inputs:
%   - im: Input image.
%   - rayon: Gaussian width.
%   - wn_x: Window size in X direction.
%   - wn_y: Window size in Y direction (square window if omitted).
%   - s_pfa: Detection threshold corresponding to a false alarm rate.
%
% Outputs:
%   - carte_MV: Hypothesis test map.
%   - liste_detect: List of detected particles [index, i, j, amplitude, local noise variance, particle radius, flag].
%   - detect_pfa: Binary detection mask.
%
% Notes:
%   - By default, the detection threshold s_pfa corresponds to PFA = 1e-7.
%   - Recommended to use an even window size if the image size is even in that direction.
%   - The detection uses a chi-square distribution with 1 degree of freedom (3 - 2).
%
%   - Recommended thresholds for s_pfa and their corresponding detection probabilities:
%
%       s_pfa       1 - PFA
%      3.77060    0.94784
%      6.63106    0.98998
%     10.79172    0.99898
%     15.00000    0.999892488823270
%     20.00000    0.999992255783569
%     24.00000    0.999999036642991  (1E-6)
%     25.00000    0.999999426696856
%     28.00000    0.999999878684549  (1E-7)
%     30.00000    0.999999956795369
%     33.00000    0.999999990784113  (1E-8)
%     37.50000    0.999999999085870  (1E-9)
%     40.00000    0.999999999746037
%
%   - Example: [cmv, liste_detect] = carte_H0H1_1vue(im, 1.3, 8, 8, 33);
%   - By default (if no threshold is provided), the PFA corresponds to approximately 1E-7.


%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: December 2017
% =========================================================================

%% Hypothesis map: no particle vs particle detection at the window center
%% Assuming Gaussian iid background noise and a Gaussian particle signal.

function [carte_MV, liste_detect, detect_pfa] = carte_H0H1_1vue(im, rayon, wn_x, wn_y, s_pfa)

if (nargin < 5) 
  s_pfa = 28.0;
end
if (nargin < 4) 
  wn_y = wn_x;
end

[N, M] = size(im);
T = wn_x * wn_y; % Number of pixels in the window

%% Hypothesis H0: no particle in the window
m = ones(wn_x, wn_y);
hm = expand_w(m, N, M);
tfhm = fft2(hm);
tfim = fft2(im);
m0 = real(fftshift(ifft2(tfhm .* tfim))) / T;

im2 = im .* im;
tfim2 = fft2(im2);
Sim2 = real(fftshift(ifft2(tfhm .* tfim2)));

%% H0 = T/2*log(2*pi*sig0^2)-T/2 ;
T_sig0_2 = Sim2 - T * m0.^2;

%% Hypothesis H1: one particle at the window center, unknown amplitude, fixed radius
%% Generate Gaussian mask with width (sigma) equal to the radius

g = gausswin2(rayon, wn_x, wn_y, 0, 0);
gc = g - sum(g(:)) / T;
Sgc2 = sum(gc(:).^2);
hgc = expand_w(gc, N, M);
tfhgc = fft2(hgc);

alpha = real(fftshift(ifft2(tfhgc .* tfim))) / Sgc2;

%%  carte_MV = -0.5*(H0 - H1) ;
% carte_MV = - T * log(1 - (Sgc2 * alpha.^2) ./ T_sig0_2) ; 
test = 1 - (Sgc2 * alpha.^2) ./ T_sig0_2;
test = (test > 0) .* test + (test <= 0);
carte_MV = - T * log(test);

%% Detection and local maxima search
detect_masque = carte_MV > s_pfa;
if (sum(detect_masque(:)) == 0)
    fprintf(display_out, 'No target detected!\n');
    liste_detect = zeros(1, 6);
    detect_pfa = zeros(size(detect_masque));
else
    detect_pfa = all_max_2d(carte_MV) & detect_masque;

    [di, dj] = find(detect_pfa);
    n_detect = size(di, 1);
    vind = N * (dj - 1) + di;
    valpha = alpha(:);
    alpha_detect = valpha(vind);

    sig1_2 = (T_sig0_2 - alpha.^2 * Sgc2) / T;
    vsig1_2 = sig1_2(:);
    sig2_detect = vsig1_2(vind);

    liste_detect = [(1:n_detect)', di, dj, alpha_detect, sig2_detect, rayon * ones(n_detect, 1), ones(n_detect, 1)];
end
end %function
