% =========================================================================
% Function Name: detect_et_estime_part_1vue.m
%
% Description:
%   This function detects and estimates particle parameters.
%
%   The function:
%     - Detects particles using the hypothesis test map.
%     - Refines particle parameters using Gauss-Newton optimization.
%     - Applies border exclusion to avoid edge artifacts.
%
% Instructions:
%   - Call this function with the input image and optionally specify detection window size,
%     initial radius, PFA threshold, and search step sizes.
%
% Inputs:
%   - input: Input image.
%   - wn: Detection and estimation window size (default: wn = 9).
%   - r0: Initial Gaussian radius for detection (default: r0 = 1.3).
%   - pfa: Detection threshold (default: pfa = 28, see carte_H0H1_1vue.m).
%   - pas_ijr: (Not used in this version) Optional search step sizes.
%
% Outputs:
%   - ldetect: Detected particles [index, i, j, amplitude, sig^2].
%   - lestime: Estimated parameters [index, i, j, amplitude, sig^2, radius, validity flag].
%   - d: Hypothesis map (output from carte_H0H1_1vue).
%   - Nestime: Number of estimated particles.
%
% Notes:
%   - The half-height radius is given by: radius_mih = sqrt(2 * ln(2)) * radius.
%   - alpha: Amplitude of the Gaussian normalized to a power of one. Signal power: P = alpha^2.
%   - Maximum amplitude of the Gaussian signal: alpha_max = alpha / (sqrt(pi) * radius).
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% =========================================================================

function [lestime, ldetect, d, Nestime] = detect_et_estime_part_1vue(input, wn, r0, pfa, pas_ijr)

if (nargin < 2)
    wn = 9;
end
if (nargin < 3)
    r0 = 1.3;
end
if (nargin < 4)
    pfa = 28;
end
% if (nargin < 5)
%     pas_ijr = [0.125 0.125 0.08];
% end

%% Optimization along the radius axis
flag_r = 1; % 0/1

[Ni, Nj] = size(input);

%% Parameter positions in the list
Nparam = 7;
detect_i = 2;
detect_j = 3;
alpha = 4;
% sig2 = 5;

%% Detection with a mean radius r0
[c, ldetect, d] = carte_H0H1_1vue(input, r0, wn, wn, pfa);

%% Estimation with search steps
%% Search intervals: [-1, 1] for i and j, [0.6, 1.4] for the radius

Ndetect = size(ldetect, 1);
if (Ndetect == 0)
    lestime = zeros(1, Nparam);
    fprintf(display_out, 'No particle detected => exit\n');
    Nestime = 0;
    return;
end

Nestime = 0;
bord = ceil(wn / 2);
lestime = zeros(Ndetect, Nparam);

fprintf(display_out, 'Number of detected particles: %d\n', Ndetect);
for n = 1:Ndetect
    test_bord = (ldetect(n, detect_i) < bord) || (ldetect(n, detect_i) > Ni - bord) || ...
        (ldetect(n, detect_j) < bord) || (ldetect(n, detect_j) > Nj - bord);
    if ((ldetect(n, alpha) > 0.0) && (~test_bord))
        Nestime = Nestime + 1;
        lestime(Nestime, :) = estim_param_part_GN(input, wn, ldetect(n, :), r0, flag_r);
    end
end

% Resize to the correct size
if (Nestime == 0)
    lestime = zeros(1, Nparam);
else
    lestime = lestime(1:Nestime, :);
end
fprintf(display_out, 'Particles detected: %d, estimated: %d, valid: %d\n', Ndetect, Nestime, sum(lestime(:, 7)));

end %function
