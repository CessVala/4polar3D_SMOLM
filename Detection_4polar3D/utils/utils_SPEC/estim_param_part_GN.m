% =========================================================================
% Function Name: estim_param_part_GN.m
%
% Description:
%   This function estimates the parameters of a particle using a Gauss-Newton
%   optimization, either with or without optimization of the Gaussian width.
%
%   The function:
%     - Extracts a local image window around the initial particle position.
%     - Performs iterative position and optional radius refinement.
%     - Stops when convergence is achieved or bounds are exceeded.
%
% Instructions:
%   - Call this function with the image, window size, and initial detection parameters.
%
% Notes:
%   - The window size 'wn' must be odd.
%   - Default displacement bounds and radius limits are used if not provided.
%
% Inputs:
%   - im: Input image.
%   - wn: Window size (must be odd).
%   - liste_info_param: detection results.
%   - r0: Initial Gaussian width (sigma).
%   - flag_r: Flag to enable/disable optimization of Gaussian width.
%   - bornes_ijr: Bounds for displacements and radius [i_min, i_max, j_min, j_max, r_min, r_max].
%
% Outputs:
%   - liste_param: Refined parameter list.
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: July 7, 2010
% =========================================================================

function liste_param = estim_param_part_GN(im, wn, liste_info_param, r0, flag_r, bornes_ijr)

%% If no displacement bounds provided, use default limits
if (nargin < 6)
    bornes_ijr(1) = -1.5;
    bornes_ijr(2) = 1.5;
    bornes_ijr(3) = -1.5;
    bornes_ijr(4) = 1.5;
    bornes_ijr(5) = 0.3;
    bornes_ijr(6) = 3.0;
    if (nargin < 5)
        flag_r = 0;
    end
end

%% Extract local window around initial position
Pi = liste_info_param(2);
Pj = liste_info_param(3);
di = (1:wn) + Pi - floor(wn / 2);
dj = (1:wn) + Pj - floor(wn / 2);
im_part = im(di, dj);

%% Initialization
r = r0;
i = 0.0;
j = 0.0;
dr = 1;
di = 1;
dj = 1;
fin = 0.005;
sig2 = inf;
cpt = 0;
test = 1;
ITER_MAX = 50;

%% Iterative Gauss-Newton estimation
while (test)

    if (flag_r)
        [r, i, j, dr, di, dj, alpha, sig2] = deplt_GN_estimation(r, i, j, im_part, sig2, dr, di, dj);
    else
        [i, j, di, dj, alpha, sig2] = deplt_GN_estimation_ssr(r0, i, j, im_part, sig2, di, dj);
    end

    cpt = cpt + 1;

    if (flag_r)
        test = max([abs(di), abs(dj), abs(dr)]) > fin;
    else
        test = max(abs(di), abs(dj)) > fin;
    end

    if (cpt > ITER_MAX)
        test = 0;
    end

    %% Stop if parameters exceed bounds
    result_ok = ~((i < bornes_ijr(1)) || (i > bornes_ijr(2)) || ...
        (j < bornes_ijr(3)) || (j > bornes_ijr(4)) || ...
        (r < bornes_ijr(5)) || (r > bornes_ijr(6)));

    test = test & result_ok;

end

%% Output structure
% liste_info_param = [num, i, j, alpha, sig^2]
% liste_param = [num, i, j, alpha, sig^2, rayon, ok]

liste_param = [liste_info_param(1), ...
    Pi + i, ...
    Pj + j, ...
    alpha, ...
    sig2, ...
    r, ...
    result_ok];

end %function
