% =========================================================================
% Function Name: deplt_GN_estimation.m
%
% Description:
%   This function performs Gauss-Newton optimization to estimate the Gaussian parameters
%   (radius, position, amplitude, and noise variance) that best fit a particle signal.
%
%   The function:
%     - Refines the Gaussian width (radius) and position estimates.
%     - Minimizes the variance of the fitting error.
%     - Adjusts the displacement step if the criterion increases.
%
% Instructions:
%   - Provide the previous Gaussian parameters, image window, and optionally previous displacements.
%
% Inputs:
%   - p_r: Previous Gaussian width (sigma).
%   - p_i: Previous i position.
%   - p_j: Previous j position.
%   - x: Image window.
%   - sig2init: Initial noise variance.
%   - p_dr: Previous displacement in r.
%   - p_di: Previous displacement in i.
%   - p_dj: Previous displacement in j.
%
% Outputs:
%   - n_r: Updated Gaussian width.
%   - n_i: Updated i position.
%   - n_j: Updated j position.
%   - dr: Displacement in Gaussian width.
%   - di: Displacement in i position.
%   - dj: Displacement in j position.
%   - alpha: Estimated amplitude.
%   - sig2: Estimated noise variance.
%   - reste: Residual error.
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS (modifications)
%
% Date: December 2017
% =========================================================================

function [n_r, n_i, n_j, dr, di, dj, alpha, sig2, reste] = deplt_GN_estimation(p_r, p_i, p_j, x, sig2init, p_dr, p_di, p_dj)

%% p_dr, p_di, p_dj: previous displacements that led to p_r, p_i, p_j

r0 = p_r;
i0 = p_i;
j0 = p_j;
prec_rel = 0.001;

if (nargin < 7)
    verif_crit = 0;
    sig2init = inf;
else
    verif_crit = 1;
    pp_r = r0 - p_dr;
    pp_i = i0 - p_di;
    pp_j = j0 - p_dj;
end

[wn_i, wn_j] = size(x);
N = wn_i * wn_j;
refi = 0.5 + (0:(wn_i - 1)) - wn_i / 2;
refj = 0.5 + (0:(wn_j - 1)) - wn_j / 2;

%% Iteratively reduce displacements as long as the new criterion increases
again = 1;
while (again)

    i = refi - i0;
    j = refj - j0;
    ii = i' * ones(1, wn_j);
    jj = ones(wn_i, 1) * j;

    %% Gaussian normalized to a power of one
    iiii = ii .* ii;
    jjjj = jj .* jj;
    iiii_jjjj = iiii + jjjj;
    g = (1 / (sqrt(pi) * r0)) * exp(-(1 / (2 * r0^2)) * (iiii_jjjj));
    gc = g - sum(g(:)) / N;
    Sgc2 = sum(gc(:).^2);
    g_div_sq_r0 = inv(r0^2) * g;

    %% Maximum likelihood estimate of alpha
    if (Sgc2 ~= 0)
        alpha = sum(sum(x .* gc)) / Sgc2;
    else
        alpha = 0;
    end

    %% Maximum likelihood estimate of m
    x_alphag = x - alpha .* g;
    m = sum(sum(x_alphag)) / N;

    err = x_alphag - m;

    %% For information
    reste = err;

    %% Criterion before displacement
    sig2 = sum(sum(err.^2)) / N;
    if ((verif_crit) && (sig2 > sig2init))
        p_dr = p_dr / 10.0;
        p_di = p_di / 10.0;
        p_dj = p_dj / 10.0;

        r0 = pp_r + p_dr;
        i0 = pp_i + p_di;
        j0 = pp_j + p_dj;
        if (max([abs(p_dr), abs(p_di), abs(p_dj)]) < prec_rel)
            n_r = p_r;
            n_i = p_i;
            n_j = p_j;
            dr = 0;
            di = 0;
            dj = 0;
            fprintf(display_out, 'EXIT !!!\n');
            return;
        end
    else
        again = 0;
    end

end

%% Compute gradients
der_g_i0 = ii .* g_div_sq_r0;
der_g_j0 = jj .* g_div_sq_r0;
der_g_r0 = (-inv(r0) + inv(r0^3) * (iiii_jjjj)) .* g;

%% Compute second derivatives
derder_g_i0 = (-1 + inv(r0^2) * iiii) .* g_div_sq_r0;
derder_g_j0 = (-1 + inv(r0^2) * jjjj) .* g_div_sq_r0;
derder_g_r0 = (1 - 3 * inv(r0^2) * iiii_jjjj) .* g_div_sq_r0 ...
    + (-inv(r0) + inv(r0^3) * (iiii_jjjj)) .* der_g_r0;

%% First derivative of cost function /2
der_J_i0 = alpha * sum(sum(der_g_i0 .* err));
der_J_j0 = alpha * sum(sum(der_g_j0 .* err));
der_J_r0 = alpha * sum(sum(der_g_r0 .* err));

%% Second derivative of cost function /2
derder_J_i0 = alpha * sum(sum(derder_g_i0 .* err)) - alpha^2 * sum(sum(der_g_i0 .^ 2));
derder_J_j0 = alpha * sum(sum(derder_g_j0 .* err)) - alpha^2 * sum(sum(der_g_j0 .^ 2));
derder_J_r0 = alpha * sum(sum(derder_g_r0 .* err)) - alpha^2 * sum(sum(der_g_r0 .^ 2));

%% Gauss-Newton displacement
dr = - der_J_r0 / derder_J_r0;
di = - der_J_i0 / derder_J_i0;
dj = - der_J_j0 / derder_J_j0;

n_r = abs(r0 + dr); %% r0 > 0
n_i = i0 + di;
n_j = j0 + dj;

end %function
