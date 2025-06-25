% =========================================================================
% Function Name: deplt_GN_estimation_ssr.m
%
% Description:
%   This function performs Gauss-Newton optimization to estimate the 2D position
%   of a Gaussian signal without optimizing the Gaussian width (sigma is fixed).
%
%   The function:
%     - Refines the position estimates of the Gaussian.
%     - Minimizes the variance of the fitting error.
%     - Iteratively reduces displacement steps if the criterion increases.
%
% Instructions:
%   - Provide the fixed Gaussian width, the initial position, the image window,
%     and optionally previous displacements.
%
% Inputs:
%   - r0: Fixed Gaussian width (sigma).
%   - p_i: Previous i position.
%   - p_j: Previous j position.
%   - x: Image window.
%   - sig2init: Initial noise variance.
%   - p_di: Previous displacement in i.
%   - p_dj: Previous displacement in j.
%
% Outputs:
%   - n_i: Updated i position.
%   - n_j: Updated j position.
%   - di: Displacement in i.
%   - dj: Displacement in j.
%   - alpha: Estimated amplitude.
%   - sig2: Estimated noise variance.
%   - reste: Residual error.
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: July 2010
% =========================================================================

function [n_i, n_j, di, dj, alpha, sig2, reste] = deplt_GN_estimation_ssr(r0, p_i, p_j, x, sig2init, p_di, p_dj)

%% p_di, p_dj: previous displacements that led to p_i, p_j

i0 = p_i;
j0 = p_j;
prec_rel = 0.001;

if (nargin < 7)
    verif_crit = 0;
    sig2init = inf;
else
    verif_crit = 1;
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
        p_di = p_di / 10.0;
        p_dj = p_dj / 10.0;
        i0 = pp_i + p_di;
        j0 = pp_j + p_dj;
        if (max(abs(p_di), abs(p_dj)) < prec_rel)
            n_i = p_i;
            n_j = p_j;
            di = 0;
            dj = 0;
            fprintf(display_out, 'convergence issue: EXIT !!!\n');
            return;
        end
    else
        again = 0;
    end

end

%% Compute gradients
der_g_i0 = ii .* g_div_sq_r0;
der_g_j0 = jj .* g_div_sq_r0;

%% Compute second derivatives
derder_g_i0 = (-1 + inv(r0^2) * iiii) .* g_div_sq_r0;
derder_g_j0 = (-1 + inv(r0^2) * jjjj) .* g_div_sq_r0;

%% First derivative of cost function /2
der_J_i0 = alpha * sum(sum(der_g_i0 .* err));
der_J_j0 = alpha * sum(sum(der_g_j0 .* err));

%% Second derivative of cost function /2
derder_J_i0 = alpha * sum(sum(derder_g_i0 .* err)) - alpha^2 * sum(sum(der_g_i0 .^ 2));
derder_J_j0 = alpha * sum(sum(derder_g_j0 .* err)) - alpha^2 * sum(sum(der_g_j0 .^ 2));

%% Gauss-Newton displacement
di = - der_J_i0 / derder_J_i0;
dj = - der_J_j0 / derder_J_j0;

n_i = i0 + di;
n_j = j0 + dj;

end %function
