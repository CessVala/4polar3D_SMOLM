% =========================================================================
% Function Name: estime_decalage_dg.m
%
% Description:
%   This function estimates the displacement (shift) between two sub-images
%   in 4polar3D data using lists of detected particles from both sub-images.
%
%   The function:
%     - Computes all possible displacements between right and left particle lists.
%     - Selects the most frequent displacement based on statistical thresholds.
%     - Returns the estimated displacement and the indices of the matched particle pairs.
%
% Inputs:
%   - vec_d_i: Row vector of particle positions along the i-axis (right image).
%   - vec_d_j: Row vector of particle positions along the j-axis (right image).
%   - vec_g_i: Row vector of particle positions along the i-axis (left image).
%   - vec_g_j: Row vector of particle positions along the j-axis (left image).
%   - sigma: Standard deviation of position errors on the i and j axes.
%   - dir_min: Minimum number of matching particle pairs required to validate the estimation.
%
% Outputs:
%   - decal_i: Estimated displacement along the i-axis.
%   - decal_j: Estimated displacement along the j-axis.
%   - nb_est: Number of particle pairs considered in the final estimation.
%   - indice_couple: Indices of the matched particle pairs used for estimation
%                    (vector of size: 2 * nb_est) formatted as [d1, g1, d2, g2, ...].
%
% Notes:
%   - This version assumes two provided lists of points (right and left).
%   - If you modify the geometry (other than right/left separation), update separe_liste.m accordingly.
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function [decal_i, decal_j, nb_est, indice_couple] = estime_decalage_dg(vec_d_i, vec_d_j, vec_g_i, vec_g_j, sigma, dir_min)

% If threshold not provided
if (nargin < 3)
    seuil = 0.5; % Defined by position estimation error
end
% If minimum number of pairs not provided
if (nargin < 4)
    dir_min = 3; % Minimum number of matched pairs
end

% Number of particles in right and left lists
nb_part_d = size(vec_d_i, 2);
nb_part_g = size(vec_g_i, 2);

% If one of the lists is empty, return zeros
if ((nb_part_d == 0) || (nb_part_g == 0))
    decal_i = 0;
    decal_j = 0;
    nb_est = 0;
    indice_couple = [0, 0];
    return;
end

% Compute displacement matrix between right and left particles
% Each row: particles from right list
% Each column: particles from left list

mat_di = repmat(vec_d_i.', 1, nb_part_g) - repmat(vec_g_i, nb_part_d, 1)  ;
mat_dj = repmat(vec_d_j.', 1, nb_part_g) - repmat(vec_g_j, nb_part_d, 1)  ;

nb_mat_dir = nb_part_d * nb_part_g;

% Represent displacement as complex numbers (i displacement + j displacement)
mat_dir = mat_di + 1i * mat_dj;
vec_dir = single(mat_dir(:)); % Flatten displacement matrix into a column vector

% Step 1: Search for similar displacement vectors to limit candidate directions
% Compute matrix of distances between all displacement vectors
vmodif = vec_dir.';
mat_dist_vec_dir = vec_dir(:, ones(nb_mat_dir, 1)) - vmodif([1:size(vmodif, 1)]' * ones(1, nb_mat_dir), :);

% Thresholding of displacement differences dist_vec_dir
% Distance follows Rayleigh distribution: sqrt(b1^2 + b2^2), with b1, b2 ~ N(0, 4 * sigma^2)
% Numerical thresholds:
% Probability > 90% => seuil = 2.1 * (2 * sigma)
% Probability > 95% => seuil = 2.4 * (2 * sigma)
% Probability > 99% => seuil = 3.0 * (2 * sigma)
seuil = 6 * sigma;
detec_dist_vec_dir = abs(mat_dist_vec_dir) < seuil;

% Select candidate displacement directions excluding self-comparisons
selec_dist_vec_dir = detec_dist_vec_dir - diag(diag(detec_dist_vec_dir));
ind_selec = find(sum(selec_dist_vec_dir)).' ;

vec_dir_selec = vec_dir(ind_selec);
norm_vec_dir_selec = abs(vec_dir_selec);

% Compute index pairs: d_i + 1i * g_j (indexing convention)
% n = d_i + nb_part_d * (g_j - 1)
ind_g = ceil(ind_selec / nb_part_d);
ind_d = ind_selec - (ind_g - 1) * nb_part_d;
ind_couple = ind_d + 1i * ind_g;

% Step 2: Search for similar directions (largest cluster)
% Sort candidate displacements by norm (ascending)
% This is a two-dimensional problem. We consider
% a sub-optimal solution by searching for
% similar directions, sorted by their magnitude.
% We choose the largest set of similar directions.
[sort_abs_vec_dir, ind_dir] = sort(norm_vec_dir_selec, 'ascend');
vec_dir_sort_abs_vec_dir = vec_dir_selec(ind_dir).' ;

% Compute successive differences between sorted vectors
diff_vec_dir_sort_abs_vec_dir = [0, vec_dir_sort_abs_vec_dir(1:end-1) - vec_dir_sort_abs_vec_dir(2:end)];

% Thresholding successive diff_abs_vec_dir
% Its norm is a random variable of Rice type: sqrt((x0 + b1)^2 + (y0 + b2)^2) with b1 and b2 ~ N(0, 2*sigma^2)
% abs((x0, y0)) >> abs((b1, b2))

% We approximate it by a random variable |d0| + b with : N(0, 2*sigma^2)
% We threshold the subtraction of two norms
% Difference is a random variable : b' with: N(0, 4 * sigma^2)

% This is a Gaussian distribution, threshold calculation:
% Probability > 95% => seuil = 2.1 * (2 * sigma)
% Probability > 99% => seuil = 2.5 * (2 * sigma)
% Probability > 99.7% => seuil = 3.0 * (2 * sigma)
% seuil = 4.2*sigma ;
seuil = 10 * sigma;
% seuil = 100*sigma ;
inf_seuil = abs(diff_vec_dir_sort_abs_vec_dir) < seuil;
inf_seuil(end) = 0;

% Identify clusters of similar displacements
diff_inf_seuil = [1, inf_seuil(2:end) - inf_seuil(1:end-1)];
ind_debut_inf_seuil = find(diff_inf_seuil == 1);
ind_fin_inf_seuil = find(diff_inf_seuil == -1);
long_inf_seuil = ind_fin_inf_seuil - ind_debut_inf_seuil;

[tmp, ind_max] = max(long_inf_seuil);
ind_debut = ind_debut_inf_seuil(ind_max);
ind_fin = ind_fin_inf_seuil(ind_max) - 1;

% Extract the largest set of similar displacements
vec_dir_simil = vec_dir_selec(ind_dir(ind_debut:ind_fin));
dir_simil = ind_couple(ind_dir(ind_debut:ind_fin));

% Number of estimated pairs
nb_est = size(vec_dir_simil, 1);

if (nb_est >= dir_min)
    % Compute mean displacement
    vec_dir_est = mean(vec_dir_simil);
    decal_i = real(vec_dir_est);
    decal_j = imag(vec_dir_est);

    % Build list of matched particle indices
    indice_couple = [real(dir_simil.'); imag(dir_simil.')];
    indice_couple = indice_couple(:).';
else
    % Not enough matches, return zeros
    decal_i = 0;
    decal_j = 0;
    nb_est = 0;
    indice_couple = [0, 0];
end

end %function
