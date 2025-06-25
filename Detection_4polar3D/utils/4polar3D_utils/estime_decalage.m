% =========================================================================
% Function Name: estime_decalage.m
%
% Description:
%   This function estimates the displacement (shift) between two sub-images
%   in 4polar3D data.
%
%   The function:
%     - Assumes two provided lists of detected particles (right and left).
%     - Computes the optimal displacement between the two sub-images.
%
% Inputs:
%   - vec_i: Row vector of particle positions along the i-axis.
%   - vec_j: Row vector of particle positions along the j-axis.
%   - sigma: Standard deviation of position errors on i and j axes.
%   - dir_min: Minimum number of matching particle pairs required per translation.
%
% Outputs:
%   - decal_i: Estimated displacement along the i-axis.
%   - decal_j: Estimated displacement along the j-axis.
%   - nb_est: Number of particle pairs considered in the final displacement estimation.
%   - indice_couple: Indices of the matched particle pairs used for estimation
%                    (vector of size 2 * nb_est) formatted as [d1, g1, d2, g2, ...].
%
% Notes:
%   - This version assumes a left/right sub-image geometry as provided by separe_liste.m.
%   - If the geometry is different, separe_liste.m must be modified accordingly.
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function [decal_i, decal_j, nb_est, indice_couple] = estime_decalage(vec_i, vec_j, sigma, dir_min)

% Separate the input list into right and left lists
[vec_d_i, vec_d_j, ind_d, vec_g_i, vec_g_j, ind_g] = separe_liste(vec_i, vec_j);

% Estimate displacement using right and left lists
[decal_i, decal_j, nb_est, indice_couple_dg] = estime_decalage_dg(vec_d_i, vec_d_j, vec_g_i, vec_g_j, sigma, dir_min);

if (nb_est >= dir_min)
    % Extract the matching indices for the right and left lists
    ind_d_simil = indice_couple_dg(1:2:end);
    ind_g_simil = indice_couple_dg(2:2:end);
    indice_couple = [ind_d(ind_d_simil); ind_g(ind_g_simil)];
    indice_couple = indice_couple(:).';
else
    indice_couple = [0, 0];
end

end %function