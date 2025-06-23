% =========================================================================
% Function: separe_liste
%
% Description:
%   This function separates a list of coordinates into two groups based on
%   their position relative to a vertical separation line placed at the center of the image.
%   The function returns the right-side and left-side coordinate lists with their indices.
%
% Inputs:
%   liste_i   - List of i (row) coordinates
%   liste_j   - List of j (column) coordinates
%
% Outputs:
%   liste_i_d - i coordinates on the right side
%   liste_j_d - j coordinates on the right side
%   ind_d     - Indices of the right-side coordinates
%   liste_i_g - i coordinates on the left side
%   liste_j_g - j coordinates on the left side
%   ind_g     - Indices of the left-side coordinates
%
% Notes:
%   The separation threshold 'sep' is read from the MATLAB base workspace.
%   This function can be adapted depending on the specific experimental setup.
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: December 2017
% =========================================================================

function [liste_i_d, liste_j_d, ind_d, liste_i_g, liste_j_g, ind_g] = separe_liste(liste_i, liste_j)

% Separation line: can be defined directly or read from base workspace
% Example: sep = 256;
sep = evalin('base', 'sep'); % Reads 'sep' from the base workspace

% Identify points on the right side of the separation line
pos_d = liste_j > sep;
liste_i_d = liste_i(pos_d);
liste_j_d = liste_j(pos_d);
ind_d = find(pos_d);

% Identify points on the left side of the separation line
pos_g = liste_j < sep;
liste_i_g = liste_i(pos_g);
liste_j_g = liste_j(pos_g);
ind_g = find(pos_g);

end % function
