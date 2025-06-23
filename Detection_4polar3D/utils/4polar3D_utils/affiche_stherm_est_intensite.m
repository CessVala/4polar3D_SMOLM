% =========================================================================
% Function: affiche_stherm_est_intensite
%
% Description:
%   This function visualizes trajectory connections by plotting colored
%   segments between trajectory pairs (left and right) on an image.
%   The color indicates the coupling type based on trajectory properties.
%
% Inputs:
%   im_t            - Input image frame
%   num_t           - Frame number
%   liste_traj_dg   - Trajectory list:
%                     [ param_traj_1_d, param_traj_1_g, param_traj_2_d, param_traj_2_g, ... ]
%
% Outputs:
%   r, v, b         - Red, green, blue image channels (RGB)
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: December 2017
% =========================================================================

function [r,v,b] = affiche_stherm_est_intensite(im_t, num_t, liste_traj_dg)

global param_t_num param_t_i param_t_j param_t_sig2 param_t_etat

nb_part = size(liste_traj_dg, 2) / 2; % Number of trajectory pairs

[M, N] = size(im_t);
dim_max = N * M;

%% Rescale image to 0:255
maxi = max(im_t(:));
mini = min(im_t(:));

r = 255 * (im_t - mini) / (maxi - mini);
v = r;
b = r;
rr = r;
vr = v;
br = b;

%% Display frame number in the image
liste_pix = mat2mat3x5(M + 1, N, M, num_t);
r(liste_pix) = 255;
v(liste_pix) = 255;
b(liste_pix) = 255;

%% Colors for different coupling types
couleur_dd_dg = [0, 255, 0];   % Both trajectories valid
couleur_nd_dg = [255, 0, 0];   % First trajectory invalid
couleur_dd_ng = [0, 0, 255];   % Second trajectory invalid
coef_nb_pix = 1;

%% Draw segments between selected trajectory pairs
for c = 1:2:2 * nb_part

    i1 = liste_traj_dg(param_t_i, c);
    j1 = liste_traj_dg(param_t_j, c);
    i2 = liste_traj_dg(param_t_i, c + 1);
    j2 = liste_traj_dg(param_t_j, c + 1);

    nb_pix = abs(floor((j2 - j1) / 2));
    [vi, vj] = get_coord_line(i1, j1, i2, j2, nb_pix);

    % Select color based on trajectory status
    if (liste_traj_dg(param_t_num, c) < 0)
        color = couleur_nd_dg;
        coef_nb_pix = 2;
    else
        if (liste_traj_dg(param_t_num, c + 1) < 0)
            color = couleur_dd_ng;
            coef_nb_pix = 2;
        else
            color = couleur_dd_dg;
            coef_nb_pix = 1;
        end
    end

    % Draw the segment
    for p = 1:coef_nb_pix:nb_pix
        mi = vi(p) + vj(p) * M;
        if ((mi > 0) && (mi < dim_max))
            r(mi) = color(1);
            v(mi) = color(2);
            b(mi) = color(3);
        end
    end
end

%% Concatenate original and modified images (display fix)
r = [rr, 200 * ones(M, 1), r];
v = [vr, 200 * ones(M, 1), v];
b = [br, 200 * ones(M, 1), b];

%% Reset random generator
rand('seed', cputime); %#ok

end % function

%% Custom ceiling function
function f = ceil_(x)
sng = sign(x);
f = sng .* ceil(sng .* x);
f = f .* (abs(x) > 1e-4);
end % function

%% Return n coordinates along the line between P1 and P2
function [vi, vj] = get_coord_line(i1, j1, i2, j2, n)

di = i2 - i1;
dj = j2 - j1;

vi = zeros(1, n);
vj = zeros(1, n);
pas = 1 / (n - 1);
coef = 0;

for p = 1:n
    vi(p) = round(i1 + coef * di);
    vj(p) = round(j1 + coef * dj);
    coef = coef + pas;
end
end % function

%% Display frame number in the image using a 3x5 font
function liste_pix = mat2mat3x5(pos, N, M, num)

n0 = [...
    1 1 1;...
    1 0 1;...
    1 0 1;...
    1 0 1;...
    1 1 1] ;%#ok
n1 = [...
    0 1 0;...
    0 1 0;...
    0 1 0;...
    0 1 0;...
    0 1 0] ;%#ok
n2= [...
    1 1 1;...
    0 0 1 ;...
    1 1 1;...
    1 0 0;...
    1 1 1] ;%#ok
n3= [...
    1 1 1;...
    0 0 1 ;...
    0 1 1;...
    0 0 1;...
    1 1 1] ;%#ok
n4= [...
    1 0 0;...
    1 0 1 ;...
    1 1 1;...
    0 0 1;...
    0 0 1] ;%#ok
n5 = [...
    1 1 1;...
    1 0 0;...
    1 1 1;...
    0 0 1;...
    1 1 1] ;%#ok
n6= [...
    1 1 1;...
    1 0 0;...
    1 1 1;...
    1 0 1;...
    1 1 1] ;%#ok
n7 = [...
    1 1 1;...
    0 0 1;...
    0 1 0;...
    0 1 0;...
    0 1 0] ;%#ok
n8 = [...
    1 1 1;...
    1 0 1;...
    1 1 1;...
    1 0 1;...
    1 1 1] ;%#ok
n9 = [...
    1 1 1;...
    1 0 1;...
    1 1 1;...
    0 0 1;...
    1 1 1] ;%#ok

strnum = num2str(num);

cmd = sprintf('[ii, jj] = find(n%s);', strnum(1));
eval(cmd);
liste_pix = pos + (ii') + jj' * M;

for c = 2:size(strnum, 2)
    cmd = sprintf('[ii, jj] = find(n%s);', strnum(c));
    eval(cmd);
    liste_pix = [liste_pix, pos + ii' + (jj' + (c - 1) * 4) * M]; %#ok
end
end % function
