% =========================================================================
% Function: affiche_stherm_debug
%
% Description:
%   This function visualizes trajectory selection and 4polar3D analysis
%   by generating a colored RGB image overlay showing the selected
%   trajectories and optionally their IDs.
%
% Inputs:
%   im_t               - Input image frame
%   num_t              - Frame number
%   liste_traj         - Trajectory list
%   coef_std           - Standard deviation coefficient
%   indice_selection   - Selected trajectory indices
%
% Outputs:
%   r, v, b            - Red, green, blue image channels (RGB)
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: December 2017
% =========================================================================

function [r,v,b] = affiche_stherm_debug(im_t, num_t, liste_traj, coef_std, indice_selection)

global param_t_num param_t_i param_t_j param_t_sig2 param_t_etat

nb_part = size(indice_selection, 2); % Number of selected trajectories

num_traj = 0; % 0/1: display trajectory number

[M,N] = size(im_t);
dim_max = N * M;

%% Define circle shape
pas_th = pi / 36;
theta = pas_th : pas_th : 2.0 * pi;
cosinus = cos(theta);
sinus = sin(theta);

%% Rescale image to 0:255
maxi = max(im_t(:));
mini = min(im_t(:));

r = 255 * (im_t - mini) / (maxi - mini);
v = r;
b = r;
rr = r;
vr = v;
br = b;

%% Display frame number
liste_pix = mat2mat3x5(M+1, N, M, num_t);
r(liste_pix) = 255;
v(liste_pix) = 255;
b(liste_pix) = 255;

%% Loop over all selected trajectories
for traj = 1:nb_part

    ind_traj = indice_selection(traj);

    %% Trajectory parameters
    ii = round(liste_traj(param_t_i, ind_traj));
    jj = round(liste_traj(param_t_j, ind_traj));

    stdij = coef_std * sqrt(liste_traj(param_t_sig2, ind_traj));
    numero = liste_traj(param_t_num, ind_traj);
    etat = liste_traj(param_t_etat, ind_traj);

    %% Trajectory color
    rand('seed', numero); %#ok
    cr = 255 * (0.6 * rand + 0.4);
    cv = 255 * (0.6 * rand + 0.4);
    cb = 255 * (0.6 * rand + 0.4);

    %% Blinking state
    if (etat > 0)
        coef = 1; % Active
    else
        coef = 0.5; % Blinking
    end

    %% Confinement circle
    pas_theta = 2;
    ci = ii + ceil_(stdij * cosinus(1:pas_theta:end));
    sj = jj + ceil_(stdij * sinus(1:pas_theta:end));

    %% Limit to image boundaries
    out_vert = ci < M;
    ci = ci(out_vert);
    sj = sj(out_vert);
    mi = ci + sj * M;

    %% Image wrap-around
    mi = 1 + mod(mi, dim_max);

    %% Display the pattern
    r(mi) = coef * cr;
    v(mi) = coef * cv;
    b(mi) = coef * cb;

    %% Optionally display trajectory number
    if (num_traj)
        liste_pix = mat2mat3x5(max(mi) + 1, N, M, numero);
        liste_pix = 1 + liste_pix(liste_pix < dim_max);
        r(liste_pix) = coef * cr;
        v(liste_pix) = coef * cv;
        b(liste_pix) = coef * cb;
    end

end

%% Draw lines between coupled trajectories
for c = 1:2:nb_part

    vec_i = liste_traj(param_t_i, :);
    vec_j = liste_traj(param_t_j, :);

    ind_traj1 = indice_selection(c);
    ind_traj2 = indice_selection(c + 1);

    i1 = liste_traj(param_t_i, ind_traj1);
    j1 = liste_traj(param_t_j, ind_traj1);
    i2 = liste_traj(param_t_i, ind_traj2);
    j2 = liste_traj(param_t_j, ind_traj2);

    nb_pix = abs(floor((j2 - j1) / 2));
    [vi, vj] = get_coord_line(i1, j1, i2, j2, nb_pix);

    for p = 1:nb_pix
        mi = vi(p) + vj(p) * M;
        r(mi) = 192;
        v(mi) = 192;
        b(mi) = 192;
    end
end

%% Concatenate original and modified images (fix display bug)
r = [rr, 200 * ones(M, 1), r];
v = [vr, 200 * ones(M, 1), v];
b = [br, 200 * ones(M, 1), b];

%% Reset random generator
rand('seed', cputime); %#ok

end

%% Custom ceiling function
function f = ceil_(x)
sng = sign(x);
f = sng .* ceil(sng .* x);
f = f .* (abs(x) > 1e-4);
end

%% Returns n coordinates along the line between P1 and P2
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
end

%% Display frame number in image
function liste_pix = mat2mat3x5(pos, N, M, num)

% Digit templates (3x5 pixel font)

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
end
