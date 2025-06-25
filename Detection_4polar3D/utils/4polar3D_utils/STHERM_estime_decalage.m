% =========================================================================
% Function Name: STHERM_estime_decalage.m
%
% Description:
%   This function processes 4polar3D data to automatically estimate
%   the displacement (shift) between two sub-image zones (H-V).
%
% Inputs:
%   - stkrange: Number of images in the stack or image index range.
%               If stkrange = N, images 1 to N are processed.
%               If stkrange = [Deb, Fin], images from Deb to Fin are processed.
%   - rep_in: Input folder containing the images.
%   - data_file: Input file name.
%
% Notes:
%   - To modify image loading, change the function read_stk4specsingle().
%   - To modify the H-V separation, change the function separe_liste.m.
%   - The term "Trajectories" is used here because the original code was created for single particle tracking,
%     even though in this context they are not real trajectories but particles estimated and coupled across frames.
%   - "Particles" refer to particles detected in the current frame.
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function [decal_i, decal_j] = STHERM_estime_decalage(stkrange, rep_in, data_file)

% Define global parameters
STHERM_param;

% If display mode is active
if (DISPLAY)
    cmd_output_i = 'outfile = sprintf(''%s%sdisplay1%sestime_decalage_%5.5d.ppm'', rep_in,char_slash,char_slash, i) ;';
    coefstd = 50;
end

% Determine image index range
if (max(size(stkrange)) == 1)
    Deb = 1;
    Fin = stkrange;
else
    Deb = stkrange(1);
    Fin = stkrange(2);
end

% Prepare image loading
name_image = [rep_in, char_slash, data_file];
initload = Fin;

% Load first image to initialize (if needed)
im_1 = read_stk4specsingle(name_image, Deb, initload);
initload = Deb;

tab_decal_ij = []; % Table of displacement estimates
nb_est_all = 0;    % Total number of matched pairs used for all estimates
nb_dir_est = 0;    % Number of displacement directions estimated

last_num_image = Fin + 1;

% Loop over each image in the stack
for i = Deb:Fin

    % Load current image
    im_i = read_stk4specsingle(name_image, i, initload);

    % Detect and estimate particles in the image
    liste_part_detect_all = detect_et_estime_part_1vue(im_i, wn, r0, pfa).';

    % Filter particles with valid detection flag
    indice_ok = (liste_part_detect_all(param_p_ok, :) == 1);

    % Limit number of selected particles to 300 if necessary
    if size(find(indice_ok), 2) > 300
        indice_oktemp = indice_ok;
        posok = find(indice_ok);
        selposok = randperm(size(posok, 2), 300);
        posok1 = posok(selposok);
        indice_ok = zeros(size(indice_ok));
        indice_ok(posok1) = 1;
        indice_ok = indice_ok == 1;
    end

    liste_part_detect = liste_part_detect_all(:, indice_ok);
    nb_part_find = size(liste_part_detect, 2);

    % Build the trajectory list
    NB_traj = size(liste_part_detect, 2);
    liste_traj = zeros(nb_param_t, NB_traj);

    liste_traj(param_t_num, :) = (1:NB_traj);
    liste_traj(param_t_i, :) = liste_part_detect(param_p_i, :);
    liste_traj(param_t_j, :) = liste_part_detect(param_p_j, :);
    liste_traj(param_t_r0, :) = liste_part_detect(param_p_r0, :);

    % Compute trajectory variance using BCR (Cramer-Rao Bound)
    % Using known r0
    %	liste_traj(param_t_sig2,:) = BCR_1G_ij_rconnnu(wn, liste_part_detect(param_p_sig2,:), r0, liste_part_detect(param_p_alpha,:)) ;
    % Test version with r0 estimated
    % we take the BCR with known r0 (whereas it is estimated => underestimation of sig2_ij)
    liste_traj(param_t_sig2, :) = BCR_1G_ij_rconnu_liste_r(wn, liste_part_detect(param_p_sig2, :), liste_part_detect(param_p_r0, :), liste_part_detect(param_p_alpha, :));

    % Estimate displacement using selected particles
    seuil_mini = 1.5; % Minimum threshold
    vec_i = liste_part_detect(param_p_i, :);
    vec_j = liste_part_detect(param_p_j, :);

    seuilstd = mean(sqrt(liste_traj(param_t_sig2, :))); % Average standard deviation

    if (seuilstd < seuil_mini)
        [decal_i, decal_j, nb_est, indice_selection] = estime_decalage(vec_i, vec_j, seuilstd, 3);

        % Store estimated displacement if valid
        if (nb_est > 0)
            tab_decal_ij = [tab_decal_ij; decal_i, decal_j];
            nb_est_all = nb_est_all + nb_est;
            nb_dir_est = nb_dir_est + 1;

            fprintf(display_out, '\nAverage direction: (%.2f, %.2f) with %d matches (#%d/#%d)\n\n', decal_i, decal_j, nb_est, nb_dir_est, i);

            % Save display image if display mode is active
            if (DISPLAY)
                [r, v, b] = affiche_stherm_debug(im_i, i, liste_traj, coefstd, indice_selection);
                eval(cmd_output_i); % Generate output file name
                ppm8write(outfile, r, v, b);
            end
        end
    end
end %for i=Deb:Fin

fprintf(display_out, 'Number of estimated directions: %d (total matches: %d)\n', nb_dir_est, nb_est_all);

% Final displacement estimation (using median)
if (nb_dir_est >= 3)
    median_all_decal = median(tab_decal_ij);
    decal_i = median_all_decal(1);
    decal_j = median_all_decal(2);

    % Display displacement statistics
    fprintf(display_out, 'Final displacement estimate: %.2f , %.2f\n', median_all_decal);
    fprintf(display_out, 'Minimum displacement: %.2f , %.2f\n', min(tab_decal_ij));
    fprintf(display_out, 'Maximum displacement: %.2f , %.2f\n', max(tab_decal_ij));

else
    decal_i = 0;
    decal_j = 0;
end

end %function
