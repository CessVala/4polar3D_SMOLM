% =========================================================================
% Function Name: STHERM_estime_decalagever3.m
%
% Description:
%   This function processes 4polar3D data to automatically estimate
%   the displacement (shift) between two sub-image zones (H-V).
%   This version works directly from a list of detected particles.
%
% Inputs:
%   - stkrange: Number of images in the stack or image index range.
%               If stkrange = N, images 1 to N are processed.
%               If stkrange = [Deb, Fin], images from Deb to Fin are processed.
%   - rep_in: Input folder containing the images.
%   - data_file: Input file name.
%   - list_trajbase: Complete list of detected particles (all frames).
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

function [decal_i, decal_j] = STHERM_estime_decalagever3(stkrange, rep_in, data_file, list_trajbase)

% Define global parameters
STHERM_param;

DISPLAY = 0; % Disable display mode by default

% Display settings if enabled
if (DISPLAY)
    cmd_output_i = 'outfile = sprintf(''%s%sdisplay1%sestime_decalage_%5.5d.ppm'', rep_in,char_slash,char_slash, i) ;';
    coefstd = 50 ;
end

% Determine image index range
if (max(size(stkrange)) == 1)
    Deb = 1;
    Fin = stkrange;
else
    Deb = stkrange(1);
    Fin = stkrange(2);
end

% No need to load stacks. This version uses the provided particle list only
name_image = [rep_in, char_slash, data_file];

tab_decal_ij = []; % Table of displacement estimates
nb_est_all = 0;    % Total number of matched pairs used for all estimates
nb_dir_est = 0;    % Number of displacement directions estimated

% Loop over each frame
for i = Deb:Fin

    % Select particles detected in frame i
    liste_part_detect = list_trajbase(:, list_trajbase(8, :) == i);

    % Search for pairs for shift estimation
    seuil_mini = 1.5;

    vec_i = liste_part_detect(param_p_i, :);
    vec_j = liste_part_detect(param_p_j, :);
    seuilstd = mean(sqrt(liste_part_detect(param_p_sig2, :)));

    if (seuilstd < seuil_mini)
        % Estimate displacement based on detected particles
        [decal_i, decal_j, nb_est, indice_selection] = estime_decalage(vec_i, vec_j, seuilstd, 3);

        % Store valid estimates
        if (nb_est > 0)
            tab_decal_ij = [tab_decal_ij; decal_i, decal_j] ;
            nb_est_all = nb_est_all + nb_est ;
            nb_dir_est = nb_dir_est + 1 ;

            fprintf(display_out, '\nAverage direction: (%.2f, %.2f) with %d matches (#%d/#%d)\n\n', decal_i, decal_j, nb_est, nb_dir_est, i);

            im_i = 0; % Dummy variable for display compatibility

            % Optional debug display if enabled
            if (DISPLAY)
                [r, v, b] = affiche_stherm_debug(im_i, i, list_trajbase, coefstd, indice_selection);
                eval(cmd_output_i); % Generate display output filename
                outfile
                ppm8write(outfile, r, v, b);
            end
        end
    end
end % for i = Deb:Fin

fprintf(display_out, 'Number of estimated average directions: %d (total matches: %d)\n', nb_dir_est, nb_est_all);

% Final displacement estimation using median if enough valid directions
if (nb_dir_est >= 3)
    median_all_decal = median(tab_decal_ij);
    decal_i = median_all_decal(1);
    decal_j = median_all_decal(2);

    % Display displacement statistics
    fprintf(display_out, 'Final displacement estimate: %.2f , %.2f\n', median_all_decal);
    fprintf(display_out, 'Minimum displacement: %.2f , %.2f\n', min(tab_decal_ij));
    fprintf(display_out, 'Maximum displacement: %.2f , %.2f\n', max(tab_decal_ij));

else
    % If insufficient valid directions, return zeros
    decal_i = 0;
    decal_j = 0;
end

end %function
