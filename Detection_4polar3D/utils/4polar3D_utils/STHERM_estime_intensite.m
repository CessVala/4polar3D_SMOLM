% =========================================================================
% Function Name: STHERM_estime_intensite.m
%
% Description:
%   This function processes 4polar3D data to automatically couple trajectories
%   detected in the two sub-image zones (H-V) and estimate their intensities.
%
% Inputs:
%   - stkrange: Number of images in the stack or image index range.
%               If stkrange = N, images 1 to N are processed.
%               If stkrange = [Deb, Fin], images from Deb to Fin are processed.
%   - rep_in: Input folder containing the images.
%   - data_file: Input file name.
%   - decal_i: Estimated displacement along the i-axis between sub-images.
%   - decal_j: Estimated displacement along the j-axis between sub-images.
%
% Outputs:
%   - liste_traj_output: List of coupled trajectories (right and left) with intensity estimation.
%                        Each column represents a particle, paired as [right, left, right, left, ...].
%	  	from STHERM_param.m
% 		param_t_num
% 		param_t_i     : i-position
% 		param_t_j     : j-position
% 		param_t_r0    : PSF width
% 		param_t_alpha : intensite
% 		param_t_sig2  : ij variance
% 		param_t_tps   : time index (used for drift correction)
%
%   - liste_traj_all: Complete list of detected trajectories (including uncoupled) for all frames.
%
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

function [liste_traj_output, liste_traj_all] = STHERM_estime_intensite(stkrange, rep_in, data_file, decal_i, decal_j)

% Define global parameters
STHERM_param;

% Display settings if enabled
if (DISPLAY)
    cmd_output_i = 'outfile = sprintf(''%s%sdisplay1%sestime_intensite_%5.5d.ppm'', rep_in,char_slash,char_slash, i) ;';
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
name_image = [rep_in, char_slash, data_file] ;
initload = Fin ;

% Load first image to initialize
im_1 = read_stk4specsingle(name_image, Deb, initload) ;
initload = Deb ;

tab_decal_ij = [] ;
nb_est_all = 0 ;
nb_dir_est = 0 ;

last_num_image = Fin + 1;
liste_traj_output = [] ;
liste_traj_all = [] ;

% Loop over each frame
for i = Deb:Fin

    % Load current image
    im_i = read_stk4specsingle(name_image, i, initload);

    % Detect and estimate particles in the image
    liste_part_detect_all = detect_et_estime_part_1vue(im_i, wn, r0, pfa).';

    % Filter particles with valid detection flag
    indice_ok = (liste_part_detect_all(param_p_ok, :) == 1);
    liste_part_detect = liste_part_detect_all(:, indice_ok);

    nb_part_find = size(liste_part_detect, 2);

    % Build trajectory list
    NB_traj = size(liste_part_detect, 2);
    liste_traj = zeros(nb_param_t, NB_traj);

    liste_traj(param_t_num, :) = (1:NB_traj);
    liste_traj(param_t_i, :) = liste_part_detect(param_p_i, :);
    liste_traj(param_t_j, :) = liste_part_detect(param_p_j, :);
    liste_traj(param_t_r0, :) = liste_part_detect(param_p_r0, :);
    liste_traj(param_t_alpha, :) = liste_part_detect(param_p_alpha, :);
    liste_traj(param_t_sig2_b, :) = liste_part_detect(param_p_sig2, :);
    liste_traj(param_t_tps, :) = i;

    % Compute trajectory variance using BCR
    % version with known r0
    % liste_traj(param_t_sig2,:) = BCR_1G_ij_rconnnu(wn, liste_part_detect(param_p_sig2,:), r0, liste_part_detect(param_p_alpha,:)) ;

    % test version with estimated r0
    % we use the BCR with known r0 (although it is estimated => underestimation of sig2_ij)
    liste_traj(param_t_sig2, :) = BCR_1G_ij_rconnu_liste_r(wn, liste_part_detect(param_p_sig2, :), liste_part_detect(param_p_r0, :), liste_part_detect(param_p_alpha, :));

    % Search for particle pairs and estimate intensities
    [liste_traj_dg, nb_couple_dec, nb_couple_est] = estime_intensite(liste_traj, im_i, decal_i, decal_j, im_i);

    if ((nb_couple_dec > 0) || (nb_couple_est > 0))
        fprintf(display_out, 'Number of coupled trajectories (detected/estimated) %d/%d (#%d)\n', nb_couple_dec, nb_couple_est, i);

        % Store coupled trajectories
        liste_traj_output = [liste_traj_output, liste_traj_dg];

        % Optional display if enabled
        if (DISPLAY)
            [r, v, b] = affiche_stherm_est_intensite(im_i, i, liste_traj_dg);
            eval(cmd_output_i); % Generate display output filename
            ppm8write(outfile, r, v, b);
        end
    end

    % Store all trajectories (including uncoupled)
    liste_traj_all = [liste_traj_all, liste_traj];
end % for i = Deb:Fin

% Assign trajectory numbers to coupled trajectories
nb_traj_couple = size(liste_traj_output, 2) / 2;
liste_traj_output(param_t_num, :) = kron(1:nb_traj_couple, [1 1]);

fprintf(display_out, 'Total number of coupled trajectories: %d\n', nb_traj_couple);

end %function
