% =========================================================================
% Function Name: STHERM_estime_intensite_ver3.m
%
% Description:
%   This function processes 4polar3D data to automatically couple trajectories
%   detected in the two sub-image zones (H-V) and estimate their intensities.
%   This version works directly from a provided list of detected particles.
%
% Inputs:
%   - stkrange: Number of images in the stack or image index range.
%               If stkrange = N, images 1 to N are processed.
%               If stkrange = [Deb, Fin], images from Deb to Fin are processed.
%   - rep_in: Input folder containing the images.
%   - data_file: Input file name.
%   - decal_i: Estimated displacement along the i-axis between sub-images.
%   - decal_j: Estimated displacement along the j-axis between sub-images.
%   - list_trajbase: Complete list of detected particles (all frames).
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

function [liste_traj_output, liste_traj_all] = STHERM_estime_intensite_ver3(stkrange, rep_in, data_file, decal_i, decal_j, list_trajbase)

% Define global parameters
STHERM_param;
DISPLAY = 0; % Disable display mode by default

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

liste_traj_output = []; % List of coupled trajectories
liste_traj_all = [];    % List of all detected trajectories

% Loop over each frame
for i = Deb:Fin

    % Select trajectories detected in frame i
    liste_traj = list_trajbase(:, list_trajbase(8, :) == i);

    im_i = 0; % Dummy variable for display compatibility

    % Estimate intensity and pair particles between sub-images
    [liste_traj_dg, nb_couple_dec, nb_couple_est] = estime_intensite(liste_traj, im_i, decal_i, decal_j, im_i) ;

    if ((nb_couple_dec > 0) || (nb_couple_est > 0))
        fprintf(display_out, 'Number of coupled trajectories (detected/estimated) %d/%d (#%d)\n', nb_couple_dec, nb_couple_est, i);

        % Store coupled trajectories
        liste_traj_output = [liste_traj_output, liste_traj_dg] ;

        % Optional display if enabled
        if (DISPLAY)
            [r, v, b] = affiche_stherm_est_intensite(im_i, i, liste_traj_dg);
            eval(cmd_output_i); % Generate display output filename
            ppm8write(outfile, r, v, b);
        end
    end

    % Store all detected trajectories
    liste_traj_all = [liste_traj_all, liste_traj];
end % for i = Deb:Fin

% Assign trajectory numbers to coupled trajectories
nb_traj_couple = size(liste_traj_output, 2) / 2;
liste_traj_output(param_t_num, :) = kron(1:nb_traj_couple, [1 1]);

fprintf(display_out, 'Total number of coupled trajectories: %d\n', nb_traj_couple);

end %function