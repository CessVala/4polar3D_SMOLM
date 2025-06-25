% =========================================================================
% Function Name: STHERMver3.m
%
% Description:
%   This function processes 4polar3D data to automatically estimate the displacement
%   between two sub-image zones (H-V) and couple the corresponding particle trajectories.
%   This version works directly from a provided list of detected particles.
%
% Inputs:
%   - stkrange_dec: Number of images or index range for displacement estimation.
%   - stkrange_int: Number of images or index range for intensity estimation.
%   - rep_in: Input folder containing the images.
%   - data_file: Input file name.
%   - rep_out: Output folder for saving results.
%   - out_file: Output file name for saving results.
%   - tab_trajsel: Provided list of detected particles (all frames).
%   - decal_i: Estimated displacement along the i-axis between sub-images (optional).
%   - decal_j: Estimated displacement along the j-axis between sub-images (optional).
%
% Outputs:
%   - tab_traj: List of coupled trajectories (right and left) with intensity estimation.
%   - tab_traj_all: Complete list of detected trajectories (including uncoupled) for all frames.
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

function [tab_traj,tab_traj_all]=STHERMver3(stkrange_dec, stkrange_int, rep_in, data_file, rep_out, out_file, tab_trajsel,decal_i, decal_j)

% Define global parameters
STHERM_param ;

% If displacement was not provided, estimate it
if (nargin <= 7)
    fprintf(display_out, 'Estimating displacement between H/V zones...\n');
    [decal_i, decal_j] = STHERM_estime_decalagever3(stkrange_dec, rep_in, data_file,tab_trajsel) ;
end

% Stop if displacement estimation failed
if ((decal_i == 0) && (decal_j == 0))
    fprintf(display_out, 'Displacement estimation failed, please perform manual alignment!\n');
    return ;
end

% Estimate and couple trajectories
fprintf(display_out, 'Coupling trajectories ...\n');
[tab_traj,tab_traj_all] = STHERM_estime_intensite_ver3(stkrange_int, rep_in, data_file, decal_i, decal_j,tab_trajsel) ;

% Save relevant parameters and results
param_wn = wn ;
param_r0 = r0 ;
param_pfa = pfa ;

name_savefile = sprintf('%s%s%s', rep_out, char_slash, out_file) ;
name_savefile

save(name_savefile, 'decal_i', 'decal_j', 'param_wn', 'param_r0', 'param_pfa', 'tab_traj', 'stkrange_dec', 'stkrange_int', 'data_file','tab_traj_all') ;

end %function
