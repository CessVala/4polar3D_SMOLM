%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: K_from_beads
%
% Description:
%   Extract the K matrix, normalized intensities, and M matrix from
%   fluorescent bead measurements.
%
% Inputs:
%   - p_angles: p angles in radians
%   - start_frame: first valid frame index to use
%   - n_frame: number of frames
%   - varargin: (optional) direct path to the beads .mat file
%
% Outputs:
%   - K_all: computed K matrix (normalized)
%   - I_all: extracted intensities
%   - M_all: M matrix
%   - p_ind: angle index
%
% Authors: Charitra S. Senthil Kumar, Miguel Sison, Cesar VALADES-CRUZ
% Date: March 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K_all, I_all, M_all, p_ind] = K_from_beads(p_angles, start_frame, n_frame, varargin)

    %% Load bead file (via dialog or direct input)
    if nargin == 3
        [fname, fpath] = uigetfile(); % Prompt user to select file
        beads_file = fullfile(fpath,fname);
    elseif nargin == 4
        beads_file = varargin{1}; % File provided as argument
    else
        errordlg('Wrong Inputs!!!') % Input check
        return
    end

    %% Load data from selected file
    load(beads_file, 'tab_traj_4x4', ...
        'param_t_r0', 'param_t_alpha', 'param_t_tps');

    % Extract data for each polarization channel
    tab_traj_0   = tab_traj_4x4(:, 1:4:end); % 0°
    tab_traj_45  = tab_traj_4x4(:, 2:4:end); % 45°
    tab_traj_90  = tab_traj_4x4(:, 3:4:end); % 90°
    tab_traj_135 = tab_traj_4x4(:, 4:4:end); % 135°

    % Extract frame numbers and adjust based on start_frame
    T = [tab_traj_0(param_t_tps, :); tab_traj_45(param_t_tps, :); ...
         tab_traj_90(param_t_tps, :); tab_traj_135(param_t_tps, :)];
    T = T(1,:) - start_frame;

    % Exclude invalid frame indices
    ex_t = ( T < 0 ) | ( T >= numel(p_angles). * n_frame );

    T = T(~ex_t); % Keep valid indices only

    p_ind = ( (T - mod(T, n_frame)) ./ n_frame ) + 1;

    %% Build model matrix M_all
    M_all = zeros(4, numel(p_ind)); % Preallocate

    for k = 1:numel(p_angles)
        temp = p_ind == k; % Select detections for angle k
        if sum(temp) == 0
            continue
        else
            M_all(:, temp) = ones(4, sum(temp)) .* [ cos( p_angles(k) ).^2; ...
                                                     sin( p_angles(k) ).^2; ...
                                                     0; ...
                                                     sin( p_angles(k) ).*cos( p_angles(k) )];
        end
    end

    %% Extract intensities for each polarization channel
    I0   = (2 * sqrt(pi)) .* tab_traj_4x4(param_t_r0, 1:4:end) .* tab_traj_4x4(param_t_alpha, 1:4:end);
    I45  = (2 * sqrt(pi)) .* tab_traj_4x4(param_t_r0, 2:4:end) .* tab_traj_4x4(param_t_alpha, 2:4:end);
    I90  = (2 * sqrt(pi)) .* tab_traj_4x4(param_t_r0, 3:4:end) .* tab_traj_4x4(param_t_alpha, 3:4:end);
    I135 = (2 * sqrt(pi)) .* tab_traj_4x4(param_t_r0, 4:4:end) .* tab_traj_4x4(param_t_alpha, 4:4:end);

    % Combine and filter out excluded frames
    I_all = [I0; I90; I45; I135];
    I_all = I_all(:, ~ex_t);

    % Normalize intensities for each detection
    I_all_norm = I_all ./ sum(I_all, 1);

    %% Calculate K matrix using pseudo-inverse
    K_all = (pinv(M_all')*I_all_norm')';

end
