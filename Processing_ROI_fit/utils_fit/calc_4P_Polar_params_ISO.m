%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: calc_4P_Polar_params_ISO
%
% Description:
%   Calculate 4-polarization parameters from raw bead data using a given
%   calibration matrix (K). This function extracts positional, intensity,
%   and polarization properties for each localization.
%
% Inputs:
%   - originalname: path to the .mat file containing 4polar3D data
%   - K: K matrix
%
% Outputs:
%   - data: structure containing:
%       * Positions, intensities, polar factors, etc.
%
% Authors: Charitra S. Senthil Kumar, Miguel Sison, Cesar Valades-Cruz
% Date: March 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = calc_4P_Polar_params_ISO(originalname, K)

% Define theoretical lambda values from delta 0 to pi (0 to 180 deg)
nD = 2500; % Number of delta points
delta_n = linspace(0.07/pi, 1, nD) * pi;
l_iso = lambda_iso(delta_n); % lambda isotropic

%% Load data
%     disp(['Processing... ' originalname(1:end-4) '...']);
load(originalname, ...
    'tab_traj_4x4', ...
    'param_t_j', 'param_t_i', 'param_t_r0', ...
    'param_t_alpha', 'param_t_sig2', 'param_t_sig2_b', ...
    'tab_trajcoupled', 'tab_traj_allcoupled', 'param_t_tps');

% Define trajectories
tab_traj_dg = tab_traj_4x4;

tab_traj_0   = tab_traj_dg(:, 1:4:end); % 0° channel
tab_traj_45  = tab_traj_dg(:, 2:4:end); % 45° channel
tab_traj_90  = tab_traj_dg(:, 3:4:end); % 90° channel
tab_traj_135 = tab_traj_dg(:, 4:4:end); % 135° channel

% Select all trajectories (in the analysis zone)
pos_traj_in_zone = 1:size(tab_traj_0, 2);

%% Extract XY positions from each polarization channel
data.X(:,1) = tab_traj_0(param_t_j, pos_traj_in_zone);
data.Y(:,1) = tab_traj_0(param_t_i, pos_traj_in_zone);

data.X(:,2) = tab_traj_45(param_t_j, pos_traj_in_zone);
data.Y(:,2) = tab_traj_45(param_t_i, pos_traj_in_zone);

data.X(:,3) = tab_traj_90(param_t_j, pos_traj_in_zone);
data.Y(:,3) = tab_traj_90(param_t_i, pos_traj_in_zone);

data.X(:,4) = tab_traj_135(param_t_j, pos_traj_in_zone);
data.Y(:,4) = tab_traj_135(param_t_i, pos_traj_in_zone);

%% Extract frame index (time) for each localization
frame = max([tab_traj_0(param_t_tps, pos_traj_in_zone)', ...
    tab_traj_45(param_t_tps, pos_traj_in_zone)', ...
    tab_traj_90(param_t_tps, pos_traj_in_zone)', ...
    tab_traj_135(param_t_tps, pos_traj_in_zone)'], [], 2)';

data.T(:,1) = frame;

%% Extract intensities and PSF radii for each channel
I0   = (2 * sqrt(pi)) .* tab_traj_0(param_t_r0, pos_traj_in_zone) .* tab_traj_0(param_t_alpha, pos_traj_in_zone);
I45  = (2 * sqrt(pi)) .* tab_traj_45(param_t_r0, pos_traj_in_zone) .* tab_traj_45(param_t_alpha, pos_traj_in_zone);
I90  = (2 * sqrt(pi)) .* tab_traj_90(param_t_r0, pos_traj_in_zone) .* tab_traj_90(param_t_alpha, pos_traj_in_zone);
I135 = (2 * sqrt(pi)) .* tab_traj_135(param_t_r0, pos_traj_in_zone) .* tab_traj_135(param_t_alpha, pos_traj_in_zone);

radius_0   = tab_traj_0(param_t_r0, pos_traj_in_zone);
radius_45  = tab_traj_45(param_t_r0, pos_traj_in_zone);
radius_90  = tab_traj_90(param_t_r0, pos_traj_in_zone);
radius_135 = tab_traj_135(param_t_r0, pos_traj_in_zone);

sigma_0    = sqrt(tab_traj_0(param_t_sig2, pos_traj_in_zone));
sigma_45   = sqrt(tab_traj_45(param_t_sig2, pos_traj_in_zone));
sigma_90   = sqrt(tab_traj_90(param_t_sig2, pos_traj_in_zone));
sigma_135  = sqrt(tab_traj_135(param_t_sig2, pos_traj_in_zone));

%% Organize intensity and radius data
data.It(:,1) = I0 + I45 + I90 + I135; % Total intensity
data.I(:,1:4) = [I0; I45; I90; I135]'; % Intensity per channel

data.radius_mean(:,1) = mean([radius_0; radius_45; radius_90; radius_135], 1);
data.radius_max(:,1) = max([radius_0; radius_45; radius_90; radius_135], [], 1);
data.radius(:,1:4) = [radius_0; radius_45; radius_90; radius_135]';

data.sigma_mean(:,1) = mean([sigma_0; sigma_45; sigma_90; sigma_135], 1);
data.sigma_max(:,1) = max([sigma_0; sigma_45; sigma_90; sigma_135], [], 1);
data.sigma(:,1:4) = [sigma_0; sigma_45; sigma_90; sigma_135]';

% Estimate localization noise
data.noise = sqrt( ...
    [tab_traj_0(param_t_sig2_b, pos_traj_in_zone); ...
    tab_traj_45(param_t_sig2_b, pos_traj_in_zone); ...
    tab_traj_90(param_t_sig2_b, pos_traj_in_zone); ...
    tab_traj_135(param_t_sig2_b, pos_traj_in_zone)]');

%% Calculate polarization factors
% Solve using pseudo-inverse of K
M = (pinv(K) * [I0; I90; I45; I135])';
A2 = abs(sum(M(:,1:3), 2));

% Calculate polarization components
Pxy = (M(:,1) - M(:,2)) ./ A2;
Puv = 2 * M(:,4) ./ A2;
Pz = M(:,3) ./ A2;
%%%%%%%%%%%%%%% DATA  STRUCTURE %%%%%%%%%%%%%%%
%   Columns:    (1) ISO, (2) TIRF 2, (3) TIRF S
%   Rows:       Indivudual PSFs

% Calculate rho (orientation), delta, and eta (polarization parameters)
rho(:,1) = 0.5 * mod(atan2(Puv, Pxy), 2 * pi);

% Isotropic illumination analysis
lambda_3    =   Pz + sqrt(Puv.^2 + Pxy.^2);
lambda      =   (1 - lambda_3)/2;

filt_l(:,1) = (lambda > max(l_iso)) | (lambda < min(l_iso));

delta(:,1) = 2 * acos(0.5 * (-1 + sqrt(9 - 24 * lambda)));
eta(:,1) = acos(sqrt((Pz - lambda) ./ (1 - 3 * lambda)));

filt_c(:,1) = any(imag([rho(:,1), delta(:,1), eta(:,1)]) ~= 0, 2);

%% Store final results in output structure
data.rho = rho;
data.eta = eta;
data.delta = delta;
data.filt_l = filt_l; % Filter: lambda out of isotropic range
data.filt_c = filt_c; % Filter: complex (invalid) solutions
data.lambdas = lambda;
data.Pxy = Pxy;
data.Puv = Puv;
data.Pz = Pz;
data.A2 = A2;
data.M = M;
data.K = K;

% Mark excluded points based on complex solutions
data.ex = any(imag([delta(:,1), eta(:,1)]) ~= 0, 2);

%     disp('4P calc done!')

end
