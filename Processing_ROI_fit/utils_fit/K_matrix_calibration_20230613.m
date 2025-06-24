%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: K_matrix_calibration.m
%
% Description:
%   Calibrate the K matrix using fluorescent beads. The script:
%     1. Loads physical and optical parameters
%     2. Calculates K from bead measurements
%     3. Fits Fourier series to measured intensities
%     4. Optimizes the BFP factor
%     5. Calculates the theoretical K matrix
%     6. Applies dipole correction
%     7. Saves the final corrected K matrix
%
% Authors: Charitra S. Senthil Kumar, Miguel Sison
% Date: March 2023
%
% Notes:
%   
%   - Requires supporting functions: K_from_beads, Polarizer_4P_Intensity, 
%     K_calc_new, K_Dipole_corr, Kzz_calib
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%% Define optical and physical parameters
NA          = 1.45;         % Numerical aperture of the objective
n_glass     = 1.515;        % Refractive index of glass/oil
n_sample    = 1.33;         % Refractive index of the sample (water)
lambda      = 0.600;        % Wavelength in micrometers
z_dipole    = 0.0;          % Dipole position (interface or z = 0)

SAF_flag    = 0;            % Exclude SAF (Supercritical Angle Fluorescence)
                            % If excluded, effective NA is limited to n_sample
                            % in this case NA_eff = 1.33

bfp_factor0  = 0.8;         % Initial BFP (Back Focal Plane) reduction factor

%% Calculate calibrated K (in plane components only)
p_angles    = deg2rad(0:10:180);    % Polarizer angles in radians
start_frame = 11;                   % First useful frame in tiff stack
n_frame     = 20;                   % Number of frames per polarizer angle

% Compute initial K from bead measurements
[K_fbeads, I_all, M_all, p_ind] = K_from_beads(p_angles, start_frame, n_frame);

%% Fit Fourier series to measured intensities for each channel
% Using the model: A0 + A2*cos(2*alpha) + B2*sin(2*alpha)

fun_I = @(f) sum((I_all(1,:)./sum(I_all,1) - (f(1) + f(2)*cos(2*p_angles(p_ind)) + f(3)*sin(2*p_angles(p_ind)))).^2);
C000  = fmincon(fun_I, [0 0 0], [], [], [], [], -[1 1 1], [1 1 1]);

fun_I = @(f) sum((I_all(2,:)./sum(I_all,1) - (f(1) + f(2)*cos(2*p_angles(p_ind)) + f(3)*sin(2*p_angles(p_ind)))).^2);
C090  = fmincon(fun_I, [0 0 0], [], [], [], [], -[1 1 1], [1 1 1]);

fun_I = @(f) sum((I_all(3,:)./sum(I_all,1) - (f(1) + f(2)*cos(2*p_angles(p_ind)) + f(3)*sin(2*p_angles(p_ind)))).^2);
C045  = fmincon(fun_I, [0 0 0], [], [], [], [], -[1 1 1], [1 1 1]);

fun_I = @(f) sum((I_all(4,:)./sum(I_all,1) - (f(1) + f(2)*cos(2*p_angles(p_ind)) + f(3)*sin(2*p_angles(p_ind)))).^2);
C135  = fmincon(fun_I, [0 0 0], [], [], [], [], -[1 1 1], [1 1 1]);

% Extract Fourier coefficients
C = [C000; C090; C045; C135];
A0beads = C(:,1); 
A2beads = C(:,2); % Cosine component
B2beads = C(:,3); % Sine component
clear C000 C090 C045 C135 C

% check phase shift or angular offset of polarized fluorescent beads
% measurements. Ideal case:
% [0 ; 90; 45; 135]
disp(rad2deg(mod(atan2(B2beads,A2beads)/2,pi)))

%% Optimize BFP reduction factor based on intensity of rotating polarizer
[~, NAred_ch] = min( mean(K_fbeads(1:2:end,1:2) + K_fbeads(2:2:end,1:2),2) );

% Estimate extinction efficiency
tempA = min(sqrt(A2beads.^2 + B2beads.^2)./A0beads, 1);

% Estimate phase offset from expected values
tempB = (mod(atan2(B2beads, A2beads)/2, pi)) - deg2rad([1.5491; 91.8193; 46.2095; 133.5999]);  

% Define fitting function for polarizer intensity
fun_K = @(f) sum( ( I_all./sum(I_all,1) - Polarizer_4P_Intensity(n_sample, n_glass, NA, f(1), p_angles(p_ind)-mean(tempB), z_dipole, lambda, SAF_flag, NAred_ch, [f(2); f(3); f(4); f(5)]) ).^2, 'all' );

% Fit BFP reduction factor and polarization extinction efficiency
factor_new = fmincon(fun_K, [bfp_factor0 tempA'], [], [], [], [], [0.5 0.8 0.8 0.8 0.8], [1 1 1 1 1]);

% Extract optimized parameters
bfp_factor = factor_new(1)
A = factor_new(2:end)'  % Polarization extinction efficiencies

clear factor_new

%% Check the BFP fitting result by plotting
figure;
plot(rad2deg(p_angles(p_ind)), (I_all./sum(I_all,1))', 'o', 'LineWidth',1.5) % Experimental data
hold on
set(gca,'ColorOrderIndex',1)
plot(rad2deg(linspace(0,pi,250)), Polarizer_4P_Intensity(n_sample, n_glass, NA, bfp_factor, linspace(0,pi,250), z_dipole, lambda, SAF_flag, NAred_ch, A)', '-', 'LineWidth', 1.5) % Fitted curve
set(gca,'Color','none')
grid on
box off

%% Calculate the theoretical K matrix
K_theo = K_calc_new(bfp_factor, n_glass, n_sample, NA, SAF_flag);
K_theo = K_theo./sum(K_theo(:,1)); % Normalize

% Temporary Kzz component
Kzz = K_theo(:,3);

%% Apply dipole correction following 'Fourier series decomposition' based on bfp_reduction factor
[K_corr, corr_factor] = K_Dipole_corr(K_fbeads, bfp_factor, n_glass, n_sample, NA, SAF_flag);

%% Calculate Kzz from isotropic fluorescent beads
[Kzz, ~, I_beads] = Kzz_calib(K_corr, K_theo);

% Final corrected K matrix
K_final = K_corr;
K_final(:,3) = Kzz;

% Save final K matrix
[fname,fpath] = uiputfile('*.mat');
save(fullfile(fpath,fname),'K_final', "bfp_factor", "-v7.3")

%% Optional: Windowed K calculation (commented block)
% These sections compute K using a detection window based on PSF size.
% It is currently disabled.

% % % Windowed K: summing intensities inside detection window
% img_pxs = 0.020; % Pixel size for K computation
% k_dim = 2^11;    % Image dimension
% cam_px_sz = 0.130; % Camera pixel size in microns
% win_pix = 11;    % Detection window size in pixels

% % Compute windowed K
% K_win1 = K_PSF_win(win_pix, cam_px_sz, bfp_factor1, NA, n_sample, n_glass, img_pxs, k_dim, z_dipole, lambda, SAF_flag);

% % Apply dipole correction for windowed K
% [K_corr1_win, corr_factor1_win] = K_Dipole_corr_psf_win(K_fbeads, bfp_factor1, n_glass, n_sample, NA, img_pxs, k_dim, win_pix, cam_px_sz, SAF_flag);

% % Calculate Kzz for windowed K
% Kzz1_win = Kzz_calib(K_corr1_win, K_win1, I_beads);

% % Final windowed K matrix
% K_test1_win = K_corr1_win;
% K_test1_win(:,3) = Kzz1_win;
