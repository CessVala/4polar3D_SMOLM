%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: calc_Polar_4I_ISO
%
% Description:
%   Calculate 4-polarization parameters from measured intensities using a
%   provided K matrix.
%
% Inputs:
%   - I: 4-polarization intensity matrix
%   - K: K matrix
%
% Outputs:
%   - data: structure containing polarization parameters and filtering info
%
% Authors: Charitra S. Senthil Kumar, Miguel Sison, Cesar Valades-Cruz
% Date: March 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = calc_Polar_4I_ISO(I, K)

%% Ensure all intensity values are non-negative
I = I .* (I >= 0); % Set negative intensities to zero

% Define theoretical lambda values from delta 0 to pi (0 to 180 deg)
nD = 2500;
delta_n = linspace(0.07/pi, 1, nD) * pi;
l_iso = lambda_iso(delta_n); % lambda isotropic

%% Calculate M matrix using pseudo-inverse of K
%     M = (K \ [I(:,1)'; I(:,3)'; I(:,2)'; I(:,4)'])';
M = (pinv(K) * I)';
A2 = sum(M(:,1:3), 2);

%% Calculate polarization components
Pxy = (M(:,1) - M(:,2)) ./ A2;
Puv = 2 * M(:,4) ./ A2;
Pz  = M(:,3) ./ A2;

%%%%%%%%%%%%%%% DATA  STRUCTURE %%%%%%%%%%%%%%%
%   Columns:    (1) ISO, (2) TIRF 2, (3) TIRF S
%   Rows:       Indivudual PSFs

rho(:,1)    =   0.5*mod(atan2(Puv,Pxy),2*pi);

% Isotropic illumination analysis
lambda_3    =   Pz + sqrt(Puv.^2 + Pxy.^2);
lambda      =   (1 - lambda_3)/2;

% Filter for lambda values out of isotropic range
filt_l(:,1) = (lambda > max(l_iso)) | (lambda < min(l_iso));

% Calculate delta
%     delta(:,1) = 2 * acos(0.5 * abs(-1 + sqrt(12 * lambda_3 - 3)));
delta(:,1) = 2 * acos(0.5 * (-1 + sqrt(9 - 24 * lambda)));

% Calculate eta
eta(:,1) = acos(sqrt((Pz - lambda) ./ (1 - 3 * lambda)));

% Filter out complex (invalid) results
filt_c(:,1) = any(imag([rho(:,1), delta(:,1), eta(:,1)]) ~= 0, 2);

%% Store results in the output structure
data.I = I;                % Raw intensities
data.It = sum(I, 1);       % Total intensity
data.rho = rho;
data.eta = eta;
data.delta = delta;
data.filt_l = filt_l;      % Filter: lambda out of isotropic range
data.filt_c = filt_c;      % Filter: complex results
data.lambdas = lambda;     % Lambda values
data.Pxy = Pxy;
data.Puv = Puv;
data.Pz = Pz;
data.A2 = A2;
data.M = M;                % M matrix
data.K = K;                % K matrix

end
