% =========================================================================
% Function Name: create_4P_images.m
%
% Description:
%   This function creates 4polar3D parameter maps (rho, eta, delta, intensity)
%   from the provided localization data.
%
%   The function:
%     - Bins the localizations into a 2D histogram based on provided X and Y edges.
%     - Computes unweighted and intensity-weighted averages for eta, delta, and rho.
%     - Builds structured output images for each parameter.
%
% Instructions:
%   - Call this function with localization positions, intensities, and parameter arrays.
%   - The function will automatically bin the data and compute the maps.
%
% Notes:
%   - rho values are circularly averaged using the circ_mean function.
%   - The binning is performed using histcounts2.
%   - This function is part of the 4polar3D MATLAB analysis pipeline.
%
% Inputs:
%   - x_loc: X positions of the localizations.
%   - y_loc: Y positions of the localizations.
%   - It: Intensity values for each localization.
%   - rho: Rho parameter for each localization.
%   - eta: Eta parameter for each localization.
%   - delta: Delta parameter for each localization.
%   - x_edges: Bin edges along the X axis.
%   - y_edges: Bin edges along the Y axis.
%
% Outputs:
%   - image_all: Structure containing binned images for intensity, rho, eta, and delta.
%       .img_strm: 2D histogram of localization counts and summed intensities.
%       .img_rho:  Unweighted and weighted rho maps.
%       .img_eta:  Unweighted and weighted eta maps.
%       .img_delta: Unweighted and weighted delta maps.
%       .Ybinned: Bin indices along Y for each localization.
%       .Xbinned: Bin indices along X for each localization.
%       .x_loc: Original X positions.
%       .y_loc: Original Y positions.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function image_all = create_4P_images(x_loc, y_loc, It, rho, eta, delta, x_edges, y_edges)

% Bin the localizations into a 2D histogram
% temp is the localization count per bin
% bnR, bnC are the bin indices of each localization
[temp, ~, ~, bnR, bnC]  = histcounts2(y_loc, x_loc, y_edges, x_edges);

% Define the size of the output images
img_px = [numel(y_edges) numel(x_edges)] - 1;

% Initialize empty arrays for rho, eta, delta
img_rho     = NaN(img_px(1), img_px(2), 2);
img_eta     = NaN(img_px(1), img_px(2), 2);
img_delta   = NaN(img_px(1), img_px(2), 2);

% Create intensity map (stream image)
% Channel 1: number of localizations per bin
% Channel 2: summed intensity per bin
img_strm(:,:,1) = temp;
img_strm(:,:,2) = zeros(img_px(1), img_px(2));

% Find all unique bins that contain localizations
temp_coord = unique([bnR, bnC], 'rows');

% Number of occupied bins
N = size(temp_coord, 1);

% Shuffle bin order (optional, just for processing order)
temp = randperm(N);

% Loop over each occupied bin
for kc = 1:N

    % Get current bin coordinates
    r = temp_coord(temp(kc), 1);
    c = temp_coord(temp(kc), 2);

    % Find all localizations belonging to the current bin
    temp_filter = (bnR == r) & (bnC == c);

    % Calculate weights based on intensity (normalize by mean intensity)
    W = It(temp_filter) ./ mean(It(temp_filter));

    % Unweighted averages
    img_eta(r, c, 1)    = mean( eta(temp_filter) );
    img_delta(r, c, 1)  = mean( delta(temp_filter) );
    img_rho(r, c, 1)    = mod( circ_mean( rho(temp_filter).* 2 ), 2*pi ) / 2;

    % Weighted averages
    img_eta(r, c, 2)    = mean( eta(temp_filter) .* W );
    img_delta(r, c, 2)  = mean( delta(temp_filter) .* W );
    img_rho(r, c, 2)    = mod( circ_mean( rho(temp_filter).* 2, W ), 2*pi ) / 2;

    % Sum intensity for this bin
    img_strm(r, c, 2) = sum( It(temp_filter) );

end

% Store all results in output structure
image_all.img_strm  = img_strm;   % Intensity maps
image_all.img_rho   = img_rho;    % Rho parameter maps
image_all.img_eta   = img_eta;    % Eta parameter maps
image_all.img_delta = img_delta;  % Delta parameter maps
image_all.Ybinned   = bnR;        % Bin indices for each localization (Y)
image_all.Xbinned   = bnC;        % Bin indices for each localization (X)
image_all.x_loc     = x_loc;      % Original X positions
image_all.y_loc     = y_loc;      % Original Y positions
end

%% Original function from Circular Statistics Toolbox (unchanged)

function [mu, ul, ll] = circ_mean(alpha, w, dim)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 3
    dim = 1;
end

if nargin < 2 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
    t = circ_confmean(alpha,0.05,w,[],dim);
    ul = mu + t;
    ll = mu - t;
end
end