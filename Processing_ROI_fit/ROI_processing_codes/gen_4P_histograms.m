% =========================================================================
% Function Name: gen_4P_histograms.m
%
% Description:
%   This function generates multiple histograms and 2D distributions to visualize
%   the 4polar3D parameters (rho, eta, delta, and photon counts) for a selected ROI.
%
%   The function:
%     - Extracts the localization data within the selected ROI.
%     - Creates polar and bar histograms for rho, eta, and delta distributions.
%     - Generates 2D histograms (heatmaps) of photon counts vs rho, eta, and delta.
%     - Displays the selected ROI mask in a separate figure.
%
% Instructions:
%   - Call this function to generate histograms for a specific ROI.
%
% Notes:
%   - The function uses the plasma colormap for eta and jet colormap for delta.
%   - The histograms adapt their bin ranges to the photon intensity distribution.
%
% Inputs:
%   - data: Structure containing the localization and parameter data.
%   - roi_struct: Structure containing ROI masks and list representations.
%   - roi_num: Index of the ROI to be analyzed (zero-based).
%   - nB: Number of bins to be used for histograms.
%
% Outputs:
%   - fig_h: Figure handle of the generated histogram panel.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function fig_h = gen_4P_histograms(data, roi_struct, roi_num, nB)

fs = 12; % Font size for plots

% Select localizations that belong to the selected ROI
roi_select = any(roi_struct.list_form(:, roi_num + 1), 2);

% Extract rho, eta, delta, and total photon counts for the ROI
rho = data.rho(roi_select);
eta = data.eta(roi_select);
delta = data.delta(roi_select);
Iph = data.I_tot(roi_select);

% Define photon count binning limits based on percentiles
I_mn = min(floor(prctile(Iph, 2.5) / 10) * 10, 1000);
I_mx = max(ceil(prctile(Iph, 97.5) / 100) * 100, 3000);
Ibin_edge = linspace(I_mn, I_mx, nB + 1);
rbin_edge = linspace(0, 180, nB + 1);
ebin_edge = linspace(0, 90, nB + 1);
dbin_edge = linspace(0, 180, nB + 1);

% Create main figure with 2x3 layout
fig_h = figure('Position', [100 200 1050 610]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% Plot polar histogram of rho
nexttile;
polarhistogram(rho, deg2rad(rbin_edge * 2), 'FaceAlpha', 0.1, ...
    'LineWidth', 1.5, 'LineStyle', '-', 'EdgeColor', [0 0.4470 0.7410] .* 0.75);
hold on
polarhistogram(rho + pi, deg2rad(rbin_edge * 2), 'FaceAlpha', 0.1, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'EdgeColor', [0.8500 0.3250 0.0980] * 0.75);
box off
set(gca, 'Color', 'none')
set(gca, 'FontSize', fs)

% Plot histogram of eta
nexttile
[b_counts, b_edge] = histcounts(rad2deg(eta), ebin_edge);
b = bar(b_edge(1:end - 1) + abs(diff(b_edge)) / 2, b_counts, ...
    'BarWidth', 1, 'FaceAlpha', 0.95, 'EdgeAlpha', 0.95, ...
    'LineWidth', 1.5, 'FaceColor', 'flat');
hclr = plasma(nB); % Use plasma colormap
for k = 1:nB
    b.CData(k, :) = hclr(k, :);
end
box off
set(gca, 'Color', 'none')
set(gca, 'FontSize', fs)

% Plot histogram of delta
nexttile
[b_counts, b_edge] = histcounts(rad2deg(delta), dbin_edge);
b = bar(b_edge(1:end - 1) + abs(diff(b_edge)) / 2, b_counts, ...
    'BarWidth', 1, 'FaceAlpha', 0.95, 'EdgeAlpha', 0.95, ...
    'LineWidth', 1.5, 'FaceColor', 'flat');
hclr = jet(nB); % Use jet colormap
for k = 1:nB
    b.CData(k, :) = hclr(k, :);
end
box off
set(gca, 'Color', 'none')
set(gca, 'FontSize', fs)

% Plot 2D histogram: Photon counts vs rho
nexttile
img_hst = histcounts2(Iph, rad2deg(rho), Ibin_edge, rbin_edge);
imagesc(img_hst, 'YData', Ibin_edge(1:end - 1) + diff(Ibin_edge) / 2, ...
    'XData', rbin_edge(1:end - 1) + diff(rbin_edge) / 2, 'AlphaData', (img_hst > 0) * 0.9 + 0.1);
axis square
colorbar
set(gca, 'YDir', 'normal')
ylabel('Photon counts', 'Rotation', 90)
xlabel('\rho', 'Rotation', 0)
set(gca, 'YTickLabelRotation', 45)
set(gca, 'Color', 'none')
set(gca, 'FontSize', fs)

% Plot 2D histogram: Photon counts vs eta
nexttile
img_hst = histcounts2(Iph, rad2deg(eta), Ibin_edge, ebin_edge);
imagesc(img_hst, 'YData', Ibin_edge(1:end - 1) + diff(Ibin_edge) / 2, ...
    'XData', ebin_edge(1:end - 1) + diff(ebin_edge) / 2, 'AlphaData', (img_hst > 0) * 0.9 + 0.1);
axis square
colorbar
set(gca, 'YDir', 'normal')
xlabel('\eta', 'Rotation', 0)
set(gca, 'YTick', [])
set(gca, 'Color', 'none')
set(gca, 'FontSize', fs)

% Plot 2D histogram: Photon counts vs delta
nexttile
img_hst = histcounts2(Iph, rad2deg(delta), Ibin_edge, dbin_edge);
imagesc(img_hst, 'YData', Ibin_edge(1:end - 1) + diff(Ibin_edge) / 2, ...
    'XData', dbin_edge(1:end - 1) + diff(dbin_edge) / 2, 'AlphaData', (img_hst > 0) * 0.9 + 0.1);
axis square
colorbar
set(gca, 'YDir', 'normal')
xlabel('\delta', 'Rotation', 0)
set(gca, 'YTick', [])
set(gca, 'Color', 'none')
set(gca, 'FontSize', fs)

% Display the mask of the selected ROI
figure;
temp = bwlabel(any(roi_struct.mask_form(:, :, roi_num + 1), 3));
imagesc(temp);
axis image
set(gca, 'YDir', 'normal')

%     setappdata(gcf, 'rho_roi', rho);
%     setappdata(gcf, 'eta_roi', eta);
%     setappdata(gcf, 'delta_roi', delta);
%     setappdata(gcf, 'It_roi', It);
%     setappdata(gcf, 'Iph_roi', Iph);
end
