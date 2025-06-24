% =========================================================================
% Function Name: create_single_4Psticks_image.m
%
% Description:
%   This function creates a single 4polar3D stick plot by overlaying colored
%   sticks (orientation lines) on top of an intensity image.
%
%   The function:
%     - Displays the intensity image using tiled layout.
%     - Draws orientation sticks based on the provided positions, angles, and colors.
%     - Uses a transparent axis for the overlay to keep the image visible.
%     - Links axes and manages color display.
%
% Instructions:
%   - Call this function to generate the stick plot for a selected image.
%   - Stick parameters (length, width, transparency, and color) must be provided.
%
% Notes:
%   - This function is part of the 4polar3D MATLAB visualization pipeline.
%   - The draw_line_patch function is used to draw individual sticks.
%
% Inputs:
%   - fig: Figure handle where the plot will be generated.
%   - img: Background image to display.
%   - xdata, ydata: X and Y scale data for the image.
%   - x, y: X and Y positions of the sticks.
%   - rho: Stick orientation angles (radians).
%   - c_line: Color assignment for each stick.
%   - sl: Stick half-length.
%   - sw: Stick half-width.
%   - sa: Stick transparency (alpha value).
%
% Outputs:
%   - ax1: Axis handle for the background image.
%   - ax_s: Axis handle for the stick overlay.
%   - fig: Figure handle of the generated plot.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function [ax1, ax_s, fig] = create_single_4Psticks_image(fig, img,  xdata, ydata, x,  y,  rho, c_line, sl, sw, sa)

%     figure('Position',[100 100 1620 520]);
figure(fig);

% Create a 1x1 tiled layout with no padding for clean display
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'loose');

% Display the background image in the first axis
ax1 = nexttile;
imagesc(ax1,img,"XData",xdata,'YData',ydata)
set(ax1,'Colormap',gray)    % Set colormap to grayscale
set(ax1,'YDir','normal')    % Correct Y-axis direction
set(ax1,'YTick',[])         % Hide Y-axis ticks
set(ax1,'XTick',[])         % Hide X-axis ticks
axis image                  % Keep aspect ratio

% Create a second transparent axis to draw the sticks
ax_s = axes;
draw_line_patch(x,y,rho,sl,sw,c_line,sa); % Draw all sticks at once

%% Optional: sequential stick drawing by color range
% rc = linspace(0,180,NC+1);    % Example binning for rho or delta
% rc = linspace(0,90,NC+1);     % Example binning for eta
% for k = 1:NC
%     temp = (c_line > rc(k)) & (c_line <= rc(k+1)); % Select sticks in the current color bin
%     draw_line_patch(x(temp),y(temp),rho(temp),sl,sw,c_line(temp),a_st);
%     hold on
% end
%%

% Configure the overlay axis
set(ax_s,'Colormap',jet);   % Set stick colormap to jet
set(ax_s,'Color','none');   % Make axis background transparent

% Add a colorbar to the overlay axis
colorbar(ax_s,'eastoutside');
axis image;                  % Keep aspect ratio

% Link the two axes to synchronize zoom and pan
linkaxes([ax1,ax_s],'xy')
ax_s.PositionConstraint = 'outerposition';
ax_s.InnerPosition = ax1.InnerPosition;

end
