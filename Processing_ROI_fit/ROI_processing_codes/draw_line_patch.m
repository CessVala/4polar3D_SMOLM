% =========================================================================
% Function Name: draw_line_patch.m
%
% Description:
%   This function draws a filled rectangle (patch) representing a line segment
%   with a specified orientation, length, width, color, and transparency.
%
%   The function:
%     - Calculates the four vertices of the rectangle based on the line orientation.
%     - Uses the MATLAB 'patch' function to draw the filled shape without edges.
%
% Instructions:
%   - Call this function to plot a single stick (line patch) at a specific location.
%
% Notes:
%   - This function is part of the 4polar3D visualization pipeline.
%   - The orientation angle is in radians.
%
% Inputs:
%   - x: X position of the center of the stick.
%   - y: Y position of the center of the stick.
%   - angle: Orientation angle in radians.
%   - sl: Stick half-length (distance from center to end along the main axis).
%   - sw: Stick half-width (distance from center to side along the minor axis).
%   - color: RGB color of the patch (e.g., [1 0 0] for red).
%   - alpha: Transparency level of the patch (0 = fully transparent, 1 = fully opaque).
%
% Outputs:
%   - None. The function directly draws on the current figure.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function draw_line_patch(x, y, angle, sl, sw, color, alpha)

% Compute projections along the length (sl) and width (sw) directions
xa = sl .* cos(angle);
xb = sw .* sin(angle);
ya = sl .* sin(angle);
yb = sw .* cos(angle);

% Calculate X coordinates of the four vertices of the rectangle
vx = [x - xa + xb, ... % Vertex 1
    x + xa + xb, ... % Vertex 2
    x + xa - xb, ... % Vertex 3
    x - xa - xb];    % Vertex 4

% Calculate Y coordinates of the four vertices of the rectangle
vy = [y - ya - yb, ... % Vertex 1
    y + ya - yb, ... % Vertex 2
    y + ya + yb, ... % Vertex 3
    y - ya + yb];    % Vertex 4

% Draw the patch (filled rectangle) without edge lines
patch(vx', vy', color, 'EdgeColor', 'none', 'FaceAlpha', alpha)

end
