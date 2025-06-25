% =========================================================================
% Function Name: alignIm.m
%
% Description:
%   This function allows the user to manually select two elliptical regions
%   in a multi-frame image, aligns and crops these regions, and outputs
%   the combined cropped images side by side.
%
%   The function:
%     - If alignment positions are not provided, displays the image and
%       prompts the user to select the ellipses.
%     - Extracts and crops the selected regions.
%     - Aligns and pads the cropped images with a defined contour.
%     - Combines the two selected image regions into a single output image.
%
% Instructions:
%   - If called with only data1 and contourIm, the user will be asked to
%     interactively select the ellipses.
%   - If called with all arguments, the provided selection is directly used.
%
% Inputs:
%   - data1: Multi-frame image (3D matrix).
%   - contourIm: Padding size added around the cropped regions.
%   - verticesleft: Coordinates of the first ellipse (optional).
%   - rectveL: Bounding box of the first ellipse (optional).
%   - verticesright: Coordinates of the second ellipse (optional).
%   - rectveR: Bounding box of the second ellipse (optional).
%
% Outputs:
%   - dataselfull: Final aligned and combined image.
%   - verticesleft: Coordinates of the first ellipse.
%   - rectveL: Bounding box of the first ellipse.
%   - verticesright: Coordinates of the second ellipse.
%   - rectveR: Bounding box of the second ellipse.
%
% Author:
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

function [dataselfull, verticesleft, rectveL, verticesright, rectveR] = alignIm(data1, contourIm, verticesleft, rectveL, verticesright, rectveR)

% If ellipses are not provided, ask user to select them interactively
if (nargin <= 2)
    h1 = figure(1); % Display the first frame
    imagesc(data1(:, :, 1));
    colormap('gray');
    axis('image');
    contrasth = imcontrast(h1); % Add contrast adjustment tool

    % Select first ellipse
    h = imellipse;
    verticesleft = wait(h); % Wait for user selection
    rectveL = round(getPosition(h)); % Get bounding box

    % Prompt for second ellipse
    msgbox('Select the second ellipse', 'Hi', 'modal');
    verticesright = wait(h); % Wait for second selection
    rectveR = round(getPosition(h)); % Get bounding box
end

% Create mask for the first ellipse
BWl = poly2mask(verticesleft(:, 1), verticesleft(:, 2), size(data1(:, :, 1), 1), size(data1(:, :, 1), 2));
BW2l = repmat(BWl, [1, 1, size(data1, 3)]); % Extend mask to all frames
data1sel = data1 .* BW2l; % Apply mask to image

% Create mask for the second ellipse
BWr = poly2mask(verticesright(:, 1), verticesright(:, 2), size(data1(:, :, 1), 1), size(data1(:, :, 1), 2));
BW2r = repmat(BWr, [1, 1, size(data1, 3)]); % Extend mask to all frames
data2sel = data1 .* BW2r; % Apply mask to image

% Crop selected regions using bounding boxes
data1selcrop = data1sel(rectveL(2):(rectveL(2) + rectveL(4)), rectveL(1):(rectveL(1) + rectveL(3)), :);
data2selcrop = data2sel(rectveR(2):(rectveR(2) + rectveR(4)), rectveR(1):(rectveR(1) + rectveR(3)), :);

% Determine maximum crop size
maxM = max(size(data1selcrop, 1), size(data2selcrop, 1));
maxN = max(size(data1selcrop, 2), size(data2selcrop, 2));

% Create padded regions with additional contour
data1selcropnew = zeros(maxM + 2 * contourIm, maxN + 2 * contourIm, size(data1selcrop, 3));
data2selcropnew = zeros(maxM + 2 * contourIm, maxN + 2 * contourIm, size(data1selcrop, 3));

% Shift cropped images to include the padding
data1selcrop = imtranslate(data1selcrop, [contourIm, contourIm], 'OutputView', 'full');
data1selcropnew(1:size(data1selcrop, 1), 1:size(data1selcrop, 2), :) = data1selcrop;

data2selcrop = imtranslate(data2selcrop, [contourIm, contourIm], 'OutputView', 'full');
data2selcropnew(1:size(data2selcrop, 1), 1:size(data2selcrop, 2), :) = data2selcrop;

% Combine both regions side by side
dataselfull = [data1selcropnew data2selcropnew];

% Fill zero-padded areas with random values to avoid pure zero regions
randmatrix = rand(size(dataselfull));
dataselfull(dataselfull == 0) = randmatrix(dataselfull == 0);

% If selection was interactive, close contrast tool and figure
if (nargin <= 2)
    close(contrasth);
    close(h1);
end

end %function
