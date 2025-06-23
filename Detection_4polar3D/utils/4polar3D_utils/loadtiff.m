function img = loadtiff(varargin)
% =========================================================================
% Function: loadtiff
%
% Description:
%   Loads a multi-page TIFF file into a 3D array. If no input argument is
%   provided, the function opens a file selection dialog. The output is 
%   a numeric array of size (height, width, n_slices).
%
% Inputs:
%   varargin - Optional input:
%              fname (char): full path to the TIFF file.
%              If not provided, a file selection dialog is opened.
%
% Output:
%   img - (numeric array) 3D array containing the image stack.
%
% Example usage:
%   stack = loadtiff('my_stack.tif');
%   stack = loadtiff();  % Opens dialog to choose file
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel  
%   Miguel Sison              - Institut Fresnel  
%   Cesar Valades-Cruz        - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

    if nargin ~= 0    
        if nargin > 1
            warndlg('Wrong number of inputs!');
            return
        end
        fname = varargin{1};  % Use provided filename
    else
        fname = uigetimagefile;  % Open file dialog
    end

    % Get metadata from TIFF
    img_inf = imfinfo(fname);

    % Read each image slice into a 3D array
    for i = 1:numel(img_inf)
        img(:,:,i) = imread(fname, i);
    end

end
