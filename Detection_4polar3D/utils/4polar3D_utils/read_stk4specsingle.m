% =========================================================================
% Function: read_stk4specsingle
%
% Description:
%   Interface function for reading STORM images.
%   You can customize this function depending on your input format
%   (STK, TIFF, multi-file stacks, etc.) while keeping the same function name and structure.
%
% Inputs:
%   fullname_file - Full file name including path and extension
%   num           - Image index in the stack or TIFF file
%   initload      - Optional initialization parameter
%                   - 'Fin': load the entire stack into RAM
%                   - 'Deb': read from pre-loaded STK_DATA in RAM
%
% Output:
%   data          - Extracted image frame
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: December 2017
% =========================================================================

function data = read_stk4specsingle(fullname_file, num, initload)

global name_imagecor
global contourIm
global verticesleft
global rectveL
global verticesright
global rectveR

%% Compatibility: replace backslashes with slashes for Windows and Linux
fullname_file(fullname_file == '\') = '/';

%% Fast implementation (17/06/2013) - Preload everything
ind100 = mod(num, 100);

if ind100 == 0
    ind100 = 100;
    numdata = (floor(num / 100) - 1) * 100 + 1;
else
    numdata = (floor(num / 100)) * 100 + 1;
end

%% Check if the data block is already loaded in the workspace
A = evalin('base', ['exist(''data' num2str(numdata) ''');']);

if A == 0
    filename = strcat(name_imagecor, num2str(numdata));
    m = matfile(filename);
    data1 = m.data;

    [data1, ~, ~, ~, ~] = alignIm(data1, contourIm, verticesleft, rectveL, verticesright, rectveR);
    assignin('base', ['data' num2str(numdata)], data1);
end

%% Extract the requested frame
data1 = evalin('base', ['data' num2str(numdata)]);
data = data1(:, :, ind100);

%% Clear data block when the last frame is read
if ind100 == 100
    evalin('base', ['clear data' num2str(numdata)])
end

end % function
