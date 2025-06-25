% =========================================================================
% Function Name: ppm8write.m
%
% Description:
%   This is a debug function to save an RGB image in PPM format (8-bit per channel).
%
%   The function:
%     - Writes a PPM (Portable Pixmap) image file using 8-bit RGB data.
%     - The RGB input data should be in the range 0-255.
%
% Instructions:
%   - Call this function with the file name and the red, green, and blue image matrices.
%
% Notes:
%   - The input image matrices must be the same size.
%   - Data must be scaled between 0 and 255.
%   - The function writes the file in binary (big-endian format).
%
% Inputs:
%   - fname: File name for the output PPM file.
%   - r: Red channel matrix.
%   - v: Green channel matrix.
%   - b: Blue channel matrix.
%
% Outputs:
%   - None. The function writes the PPM image to disk.
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: December 2017
% =========================================================================

%% Debug function: write 8-bit PPM image
%% Data expected in the range 0-255 for RGB channels

function ppm8write(fname, r, v, b)

[L, C] = size(r);

r = r.';
v = v.';
b = b.';
data = reshape([r(:).' ; v(:).' ; b(:).'], 1, L*C*3) ;

% Alternative reshape versions (commented)
% data = reshape([(r')(:)' ; (v')(:)' ; (b')(:)'], 1, L*C*3);
% data = reshape([((r.')(:)).' ; ((v.')(:)).' ; ((b.')(:)).'], 1, L*C*3);

fid = fopen(fname, 'wb');

%% Header
fprintf(fid, 'P6\n');
fprintf(fid, '%d %d\n', C, L);
fprintf(fid, '255\n');

%% Data
fwrite(fid, data, 'uint8', 0, 'ieee-be');

fclose(fid);

end %function
