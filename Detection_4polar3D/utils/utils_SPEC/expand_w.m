% =========================================================================
% Function Name: expand_w.m
%
% Description:
%   This function resizes a matrix to the target size (N, M) by embedding
%   the input matrix at the center of the new matrix.
%
%   The function:
%     - Takes a smaller matrix as input.
%     - Creates a larger matrix of the specified size.
%     - Places the input matrix at the center of the output matrix.
%
% Instructions:
%   - Call this function with an input matrix and the desired output size (N, M).
%
% Notes:
%   - The output matrix is initialized to zeros.
%   - The input matrix is placed exactly at the center.
%
% Inputs:
%   - in: Input matrix to be resized.
%   - N: Desired number of rows in the output matrix.
%   - M: Desired number of columns in the output matrix.
%
% Outputs:
%   - out: The resized matrix with the input centered inside.
%
% Authors:
%   Nicolas Bertaux - Institut Fresnel
%   Cesar Valades-Cruz - Institute of Hydrobiology (IHB), CAS
%
% Date: December 2017
% =========================================================================

function out = expand_w(in, N, M)

% Get size of input matrix
[N_in, M_in] = size(in);

% Initialize output matrix with zeros
out = zeros(N, M);

% Calculate center position for placement
% nc = 1+floor(N/2 - N_in/2);
% mc = 1+floor(M/2 - M_in/2);
% out(nc:(nc+N_in-1), mc:(mc+M_in-1)) = in;

nc = floor(N/2 - N_in/2);
mc = floor(M/2 - M_in/2);
out((nc+1):(nc+N_in), (mc+1):(mc+M_in)) = in;

end %function
