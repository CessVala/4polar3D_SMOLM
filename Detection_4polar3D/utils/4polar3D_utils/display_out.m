% =========================================================================
% Function Name: display_out.m
%
% Description:
%   This function provides a file identifier (fid) to direct output to the standard error stream.
%
%   The function:
%     - Is used to ensure compatibility between Octave and MATLAB.
%     - Returns fid = 2, which points to stderr in Octave and displays text in red in MATLAB.
%
% Notes:
%   - fid = 0 would correspond to no output (NULL).
%   - fid = 1 would correspond to stdout (standard output).
%   - fid = 2 is recommended for cross-platform debug/error messages.
%
% Inputs:
%   - None.
%
% Outputs:
%   - fid: File identifier for standard error stream.
%
% Author:
%   Nicolas Bertaux - Institut Fresnel
% 
% Date: December 2017
%
% =========================================================================

% For Octave/MATLAB compatibility
function fid = display_out()

% fid = 0 ; % No output (NULL)
% fid = 1 ; % Standard output (stdout)
fid = 2 ; % Standard error stream for Octave / red text display in MATLAB

end %function
