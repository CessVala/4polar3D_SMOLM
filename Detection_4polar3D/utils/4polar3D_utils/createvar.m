function createvar(num, data1)
% =========================================================================
% Function: createvar
%
% Description:
%   Dynamically creates a variable named "data<num>" in the base workspace,
%   where <num> is the numeric input. This is useful for programmatically
%   assigning labeled variables when working interactively in MATLAB.
%
% Inputs:
%   num    - (numeric) A number used to construct the variable name.
%   data1  - (any type) The value to assign to the created variable.
%
% Output:
%   A new variable named data<num> appears in the base workspace.
%
% Example usage:
%   createvar(3, rand(10));  % Creates variable 'data3' in base workspace
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel  
%   Miguel Sison              - Institut Fresnel  
%   Cesar Valades-Cruz        - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

    assignin('base', ['data' num2str(num)], data1);

end
