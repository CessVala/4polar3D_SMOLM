% =========================================================================
% Script Name: STHERM_param.m
%
% Description:
%   This script sets global variables used for 4polar3D STHERM analysis.
%   It defines user-configurable parameters such as window size and 
%   Gaussian width, as well as internal indexing constants for data structures.
%
%   Only modify values in the "USER PARAMETERS" section.
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel  
%   Miguel Sison              - Institut Fresnel  
%   Cesar Valades-Cruz        - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

%% =========================================================================
%% DO NOT MODIFY THE FOLLOWING LINES
%% =========================================================================

global wn r0 pfa 
global NumImg NumStart NumImgVector
global DEBUG DISPLAY 
global char_slash

%% =========================================================================
%% USER PARAMETERS — CAN BE MODIFIED
%% =========================================================================

% Uncomment depending on your operating system
% Windows:
% char_slash = '\' ;
% Unix/Linux/MacOS:
char_slash = '/' ;

% Window size (in pixels) for testing
% Must be an odd number
wn = 11;

% Half-width of Gaussian (in pixels)
% 1.3 for EMCCD, 2.6 for sCMOS
% Example: TL = 1.0x --> r0 = 1.3
%          TL = 1.5x --> r0 = 1.3 * 1.5 = 1.95
r0 = 1.3;

% False alarm probability threshold
% 28 => 1E-7, 24 => 1E-6
pfa = 28;

% Total number of images to analyze
NumImg = 8000;

% Starting frame index
NumStart = 1;

% Number of images used to estimate the vector
NumImgVector = 500 + NumStart;

% Debug mode: 0 = off, 1 = on
DEBUG = 0;

% Display mode: 0 = no visualization, 1 = show detection/localization output
DISPLAY = 1;

%% =========================================================================
%% DO NOT MODIFY THE FOLLOWING LINES
%% =========================================================================

global param_p_num param_p_i param_p_j param_p_alpha param_p_sig2 param_p_r0 param_p_ok
global param_t_num param_t_i param_t_j param_t_r0 param_t_sig2  param_t_tps
global nb_param_t 

% Row indices for matrices:
% liste_part, compat, detect_reconnex
param_p_num         = 1;
param_p_i           = 2;
param_p_j           = 3;
param_p_alpha       = 4;
param_p_sig2        = 5;  % noise variance
param_p_r0          = 6;
param_p_ok          = 7;

% liste_traj
nb_param_t = 8;
param_t_num         = 1;
param_t_i           = 2;  % trajectory i-position
param_t_j           = 3;  % trajectory j-position
param_t_r0          = 4;  % PSF width
param_t_sig2        = 5;  % ij variance
param_t_alpha       = 6;  % intensity
param_t_sig2_b      = 7;  % noise variance (for SNR)
param_t_tps         = 8;  % time index (used for drift correction)
