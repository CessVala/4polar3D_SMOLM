% =========================================================================
% Script Name: Merge_tforms.m
%
% Description:
%   This script merges the three transformation files (tforms) generated 
%   by the 'Correct_distortion' script into a single .mat file called 
%   'tformbeads'. These transformations correspond to the corrections 
%   applied between different polarized projections.
%
% Instructions:
%   1. Run this script.
%   2. When prompted, select the three tform files generated from 
%      the 'Correct_distortion' runs.
%   3. Save the merged transformations into a new file (default: tformbeads.mat).
%
% Authors:
%   Miguel Sison              - Institut Fresnel  
%   Cesar Valades-Cruz        - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

%% Merge tform files into a single structure
clear all
close all
clc

%% Prompt user to select the three tform files
[filename, pathname, filterindex] = uigetfile( ...
    { 'tform*.mat', 'Tform MAT-files (*.mat)' }, ...
    'Pick three tform files', ...
    'MultiSelect', 'on');

%% Load each selected tform file
for i = 1:size(filename, 2)
    load([pathname filename{i}])
end

%% Save all loaded tforms into one final file
uisave(who('tform*'), 'tformbeads')
