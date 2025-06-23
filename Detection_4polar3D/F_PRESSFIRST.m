% =========================================================================
% Script Name: PRESSFIRST.m
%
% Description:
%   This script sets the working directory and automatically adds all 
%   subdirectories to the MATLAB search path. It is an optional setup 
%   step that can also be done manually by the user if preferred.
%
%   It scans the current folder, identifies all directories, and 
%   recursively adds them to the MATLAB path using `genpath`.
%
% Instructions:
%   1. Place this script in the root processing folder.
%   2. Run the script to set up the environment.
%   3. You should see 'Set path finished!!!' when it completes.
%
% Notes:
%   This step simplifies access to all functions in subfolders.

%
% Authors:

%   Miguel Sison              - Institut Fresnel  
%   Cesar Valades-Cruz        - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

%% Initialize
clear all
close all
clc

%% List all directories in the current folder
listing = dir;

%% Keep only directories
val = [listing.isdir];
listing = listing(val);
mat = struct2cell(listing);
matt = find(cell2mat(mat(4, :)));

%% Add all subdirectories to MATLAB path
for i = 3:size(mat, 2)  % Skip '.' and '..'
    stringss = char(mat(1, i));
    addpath(genpath([pwd '/' stringss]));
end

%% Confirm completion
disp('Set path finished!!!');