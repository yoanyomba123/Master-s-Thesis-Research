% This is the main program.
% On segmentation of data set provided to my by Maomao.
% This program utilizes the level set method (multi-phase).

% Provided by Pengwei Wu, under the supervision of Tianye Niu.

clear; close all; clc
addpath(genpath('../'));
addpath(path, '../Function_Misc/');
addpath(path, '../Function_Denoise/');
addpath(path, '../Function_Segmentation/');

booDebug = 1;

%% Read the raw data

pathBase = '..\Data\'; % Do not forget the final back slash.
pathName = 'Demo.bin';
imgRaw = read_raw_data(strcat(pathBase, pathName), [512 512]);

if(booDebug); figure, imshow(imgRaw, [0.015 0.025]); title('Original image'); end

imgUnSeg = imadjust_ya(imgRaw, [0.013 0.020], 0);
if(booDebug); figure, imshow(imgUnSeg, [0 255]); title('Unsegmented image (after normalization)'); end

imgUnSegSmooth = smooth_2d(imgUnSeg, 2, 5);
if(booDebug); figure, imshow(imgUnSegSmooth, [0 255]); title('Unsegmented image (after normalization and smoothing'); end

ParaSeg = struct('iterOuter', 100, 'sigma', 3,  ...
    'timeStep', 0.1, 'muBase', 0.1, 'nuBase', 0.01, 'epsilon', 1, ...
    'scale', 0.5, 'save', []);
tic
[M1, M2, M3] = imseg_levelset(imgUnSegSmooth, ParaSeg, 1);
toc