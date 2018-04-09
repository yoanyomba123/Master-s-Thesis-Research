% Goal: Adequately Segment OCT Images
% Author: Yoan Yomba

% clear workspace
clc; clear all; close all;

% add all path
addpath(genpath('./'));

% declare variables
FILENAME = 'filenames.list';
FILEPATH = 'paths.txt';
FILTER_OP = 'filter';

% load in images
[images, file, n, paths] = getAllfolderimages('C:/Users/undergrad/Desktop/Publication_Dataset/DME6/TIFFs/8bitTIFFs/*.tif', FILENAME, FILEPATH);

% remove all empty cell arrary contents
images = images(~cellfun('isempty', images));

% Read BScans as a volume 
BScans_UF = Convert_to3d( images);

% Obtained filtered bscans
[Descriptors2.Header,Descriptors2.BScanHeader, slo2, BScans_F] = openOctListHelp(strcat('./',FILENAME), paths{1}, 'metawrite', FILTER_OP);
dimensions = [5.51, 9.51, 501.98, 478.98];

%% Obtain preprocess the unfiltered image
BScans_UF = imagePreprocessing(BScans_UF);
BScans_F = imagePreprocessing(BScans_F);

%% Load Parameters
PARAMETER_FILENAME = 'octseg.param';
med_params = loadParameters('MEDLINELIN', PARAMETER_FILENAME);
onh_params = loadParameters('ONH', PARAMETER_FILENAME);
rpe_params = loadParameters('RPELIN', PARAMETER_FILENAME);
bv_params = loadParameters('BV', PARAMETER_FILENAME);
infl_params = loadParameters('INFL', PARAMETER_FILENAME);
inner_params = loadParameters('INNERLIN', PARAMETER_FILENAME);
onfl_params = loadParameters('ONFLLIN', PARAMETER_FILENAME);


%% Obtain the middle line present between the two largest components in OCT image
% extract the median line between the two highest pixel intensity regions
 medlinevol = segmentMedlineVolume(BScans_UF, med_params);
 
% visualize progress
VisualizeMed(BScans_UF, medlinevol, 'med');
 
% extract the rpe layer present in the image
rpelinevol = segmentRPEVolume(BScans_UF, rpe_params, medlinevol);

% visualize progress
VisualizeMed(BScans_UF, rpelinevol, 'rpe');

% Segment onh line
[onh, onhCenter, onhRadius] = segmentONHVolume(BScans_UF, onh_params, rpelinevol);

% extract infl layer present in the image
infllinevol = segmentINFLVolume(BScans_UF, infl_params, onh, rpelinevol, medlinevol);

% visualize progress
VisualizeMed(BScans_UF, infllinevol, 'infl');

% extract blood vessel regions
bvlinevol = segmentBVVolume(BScans_UF, bv_params,onh, rpelinevol);

% visualize progress
VisualizeMed(BScans_UF, bvlinevol, 'bv');

% segment inner line
[ipl,opl, icl] = SegmentInnnerStack(BScans_UF,rpelinevol, infllinevol, medlinevol,bvlinevol, inner_params)

% segment onfl line
[onfl] = SegmentONflStack(BScans_UF,infllinevol,ipl,icl,opl,bvlinevol,onfl_params);
%%

order = 120;
figure;
for i = 1:size(BScans_UF, 3)
    imshow(BScans_UF(:,:,i)); hold on;
%     plot(medlinevol(i, :), 'r'); hold on;
    plot(rpelinevol(i, :), 'b'); hold on;
    plot(medfilt1(infllinevol(i, :),order), 'c');hold on;
    plot(medfilt1(onfl(i, :),order), 'w'); hold on;
    plot(medfilt1(ipl(i, :),order), 'r'); hold on;
    plot(medfilt1(opl(i, :),order), 'g'); hold on;
    plot(medfilt1(icl(i, :),order), 'y'); hold on;
    drawnow
end
hold off;

%%
I = BScans_UF(:,:,1);
% constrain the image based on upper and lower bounds
infline = infllinevol(1,:); % upper bound
medline = medlinevol(1,:);
rpeline = rpelinevol(1,:); % lower bound
bv = bvlinevol(1,:);
Params = inner_params;
[ipl,opl, icl] = SegmentInner(I,rpeline, infline, medline,bv, Params);
% 1) Normalize intensity values and align the image to the RPE
rpe = round(rpeline);
infl = round(infline);
[alignedBScan, flatRPE, transformLine] = alignAScans(I, inner_params, [rpe; infl]);
flatINFL = infl - transformLine;
medline = round(medline - transformLine);
rpeln = round(rpeline - transformLine);
alignedBScanDSqrt = sqrt(alignedBScan); % double sqrt for denoising

% 3) Find blood vessels for segmentation and energy-smooth 
idxBV = find(extendBloodVessels(bv, Params.INNERLIN_EXTENDBLOODVESSELS_ADDWIDTH, ...
                                    Params.INNERLIN_EXTENDBLOODVESSELS_MULTWIDTHTHRESH, ...
                                    Params.INNERLIN_EXTENDBLOODVESSELS_MULTWIDTH));

averageMask = fspecial('average', Params.INNERLIN_SEGMENT_AVERAGEWIDTH);
alignedBScanDen = imfilter(alignedBScanDSqrt, averageMask, 'symmetric');

idxBVlogic = zeros(1,size(alignedBScan, 2), 'uint8') + 1;
idxBVlogic(idxBV) = idxBVlogic(idxBV) - 1;
idxBVlogic(1) = 1;
idxBVlogic(end) = 1;
idxBVlogicInv = zeros(1,size(alignedBScan, 2), 'uint8') + 1 - idxBVlogic;
alignedBScanWoBV = alignedBScanDen(:, find(idxBVlogic));
alignedBScanInter = alignedBScanDen;
runner = 1:size(alignedBScanDSqrt, 2);
runnerBV = runner(find(idxBVlogic));
for k = 1:size(alignedBScan,1)
    alignedBScanInter(k, :) = interp1(runnerBV, alignedBScanWoBV(k,:), runner, 'linear');
end

alignedBScanDSqrt(:, find(idxBVlogicInv)) = alignedBScanInter(:, find(idxBVlogicInv)) ;
averageMask = fspecial('average', Params.INNERLIN_SEGMENT_AVERAGEWIDTH);
alignedBScanDenAvg = imfilter(alignedBScanDSqrt, averageMask, 'symmetric');


% 4) We try to find the CL boundary.
% This is pretty simple - it lies between the medline and the RPE and has
% rising contrast. It is the uppermost rising border.
extrICLChoice = findRetinaExtrema(alignedBScanDenAvg, Params,2, 'max', ...
                [medline; flatRPE - Params.INNERLIN_SEGMENT_MINDIST_RPE_ICL]);
extrICL = min(extrICLChoice,[], 1);
extrICL(idxBV) = 0;
extrICL = linesweeter(extrICL, Params.INNERLIN_SEGMENT_LINESWEETER_ICL);
flatICL = round(extrICL);
                         
% 5) OPL Boundary: In between the ICL and the INFL
oplInnerBound = flatINFL;

extrOPLChoice = findRetinaExtrema(alignedBScanDenAvg, Params,3, 'min', ...
                [oplInnerBound; flatICL - Params.INNERLIN_SEGMENT_MINDIST_ICL_OPL]);
extrOPL = max(extrOPLChoice,[], 1);
extrOPL(idxBV) = 0;
extrOPL = linesweeter(extrOPL, Params.INNERLIN_SEGMENT_LINESWEETER_OPL);
flatOPL = round(extrOPL);


% 5) IPL Boundary: In between the OPL and the INFL
iplInnerBound = flatINFL;
extrIPLChoice = findRetinaExtrema(alignedBScanDenAvg, Params,2, 'min pos', ...
                [iplInnerBound; flatOPL - Params.INNERLIN_SEGMENT_MINDIST_OPL_IPL]);
extrIPL = extrIPLChoice(2,:);
extrIPL(idxBV) = 0;
extrIPL = linesweeter(extrIPL, Params.INNERLIN_SEGMENT_LINESWEETER_IPL);

figure; imshow(alignedBScanDenAvg); hold on; plot(flatINFL,'r'); hold on; plot(flatICL, 'c'); hold on; plot(extrOPL, 'y');
hold on; plot(extrIPL, 'g');hold on; plot(rpeln, 'b');

icl = extrICL + transformLine;
opl = round(extrOPL + transformLine);
ipl = round(extrIPL + transformLine);

icl(icl < 1) = 1;
opl(opl < 1) = 1;
ipl(ipl < 1) = 1;

ipl(ipl < infl) = infl(ipl < infl);
icl(icl > rpe) = rpe(icl > rpe);
opl(opl > icl) = icl(opl > icl);
opl(opl < ipl) = ipl(opl < ipl);
icl(icl < opl) = opl(icl < opl);

figure; imshow(I);hold on; plot(ipl,'r'); hold on; plot(icl, 'b'); hold on; plot(opl, 'c');
hold on; plot(rpe, 'y'); hold on; plot(infl, 'w');

[onflAuto] = SegmentOnfl(I,infl,ipl,icl, opl,bv,onfl_params)


% find ONFL now
bscanDSqrt = sqrt(I);
% 2) Find blood vessels for segmentation and energy-smooth 
[alignedBScanDSqrt flatICL transformLine] = alignAScans(bscanDSqrt, onfl_params, [icl; round(infl)]);
flatINFL = round(infl - transformLine);
flatIPL = round(ipl - transformLine);

averageMask = fspecial('average', [3 7]);
alignedBScanDen = imfilter(alignedBScanDSqrt, averageMask, 'symmetric');

idxBVlogic = zeros(1,size(alignedBScanDSqrt, 2), 'uint8') + 1;
idxBVlogic(bv(1,:) == 1) = idxBVlogic(bv(1,:) == 1) - 1;
idxBVlogic(1) = 1;
idxBVlogic(end) = 1;
idxBVlogicInv = zeros(1,size(alignedBScanDSqrt, 2), 'uint8') + 1 - idxBVlogic;
alignedBScanWoBV = alignedBScanDen(:, find(idxBVlogic));
alignedBScanInter = alignedBScanDen;
runner = 1:size(alignedBScanDSqrt, 2);
runnerBV = runner(find(idxBVlogic));
for k = 1:size(alignedBScanDSqrt,1)
    alignedBScanInter(k, :) = interp1(runnerBV, alignedBScanWoBV(k,:), runner, 'linear', 0);
end

alignedBScanDSqrt(:, find(idxBVlogicInv)) = alignedBScanInter(:, find(idxBVlogicInv)) ;
% 2) Denoise the image with complex diffusion
noiseStd = estimateNoise(alignedBScanDSqrt, onfl_params);

% Complex diffusion relies on even size. Enlarge the image if needed.
if mod(size(alignedBScanDSqrt,1), 2) == 1 
    alignedBScanDSqrt = alignedBScanDSqrt(1:end-1, :);
end
onfl_params.DENOISEPM_SIGMA = [(noiseStd * onfl_params.ONFLLIN_SEGMENT_DENOISEPM_SIGMAMULT) (pi/1000)];

if mod(size(alignedBScanDSqrt,2), 2) == 1
    temp = alignedBScanDSqrt(:,1);
    alignedBScanDSqrt = alignedBScanDSqrt(:, 2:end);
    alignedBScanDen = real(denoisePM(alignedBScanDSqrt, onfl_params, 'complex'));
    alignedBScanDen = [temp alignedBScanDen];
else
    alignedBScanDen = real(denoisePM(alignedBScanDSqrt, onfl_params, 'complex')); 
end

% Find extrema highest Min, 2 highest min sorted by position
extr2 = findRetinaExtrema(alignedBScanDen, onfl_params, 2, 'min pos th', ...
    [flatINFL + 1; flatIPL - onfl_params.ONFLLIN_SEGMENT_MINDIST_IPL_ONFL]); 
extrMax = findRetinaExtrema(alignedBScanDen, onfl_params, 1, 'min', ...
    [flatINFL + 1; flatIPL - onfl_params.ONFLLIN_SEGMENT_MINDIST_IPL_ONFL]);

% 6) First estimate of the ONFL:
dist = abs(flatIPL - flatINFL);
onfl = extrMax(1,:); 
idx1Miss = find(extr2(1,:) == 0); 
idx2Miss = find(extr2(2,:) == 0); 
onfl(idx2Miss) = flatINFL(idx2Miss); 
onfl(bv(1,:) == 1) = flatIPL(bv(1,:) == 1);
%$onfl(onh(onh == 1)) = flatIPL(onh(onh == 1));
onfl= linesweeter(onfl, onfl_params.ONFLLIN_SEGMENT_LINESWEETER_INIT_INTERPOLATE);

% 7) Do energy smoothing
% Forget about the ONFL estimate and fit a poly trough it. Then Energy!
onfl = linesweeter(onfl, onfl_params.ONFLLIN_SEGMENT_LINESWEETER_INIT_SMOOTH);
onfl = round(onfl);

gaussCompl = fspecial('gaussian', 5 , 1);
smoothedBScan = alignedBScanDen;
smoothedBScan = imfilter(smoothedBScan, gaussCompl, 'symmetric');
smoothedBScan = -smoothedBScan;
smoothedBScan = smoothedBScan ./ (max(max(smoothedBScan)) - min(min(smoothedBScan))) .* 2 - 1;

onfl = energySmooth(smoothedBScan, onfl_params, onfl, find(single(bv(1,:))), [flatINFL; flatIPL]);

% Some additional constraints and a final smoothing
onfl(idx2Miss) = flatINFL(idx2Miss);

onfl = linesweeter(onfl, onfl_params.ONFLLIN_SEGMENT_LINESWEETER_FINAL);

diffNFL = onfl - flatINFL;
onfl(find(diffNFL < 0)) = flatINFL(find(diffNFL < 0));

onflAuto = onfl + transformLine;

figure; imshow(alignedBScanDen); hold on; plot(flatINFL,'r'); hold on; plot(flatICL, 'c'); hold on; plot(extrOPL, 'y');
hold on; plot(extrIPL, 'g');hold on; plot(rpeln, 'b'); hold on; plot(onfl, 'r');


figure; imshow(I);hold on; plot(ipl,'r'); hold on; plot(icl, 'b'); hold on; plot(opl, 'c');
hold on; plot(rpe, 'y'); hold on; plot(infl, 'w');hold on; plot(onflAuto,'g');

%% extract inner lines
[icl, opl, ipl] = segmentInnerLayersVolume(BScans_F, inner_params, onh, rpelinevol, infllinevol, medlinevol, bvlinevol);
 
% visualize progress
VisualizeInnerLines(BScans_UF, icl, opl, ipl);

%% Extract ONFL lines
onfllinevol = segmentONFLVolume(BScans_F(:,:, 1:10), onfl_params, onh, rpelinevol(1:10,:), icl(1:10,:), ipl(1:10,:), infllinevol(1:10,:), bvlinevol(1:10,:));
