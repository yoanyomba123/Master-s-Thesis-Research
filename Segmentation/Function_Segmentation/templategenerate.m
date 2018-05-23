function [imgTemplate, energyField, stopNormSeg] = templategenerate(imgRaw, Parameter)
dim1 = Parameter.infoDimension(1);
dim2 = Parameter.infoDimension(2);

truncate = Parameter.truncate;
numPadding1 = Parameter.numPadding1;
numPadding2 = Parameter.numPadding2;
upperBound = truncate(2) * ones(dim1, dim2);
lowerBound = truncate(1) * ones(dim1, dim2);
if (strcmp(Parameter.method, 'simple'))
    
%     mask = Parameter.mask;
    logictemp = (imgRaw >= lowerBound) .* (imgRaw <= upperBound);
%     logictemp = logictemp .* mask;
    tempPure = numPadding1 .* ones(dim1, dim2);
    imgTemplate = imgRaw .* ~logictemp + tempPure .* logictemp;
    energyField = logictemp; % no meaning, just to keep a fluent move
    stopNormSeg = 0;
elseif (strcmp(Parameter.method, 'LSF'))
%     These comments are another way of setting the initial value
%     initialLSF = ones(dim1, dim2);
%     CoorLSFX = [62.8797 61.2595 213.5633 218.4241];
%     CoorLSFY = [196.550632911393,347.234177215190,356.955696202532,175.487341772152];
%     roi = roipoly(ones(dim1, dim2), CoorLSFX, CoorLSFY); 
%     roi = roi + fliplr(roi);
%     initialLSF = initialLSF .* roi .* -2 + initialLSF;
%     Parameter.initialLSF = initialLSF;
    
%     This is the calibration
%     So that parameter will not need to be adjusted

    imgRawBack = imgRaw;
    if(Parameter.booCali)
        imgRaw = polyval(Parameter.p, imgRaw);
    end
    Parameter = segparainitial(Parameter, zeros(size(imgRaw)));
    imgInput = 255 * imadjust(imgRaw, Parameter.bandWidth, [0 1]);
    
%     Do it with automatic stop criteria
%     [energyField, stopNormSeg] = itersegmentationRSFauto(imgInput, Parameter);

    energyField = itersegmentationRSF(imgInput, Parameter); stopNormSeg = 0;
    
    coorROI = Parameter.afterROI;
    logicMuscleMisc = LSFafter(energyField, coorROI);
    logicBone = imgRaw >= upperBound;
    logicAir =  imgRaw <= lowerBound;
    logicMuscle = logicMuscleMisc .* (~logicBone) .* (~logicAir);
    
%     Can be used: dilation
%     koreAverage = fspecial('average', 3);
%     logicMuscle = (imfilter(logicMuscle, koreAverage) > 0);
%     For now, dilation is not used

    logicSoftTissue = ones(dim1, dim2) - logicBone - logicAir - logicMuscle;
    imgTemplate = imgRawBack .* (logicBone + logicAir) + numPadding1 * ones(dim1, dim2) .* logicSoftTissue ...
        + numPadding2 * ones(dim1, dim2) .* logicMuscle;
    
elseif (strcmp(Parameter.method, 'LSFROI'))
    imgInput = 255 * imadjust(imgRaw, [0.015 0.025], [0 1]);
    [energyField, stopNormSeg] = itersegmentationRSFauto(imgInput, Parameter);
    coorROI = Parameter.afterROI;
    preROI = Parameter.preROI;
    logicMuscleMisc = LSFafter(energyField, coorROI);
    logicBone = (imgRaw >= upperBound);
    logicAir =  (imgRaw <= lowerBound);
    logicMuscle = (logicMuscleMisc .* (~logicBone) .* (~logicAir));
    logicSoftTissue = ones(dim1, dim2) - logicBone - logicAir - logicMuscle;
    logicAirBack = logicAir .* (~preROI);
    logicAirTrue = logicAir .* preROI;
    imgTemplate = imgRaw .* (logicAirBack) + numPadding1 * ones(dim1, dim2) .* logicSoftTissue ...
        + numPadding2 * ones(dim1, dim2) .* logicMuscle + (-0.02) .* ones(dim1, dim2) .* logicAirTrue ...
        + 0.04 * ones(dim1, dim2) .* logicBone;

elseif (strcmp(Parameter.method, 'shape'))
    logicBone = imgRaw >= upperBound;
    logicAir =  imgRaw <= lowerBound;
    logicMuscleMisc = zeros(dim1, dim2);
    for i = 1:1:dim1
        for j = 1:1:dim2
            if( ((i-dim1/2)^2 + (j-dim2/2)^2) <= Parameter.radius ^ 2 ); logicMuscleMisc(i,j) = 1; end;
        end
    end
    logicMuscle = logicMuscleMisc .* (~logicBone) .* (~logicAir);
    logicSoftTissue = ones(dim1, dim2) - logicBone - logicAir - logicMuscle;
    imgTemplate = imgRaw .* (logicBone + logicAir) + numPadding1 * ones(dim1, dim2) .* logicSoftTissue ...
        + numPadding2 * ones(dim1, dim2) .* logicMuscle;
    energyField = zeros(512, 512);
    stopNormSeg = 0;
end

end

function logicMuscleMisc = LSFafter(energyField, coorROI)
if( size(coorROI,1) == 2 )
    coorROIX = coorROI(1,:); % [182,55,42,444,453,283];
    coorROIY = coorROI(2,:); % [142,215,363,366,213,135];
    energyField = (roipoly(energyField,coorROIX,coorROIY)) .* energyField;
    logicMuscleMisc = (energyField < 0);
elseif( size(coorROI,1) == 1 )
    energyFieldLogic = zeros(size(energyField,1), size(energyField,2));
    energyFieldLogic(coorROI(2):coorROI(4),coorROI(1):coorROI(3)) = 1;
    energyField = energyField .* energyFieldLogic;
    logicMuscleMisc = (energyField < 0);
elseif( size(coorROI,1) == 0 )
    energyFieldLogic = ones(size(energyField));
    energyField = energyField .* energyFieldLogic;
    logicMuscleMisc = (energyField < 0);
end
    
end 
