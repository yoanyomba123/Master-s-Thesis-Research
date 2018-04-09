function res = estimateNoise(octimg, Params)
    octimg = octimg - mean(mean(octimg)); 
    mimg = medfilt2(octimg, Params.ONFLLIN_SEGMENT_NOISEESTIMATE_MEDIAN);
    octimg = mimg - octimg;

    octimg = abs(octimg);

    line = reshape(octimg, numel(octimg), 1);
    res = std(line);
end
