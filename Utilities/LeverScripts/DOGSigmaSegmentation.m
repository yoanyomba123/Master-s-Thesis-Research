function result = DOGSigmaSegmentation(BScans,sigma)
%DOGSIGMASEGMENTATION Summary of this function goes here
%   Detailed explanation goes here
    % First Gaussian Operation
    imd1 = imgaussfilt3(BScans,sigma/sqrt(2));
    imd2 = imgaussfilt3(BScans,sigma*sqrt(2));
    imdog = imd1 - imd2;
    % find image borders

    lm=multithresh(imdog,3);
    q=imquantize(imdog,lm);
    imp=max(imdog,[],3);
    q(1:100, :,:) = min(q(:)); 
    q(300:end, :,:) = min(q(:));
    imp=max(q,[],3);
    bw=logical(q>=3);
    CC = bwconncomp(bw);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    % find indexes of all but 4 largest connected components
    [biggest,idx] = mink(numPixels, numel(numPixels)-2);
    for value = 1: numel(idx)
       bw(CC.PixelIdxList{idx(value)}) = 0; 
    end
    result = bw;
end

