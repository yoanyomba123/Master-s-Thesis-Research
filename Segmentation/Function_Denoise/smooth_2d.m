function imgOut = smooth_2d(img, nOrder, nSize)
% input the image and nOrder and nSize

[NamesCoefs, NamesTerms, XPow, YPow, SG] = SavGol(nOrder,nSize);
imgOut = conv2(img, SG(:,:,1), 'same');

end