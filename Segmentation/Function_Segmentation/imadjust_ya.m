function imgOut = imadjust_ya(img, bandWidth, booNormalize)

    if(nargin == 2)
        booNormalize = 1;
    end

    % Trying to adjust the brightness value range of specific image.
    % Result will be a image in uint8 (0 - 255)
    if(booNormalize)
        img = normalize01(img);
    end

    img = imadjust(img, [bandWidth], [0 1]);
    A = 255;
    imgOut = img .* A;

end

function imgOut = normalize01(img)

    % Normalize to the range of [0,1]
    imgmin = min(img(:));
    imgmax = max(img(:));
    imgOut = (img-imgmin) / (imgmax-imgmin);

end