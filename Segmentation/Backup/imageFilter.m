function imagesfile = imageFilter(I)
% Returns all images files specified at x in image struct
% Performs image filtration

% convert to regular array type
I = gather(I);
% implement level set technique and extract image bias
% I_ls = levelSet(I, I);
% 
% % releave image from bias
% I = I .* I_ls;

% threshold the image
thresh = multithresh(gather(I),15);

% apply image threshold
Imgs = imquantize(I, thresh);

% convert to rgb
RGB = label2rgb(Imgs);

% conver rgb image to gray scale
I = rgb2gray(RGB);
% apply michelle filter
I2 = medfilt2(I, [3,3]) - imgaussfilt(I, 256); 


% Imgo = I > 0.0001;
% 
% 
% Img0 = Imgs .* Imgo;
% % removes redundant noise but also removes edge intensities
% I = (im2double(I) .* Img0) .* Imgs;
% 
% figure;
% imshow(I); colormap gray;


% apply michelle filter
I2 = im2double(I2) - stdfilt(im2double(I2));
% return image
imagesfile = gpuArray(I2);

end