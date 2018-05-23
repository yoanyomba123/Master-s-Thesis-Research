
% This code demonstrates the level set evolution (LSE) and bias field estimation
% proposed in the following paper:
%      C. Li, R. Huang, Z. Ding, C. Gatenby, D. N. Metaxas, and J. C. Gore, 
%      "A Level Set Method for Image Segmentation in the Presence of Intensity
%      Inhomogeneities with Application to MRI", IEEE Trans. Image Processing, 2011
%

% Note: 
%    This code implements the two-phase formulation of the model in the above paper.
%    The two-phase formulation uses the signs of a level set function to represent
%    two disjoint regions, and therefore can be used to segment an image into two regions,
%    which are represented by (u>0) and (u<0), where u is the level set function.
%
%    All rights researved by Chunming Li, who formulated the model, designed and 
%    implemented the algorithm in the above paper.
%
% E-mail: lchunming@gmail.com
% URL: http://www.engr.uconn.edu/~cmli/
% Copyright (c) by Chunming Li
% Author: Chunming Li

close all;
clear all;
clc;
%% Add all paths
addpath(genpath('./'));

% compile c files for anisotropic diffusion
compile_c_files;
%% Get Image files
%image = getAllfolderimages('C:/Users/undergrad/Desktop/Publication_Dataset/AMD1/VOLUMEDIR/RawImages/*.tif');
image = getAllfolderimages('C:/Users/undergrad/Desktop/Publication_Dataset/DME6/TIFFs/8bitTIFFs/*.tif');

% remove empty cell arrays
image = image(~cellfun(@isempty, image));
% acquire image set dimension
[mpost, lpost] = size(image);

% extract images and store in img object
datapost = [];
figure;
for item=1:150%numel(image)
    % acquire orignal image for future use
    image_orig = image{1, item};
    
    % store image data into some image struct
    image_struct{1, item} = image_orig;
    
    % filter the image and add filtered image to post_processed
    % data-structure
    data_post_processed{1,item} = imageFilter(image_orig);

end


% Perform anisotropic diffusion
figure;
diffused_image = [];
for item = 1:150%numel(image);
    filtered_image = gather(data_post_processed{1,item});
    JO_EED = CoherenceFilter(filtered_image,struct('T',10,'dt',1,'sigma',3,'rho',1,'Scheme','O','eigenmode',3));
    % sharpen image boundaries through imsharpen
    sharpened_image = imsharpen(JO_EED);
    diffused_image{1, item} = sharpened_image;
end

%% obtaine diffused image
dimensions = [2.51, 3.51, 503.98, 491.98];

IMAGE_VAL = 9;

% crop image
sample_image = imcrop(diffused_image{1, IMAGE_VAL}, dimensions);


% sharpen the image boundaries
sample_image = imsharpen(sample_image, 'Radius', 10);
figure;imagesc(sample_image); colormap gray, colorbar;

% threshold image and convert from grayscale to black and white
BW_image=im2bw(sample_image,graythresh(sample_image));
BW_image = imopen(BW_image, ones(1,4));
BW_image = imfill(BW_image, 'holes');
figure; imshow(BW_image);

% extract connected components from image
CC = bwconncomp(BW_image);

% extract only the largest connected components
numPixels = cellfun(@numel, CC.PixelIdxList);
% sort obtained pixels in descending order
pixel_list = sort(numPixels,'descend');
% extract the 10 biggest connected components
largest_CC = [];
for item = 1:5
    [biggestCC, idx] = max(numPixels);
    largest_CC{1, item} = idx;
    numPixels(1,idx) = 0;
end

%% Remove all smaller connected components from the image
for value = 1:numel(numPixels)
    if(numPixels(value) ~= 0)
       BW_image(CC.PixelIdxList{value}) = 0;
    end
end
figure; imshow(BW_image);

[Gx, Gy] = imgradientxy(BW_image);
[Gmag, Gdir] = imgradient(Gx, Gy);

figure, imshow(Gmag, []), title('Gradient magnitude')
figure, imshow(Gdir, []), title('Gradient direction')
figure, imshow(Gx, []), title('Directional gradient: X axis')
figure, imshow(Gy, []), title('Directional gradient: Y axis')

Gy = imsharpen(Gy, 'radius', 3);
imf =im2bw(Gy, graythresh(Gy));
figure; imshow(imf);
BW = imf;

BW = imfill(BW_image, 'holes');
strel = ones(3,1);
BW = imopen(BW, ones(1, 10));
BW = bwareaopen(BW, 10000);
figure; imshow(BW);
BW2 = BW;
BW2 = imdilate(BW2, strel) & ~BW;
figure; imshow(BW2);
%BW2 = bwareaopen(BW2, 5);
BW2 = imdilate(BW2, ones(5,5));
BW2 = imopen(BW2, ones(3,3));
BW2 = bwareaopen(BW2, 100);
figure; imshow(BW2);

draw = cat(3, BW, BW2, BW_image);
figure;
imagesc(draw); colormap gray;

L1 = bwlabel(BW2);
figure;imagesc(L1);

L1_max = max(L1(:));
figure; imshow(image{1,IMAGE_VAL});
hold on;
color = ['g', 'r', 'y', 'c', 'w','b'];
for i = 1:L1_max
    lt = L1 == i;
    [x,y] = find(lt == 1);
%     if(max(y) < 500) 
%         % compute gradients
%        grad_x = gradient(x);
%        grad_y = gradient(y);
%        %extract gradient at last point
%        size_arr = length(x);
%        x1 = x(size_arr,1);
%        y1 = y(size_arr, 1);
%        x0 = x(size_arr-10, 1);
%        y0 = y(size_arr-10, 1);
%        % compute the slope
%        m = (y0 - y1)/(x0 - x1); 
%        x = x % specifying x limits
%        y = x1+m*(x-x0);
%        figure;
%        plot(y, x,'g');
% %        t = (x == x0) & (y == y0);
% %        indt = find(t);
% %        f_grad = [grad_x(indt), grad_y(indt)]
%     end
    x = medfilt1(x,15);
    y = medfilt1(y,15);
    plot(y,x, color(i), 'lineWidth' , 0.5);
end

hold on;
%% SEGMENT & ACQUIRE INTRA-RETINAL BOUNDARIES
% crop image
im = imcrop(gather(image{1, IMAGE_VAL}), dimensions);
im = imgaussfilt(im, 4);
thresh = graythresh(im);
im_tmp = im2bw(im, thresh);
im_tmp = bwareaopen(im_tmp, 10000);
%figure; imshow(im_tmp);
strel = ones(7, 4);
im_tmp = imdilate(im_tmp,strel) & ~im_tmp;
im_tmp = bwareaopen(im_tmp, 100);
%figure; imshow(im_tmp);

L1 = bwlabel(im_tmp);
L1_max = max(L1(:));
color = ['c', 'w','b','g', 'r', 'y'];
for i = 1:L1_max
    lt = L1 == i;
    [x,y] = find(lt == 1);
%     if(max(y) < 500) 
%         % compute gradients
%        grad_x = gradient(x);
%        grad_y = gradient(y);
%        %extract gradient at last point
%        size_arr = length(x);
%        x1 = x(size_arr,1);
%        y1 = y(size_arr, 1);
%        x0 = x(size_arr-10, 1);
%        y0 = y(size_arr-10, 1);
%        % compute the slope
%        m = (y0 - y1)/(x0 - x1); 
%        x = x % specifying x limits
%        y = x1+m*(x-x0);
%        figure;
%        plot(y, x,'g');
% %        t = (x == x0) & (y == y0);
% %        indt = find(t);
% %        f_grad = [grad_x(indt), grad_y(indt)]
%     end
    x = medfilt1(x,15);
    y = medfilt1(y,15);
    plot(y,x, color(i), 'lineWidth' , 0.5);
end



%%
lt = L1 == 4;
figure; imagesc(lt);
[x,y] = find(lt == 1);


S = regionprops(BW2, 'PixelList');

for i = 1: length(S)
   point_vect = S(i).PixelList;
   for ii = min(point_vect(:, 2)):max(point_vect(:,2))
       
   end
    
end



%%
points = S(1).PixelList;

ii = 1;
for i = min(points(:, 2)):max(points(:,2))
new_points(ii,2) = i;    
new_points(ii,1) = mean(points(points(:,1) == i));
    ii = ii + 1;
end
figure;
plot(new_points(:,1),new_points(:,2))



%%
% extract image boundaries
boundaries = bwboundaries(BW_image);
bwCC = bwconncomp(BW_image);

% find the biggest boundary portion
bsize=cellfun(@length,bwCC.PixelIdxList);
% find maximal connected component
[~,idxLargest]=max(bsize);
% acquire boundary
b1=boundaries{idxLargest};



% plit this curve into two different boundaries
for item = 1:(length(b1)-1)
    slope = (b1(item+1, 1) - b1(item, 1));
    if(slope == 0)
        b1(item, :) = 0;
    end
end

for i = 1:length(b1)
    if(b1(i, 1) == 0 && b1(i, 2) == 0)
        b1(i,:) = [];
    end
end
figure;
plot(b1);
% visualize the boundary
figure;
plot(b1(:,2),b1(:,1),'-c')

%%
% apply polifit to the boundary
p=polyfit(b1(:,2),b1(:,1),3);

[r,c]=find(BW_image);

xx=[0:520];
yy=[];

for i=1:length(xx)
     yy(i)=p(4)+p(3)*xx(i)+p(2)*xx(i)^2+p(1)*xx(i)^3;
      %yy(i) = p(3) + p(2)*xx(i) + p(3)*xx(i)^2;
end
figure;
%imshow(image{1,1}); colorbar;
hold on;
plot(xx,yy,'-m','linewidth',1)
hold on;
plot(b1(:,2),b1(:,1),'-c')



%%
figure;
for i = 1:1%numel(B);
    drawnow;
    [Fx, Fy] = gradient(B{i});
    plot(B{i}(:, 1),B{i}(:, 2));
    hold on;
    plot(Fx, Fy);
end


%%


bounds=bwboundaries(bwCurves);
bwcc=bwconncomp(bwCurves);
% find the biggest boundary portion
bsize=cellfun(@length,bwcc.PixelIdxList);
[~,idxLargest]=max(bsize);
% this was the artifact, get rid of it
bounds(idxLargest)=[];

bsize=cellfun(@length,bounds);
[~,idxLargest]=max(bsize);
b1=bounds{idxLargest};
plot(b1(:,2),b1(:,1),'-c')


p=polyfit(b1(:,2),b1(:,1),3);
figure;imagesc(im);colormap gray;hold on
[r,c]=find(bwCurves);
% plot(c,r,'.r')

xx=[6:504];
yy=[];

v = max(bsize);
for i = 1 :numel(bsize)
    if v == bsize(i,1)
        bsize(i,1) = 0;
        disp('set');
    end
end
[~,idxLargest2]=max(bsize);
b2=bounds{idxLargest2};
plot(b2(:,2),b2(:,1),'-c')


z=polyfit(b2(:,2),b2(:,1),3);


figure;imagesc(im);colormap gray;hold on

[r,c]=find(bwCurves);
% plot(c,r,'.r')

xx=[6:504];
yy=[];

v = max(bsize);
for i = 1 :numel(bsize)
    if v == bsize(i,1)
        bsize(i,1) = 0;
        disp('set');
    end
end
[~,idxLargest3]=max(bsize);
b3=bounds{idxLargest3};
zz=polyfit(b3(:,2),b3(:,1),3);
for i=1:length(xx)
    yyz(i)=z(4)+z(3)*xx(i)+z(2)*xx(i)^2+z(1)*xx(i)^3;
end

for i=1:length(xx)
    yy(i)=p(4)+p(3)*xx(i)+p(2)*xx(i)^2+p(1)*xx(i)^3;
end

for i=1:length(xx)
    yyzz(i)=zz(4)+zz(3)*xx(i)+zz(2)*xx(i)^2+zz(1)*xx(i)^3;
end

plot(xx,yy,'-m','linewidth',1)
plot(xx,yyzz,'-r','linewidth',1)

plot(xx,yyz,'-b','linewidth',1)

%% extract 10 largest boudaries
for item = 1:10
    [largest_boundary, index] = max(boundaries);
    boundaries(index) = 0;
end


figure; imshow(BW_image);
DrawBounds(BW_image,gather(image{1,1}));

%%
figure;

for i = 1:1
    img = gather(datapost{1,i});
    datapost{3, i} = levelSet(gather(datapost{2,1}),image{1,i});
    imagesc(datapost{3,i}); colormap(gray);
    drawnow;
    fr3 = getframe();
end
%
%% subtract bias from segmented image
clear image_analysis;
image_analysis1 = gather(datapost{2,1});
image_analysis = medfilt2(image_analysis1, [10,10]);
figure;
imshowpair(image_analysis1, image_analysis, 'montage');
%%

sample_image=im2bw(image_analysis1,graythresh(image_analysis1));
figure;
imshow(sample_image);
% % acquire old image
% imx = image_analysis;
% % get all bright spots
% bBig = imx==1;
% % remove values less than 200 pixels from image
% bBig = bwareaopen(bBig,500);
% % acquire the minimal value for bBig
% imx(bBig) = min(imx(~bBig));
% 
% % apply a median filter to the image on a 10 x 10 neighborhood
% imx2 = medfilt2(imx,[10,10]);
% % subtract the gaussian filtered image from the original
% imx2 = imx2 - imgaussfilt(imx2,[30,1]);
% 
% % enhance edges of the pre processed image
% %imx2 =  CoherenceFilter(imx2,struct('T',10,'dt',0.1,'rho',0.2,'Scheme','O','eigenmode',1));
% 
% %  imx2 = max(imx2,0);
% imx2 = mat2gray(imx2);
%  
% % threshold and partition the image domain
% imy = imquantize(imx2,multithresh(imx2,4));
% imy2 = imy*0;
% % perform morphological operation on the image 
% for i = 1:max(imy(:))
% bY = imy==i;
% bY = bwareaopen(bY,300);
% bY = imclose(bY,ones(10));
% imy2(bY) = i; 
% end
% %imy = mat2gray(imy);
% %imy2 = mat2gray(imy2);
% %imy2(imy2 <= 0.4) = 0;
% figure;
% imagesc([imx,imx2;imy,imy2]); colorbar; colormap gray;

thresh = multithresh(image_analysis,3);

% apply image threshold
Imgs = imquantize(image_analysis, thresh);

% convert to rgb
RGB = label2rgb(Imgs);

% conver rgb image to gray scale
I = rgb2gray(RGB);
I(I >120 & I < 200) = 0;
figure;
imagesc(I); colorbar;
%% extract connecterd components from the image
image_analysis = I;
thresh = graythresh(image_analysis);
BW = im2bw(image_analysis, thresh);
%CC = bwconncomp(im2bw(imy2));
%numPixels = cellfun(@numel,CC.PixelIdxList);
%[biggest,idx] = min(numPixels);
%BW(CC.PixelIdxList{idx}) = 0;
figure;
imshow(BW);

% bw boundarie time
[B,L] = bwboundaries(BW);
imshowpair(image_analysis, gather(image{1,1}), 'montage');
figure;
imshow(gather(image{1,1}));
hold on;
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5)
end
 
BW2 = bwmorph(BW,'remove');
BW3 = bwmorph(BW2,'skel',Inf);
figure;
imshow(BW3);
%%
%
% get a zeri matrix of size q
grps = 6;
thresh = multithresh(gather(imsharpen(datapost{3,20})), grps);
q = imquantize(gather(datapost{3,20}), thresh);

imy2 = q*0;
% perform morphological operation on the image 
for i = 1:max(imy(:))
bY = imy==i;
bY = bwareaopen(bY,1000);
bY = imclose(bY,ones(10));
imy2(bY) = i; 
end
%imy = mat2gray(imy);
%imy2 = mat2gray(imy2);
imy2(imy2 < 3.5) = 0;
figure;
imagesc([q,imy2]); colorbar; colormap gray;

[B,L] = bwboundaries(imy2);
figure;
imshow(gather(image{1,10}));hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end
%%
L1 = (q*0);

for i = 2:max(q(:))
   temp = q == i;
   % perform morphological closing ( dilate then erode ) -> remove gaps
   temp = imclose(temp,ones(2,7));
   % perform morphological opening ( erode then dilate ) 
   temp = imopen(temp,ones(2,3));
   % remove all connected components that have fewere pixels items 
   temp = bwareaopen(temp,50);
   % label the connected components within an image
   L = bwlabel(temp);
   L1(temp) = L(temp) + max(L(:));   
end
figure; imagesc(L1); colorbar;
hold on 
bound_set = [];
for k = 2:max(L1(:))
    if nnz(L1==k)==0; continue, end 
    B = bwboundaries(L1==k);
    bound_set = [bound_set; B];
    for kk = 1:length(B)
       boundary = B{kk};
       plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
    end
end

row_set = [];

for item = 1: numel(bound_set)
    [row, col] = size(bound_set{item, 1});
    bound_set{item, 2} = row;
    row_set = [row_set row];
end

% sort based on rows;
bs_sort = sort(row_set(:), 'descend');
bs_sort = bs_sort(1:grps);
bs_sort = sort(bs_sort(:), 'ascend');
for i = 1: length(bound_set)
    for j = 1:length(bs_sort)
       if(bs_sort(j, :) == bound_set{i, 2})
          bound_set{i, 3} = bound_set{i, 1}; 
       end
    end
end
%
bou = bound_set(:, 3);
bou = bou(~cellfun(@isempty, bou));

for i = 1: numel(bou)
   x = smooth(bou{i,1}(:,1),'sgolay'); 
   y = smooth(bou{i,1}(:,2),'sgolay');

   bou{i, 2}(:,1) = x;
   bou{i, 2}(:,2) = y;
end
%
figure; imshow(image{1,20}); colorbar;
hold on 
color = ['g', 'r', 'c','m', 'y', 'w', 'b'];
for u = 1:length(bou)
       boundary = bou{u, 2};
       data = polyfit(bou{u, 2}(:,1), bou{u, 2}(:,2), 20);
       f = polyval(data,bou{u, 2}(:,1));
       plot(boundary(:,2), boundary(:,1), color(u), 'LineWidth', 1);
end
% Acquire only array matching the first ith maximal lines

