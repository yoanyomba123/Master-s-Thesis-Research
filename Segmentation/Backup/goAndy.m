%im=imread('./images/OCT-Images/042.tif');
im =I;
% mask off blank region
imtexture=entropyfilt(im);
bwTexture=im2bw(imtexture,graythresh(imtexture));
bwTexture=bwmorph(bwTexture,'erode',5);

% imf=medfilt2(im)-imgaussfilt(im,[256,256]);
h=fspecial('log',100,3);
imlog=imfilter(im,h);

bwLog=im2bw(imlog,graythresh(imlog));
figure; imshow(bwLog);
% bwCurves=bwmorph(bwLog,'skel',Inf);
bwCurves=bwLog;

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


% Extract centroids
%% distance from  a cc to p
targetCC=10;
[xCC,yCC]=ind2sub(size(im),bwcc.PixelIdxList{targetCC});
pY=[];
for i=1:length(xCC)
    pY(i)=p(4)+p(3)*xCC(i)+p(2)*xCC(i)^2+p(1)*xCC(i)^3;
end
dy=(pY-yCC')
min(dy)

plot(xx,yy,'-m','linewidth',2)
