% %% Write image raw


% adding specific folders to current path
addpath(genpath('.'));
addpath(genpath('D:\Yoan\Git\leverUtilities\src\MATLAB'));
addpath(genpath('D:\Yoan\Git\leverjs'));
addpath('D:\Yoan\Git\leverjs\matlab\');
Normal = BScans_stack{1, 3};
% % write to h5 lever format in specified folder
curr_path = 'D:\Yoan\DATA\NORMAL\';
vol_name = 'NORMALRAW4';
volumeToH5(Normal, char(vol_name), char(curr_path));
Import.leverImport(curr_path, char(curr_path), char(vol_name), '');
% 
% 
% 
strDB ='D:\Yoan\DATA\NORMAL\NORMALRAW4.LEVER'
AddSQLiteToPath();
conn = database(strDB, '','', 'org.sqlite.JDBC', 'jdbc:sqlite:');
CONSTANTS=Read.getConstants(conn);
% % sigma = 7;
% type = ""
% volume = dog3D(volume, sigma, type);
Normal(1:100, :,:) = 0; 
Normal(400:end, :,:) = 0;

tic
sigma=10;
% First Gaussian Operation
imd1 = imgaussfilt3(Normal,sigma/sqrt(2));
imd2 = imgaussfilt3(Normal,sigma*sqrt(2));
imdog = imd1 - imd2;

% find image borders

lm=multithresh(imdog,3);
q=imquantize(imdog,lm);
imp=max(imdog,[],3);
q(1:100, :,:) = min(q(:)); 
q(300:end, :,:) = min(q(:));
imp=max(q,[],3);figure;imagesc(imp);

for i = 400:size(imp,2)
   imp(i, :) = 0;
end

figure;imagesc(imp); colorbar
bw=logical(q>=3);

[faces,verts]=isosurface(bw,eps);
norms = isonormals(bw,verts);
edges=Segment.MakeEdges(faces);

% make edges and faces zero indexed
edges=edges-1;
faces=faces-1;

% subtract extract for padding
maxRad = verts-repmat(mean(verts,1),size(verts,1),1);
maxRad=sum(abs(maxRad),2);
maxRad=max(maxRad);
newCell=[];
newCell.time=1;
newCell.centroid=[mean(verts,1)]-1;
newCell.edges=edges;
newCell.faces=faces;
newCell.verts=verts;
newCell.normals=norms;
xyzPts=[];

idxPixels=find(bw);
% note that xyzPts are on the correct range, as idxPixels is unpadded...
[xyzPts(:,2),xyzPts(:,1),xyzPts(:,3)]=ind2sub(size(bw),idxPixels);

newCell.pts=uint16(xyzPts);
newCell.maxRadius=maxRad;
newCell.channel=1;

Write.CreateCells_3D(conn,newCell)
toc

%%
clear;
close all
% adding specific folders to current path
addpath(genpath('.'));
addpath(genpath('D:\Yoan\Git\leverUtilities\src\MATLAB'));
addpath(genpath('D:\Yoan\Git\leverjs'));
addpath('D:\Yoan\Git\leverjs\matlab\');
load  imagedb.mat 

specified_path = 'D:\Yoan\IMAGING';
%%
for i = 1 : size(BScans_stack, 2)   
    for j = 1: size(BScans_stack, 1)
        Normal = BScans_stack{j, i};
        if(i == 1)
            path_ext = "\AMD\";
            vol_name = strcat("AMD", num2str(j));
        end
        if(i == 2)
            path_ext = "\DME\";
            vol_name = strcat("DME", num2str(j));
        end
        if(i == 3)
            path_ext = "\NORMAL\";
            vol_name = strcat("NORMAL", num2str(j));
        end
        curr_path = strcat(specified_path, path_ext);
        
        

        % mask the image
        Normal(1:50, :,:) = 0; 
        Normal(450:end, :,:) = 0;
        Normal(:,1:40, :) = 0;
        Normal(:,end-40:end, :) = 0;
        volume = Normal;

        sigma=10;
        % First Gaussian Operation
        imd1 = imgaussfilt3(Normal,sigma/sqrt(2));
        imd2 = imgaussfilt3(Normal,sigma*sqrt(2));
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
        Normal(:,:,:,1) = Normal;
        Normal(:,:,:,2) = bw;
        Normal(:,:,:,3) = imdog;

% Cannot write to '/Images/Original' because the edge of the dataset would be exceeded, and '/Images/Original' is not extendable along all affected extents.

        % write to h5 lever format in specified folder
        volumeToH5(Normal, char(vol_name), char(curr_path));
        Import.leverImport(curr_path, char(curr_path), char(vol_name), '');
        
        
        [faces,verts]=isosurface(bw,eps);
        norms = isonormals(bw,verts);
        edges=Segment.MakeEdges(faces);

        % make edges and faces zero indexed
        edges=edges-1;
        faces=faces-1;

        % subtract extract for padding
        maxRad = verts-repmat(mean(verts,1),size(verts,1),1);
        maxRad=sum(abs(maxRad),2);
        maxRad=max(maxRad);
        newCell=[];
        newCell.time=1;
        newCell.centroid=[mean(verts,1)]-1;
        newCell.edges=edges;
        newCell.faces=faces;
        newCell.verts=verts;
        newCell.normals=norms;
        xyzPts=[];

        idxPixels=find(bw);
        % note that xyzPts are on the correct range, as idxPixels is unpadded...
        [xyzPts(:,2),xyzPts(:,1),xyzPts(:,3)]=ind2sub(size(bw),idxPixels);

        newCell.pts=uint16(xyzPts);
        newCell.maxRadius=maxRad;
        newCell.channel=1;

        conn = database(curr_path + vol_name +'.LEVER', '','', 'org.sqlite.JDBC', 'jdbc:sqlite:');
        AddSQLiteToPath();
        Write.CreateCells_3D(conn,newCell)
    end
end