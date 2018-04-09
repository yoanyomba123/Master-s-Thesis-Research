function [imagesfile, fileID, n, paths] = getAllfolderimages(x, file, file_path)
% Returns all images files specified at x in image struct
imagefilesPatient1 = dir(x);      
nfiles = length(imagefilesPatient1);    % Number of files found
images = [];
x = nfiles;

% open file;
fileID = fopen(file, 'w');
filePaths = fopen(file_path, 'w');
for ii=1:nfiles
   currentfilename = imagefilesPatient1(ii).name;
   currentfolder = imagefilesPatient1(ii).folder;
   sample = strcat(currentfolder, '\');
   paths{ii} = sample;
   fprintf(fileID, '%s\n',currentfilename); 
   fprintf(filePaths, '%s\n',currentfolder); 
   fullFileName = fullfile(imagefilesPatient1(ii).folder, currentfilename);
   %currentimage = gpuArray(imcrop(imread(fullFileName), [2.51, 134.51, 509.98, 238.98]));
   currentimage = gpuArray(imread(fullFileName));

   I = mat2gray(currentimage);
   imageset{ii} = I;

end
imagesfile = imageset;
fclose(fileID);
fclose(filePaths);
n = x;
end