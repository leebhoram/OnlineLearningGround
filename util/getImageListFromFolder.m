function [ ImageSequence] = getImageListFromFolder( folderName )

ImageSequence = [];

if nargin == 0    
    error('No folder specified');
end

olderFolder = cd(folderName);
ImageFiles = dir('*.png');
   
if isempty(ImageFiles)
    disp('No image is found in the folder.');
    cd(olderFolder);
    return;
end

numfiles = length(ImageFiles);
ImageSequence = cell(numfiles,1);

if numfiles > 0
    disp(folderName);
else
    cd(olderFolder);
    return;
end

if strcmp(folderName(end),'/')
    for n=1:numfiles    
        ImageSequence{n} = strcat(folderName,ImageFiles(n).name);
    end
else
    for n=1:numfiles
        ImageSequence{n} = strcat(folderName,'/', ImageFiles(n).name);
    end
end

cd(olderFolder);
