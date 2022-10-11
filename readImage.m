function [Data,Header,baseFileName] = readImage(acqDirectory,faArray,nimages,imgType)
%readImage Reads Dicom images
%
%     Inputs:
%         acqDirectory - directory where data lives
%         faArray - nominal FA prescribed at the scanner
%         nimages - number of slices
%         imgType - Dixon can have IP, OP, Fat, Water; input not used as
%         the water images were extracted into their own directory.

%Note: All slices from each acquisition (B1, B0, SPGR) should be saved in
%their own folders: example folder directory:
%2DGRE_EPI>FS>FA65
%3DVIBE_VFA_Dixon>FA2
%B0_2DGREMultiSlice>Phase>TE1>FA15
%
%     Outputs:
%         Data: Dicom data
%         Header: Dicom Header
%         baseFileName: used for debug
%
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

Data = cell(nimages,length(faArray));
Header = cell(nimages,length(faArray));

if exist(acqDirectory,'dir') ~= 7
    fprintf('Sorry, directory does not exist. Check your directory path. \n')
end

for ifa = 1:length(faArray)
    
    faDirectory = ['FA',faArray{ifa},'/'];
    if nargin > 3
      faDirectory = [faDirectory,imgType];
    end
    imageDirectory = fullfile(acqDirectory,faDirectory);   
    filePattern = fullfile(imageDirectory, '*.IMA'); % Change to whatever pattern you need.
    filesDirectory = dir(filePattern);
    
    if isempty(filesDirectory) %Nottingham does not append file extension .IMA (uses Unix)
        filePattern = fullfile(imageDirectory, 'IM_*'); % Change to whatever pattern you need.
        filesDirectory = dir(filePattern);
    end
    
    cnt = 0;
    
    for img = 1:1:length(filesDirectory)
        
        baseFileName = filesDirectory(img).name;
        fullFileName = fullfile(imageDirectory, baseFileName);
        fprintf(1, 'Now reading %s\n', baseFileName);
        cnt = cnt + 1;
        Data{cnt,ifa} = double(dicomread(fullFileName));
        %  Data{cnt,ifa} = (dicomread(fullFileName));
        Header{cnt,ifa} = dicominfo(fullFileName);
    end
    
end



end


