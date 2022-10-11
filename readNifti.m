function [img_dcmOrient] = readNifti(imgDirectory,varargin)

%readNifti read .nii images
%
%     Inputs:
%         imgDirectory - directory where data lives
%         varargin - boolean flag: 0=axial; 1= coronal orientation
%     Outputs:
%         img_dcmOrient: data in .mat 
%
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


if nargin>1
flagCoronal = varargin{1,1};
end

img = niftiread(imgDirectory);

if size(size(img),2)== 3
    img_dcmOrient = permute(img,[2 1 3]);
elseif size(size(img),2)== 4
    img_dcmOrient = permute(img,[2 1 3 4]);
elseif size(size(img),2)== 2
    img_dcmOrient = permute(img,[2 1]);
end

if nargin>1 && flagCoronal==1
    img_dcmOrient = double(flip(flip(flip(img_dcmOrient,1),2),3));

else %default = axial
img_dcmOrient = double(flip(flip(img_dcmOrient,1),3));
end
%flip(,1) flips along the columns
%flip(,2) flips along the rows

%Changes of nii relative to dcm:
% nii permutes x and y 
% nii reverses the z direction axis

end

