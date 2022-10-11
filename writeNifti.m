function writeNifti(img2Save2Nii,niiImgName,niiFolder,niiHeaderFolder,niiImgHeaderName,flagCoronal,varargin)
%writeNifti write images to .nii

% Inputs:
    %img2Save2Nii: image to write to .nii
    %niiImgName: nifti image name
    %niiFolder: folder to save nifti image
    %niiHeaderFolder: folder directory with image to extract nifti header
    %niiImgHeaderName: name of nifti image to extract nifti header
    %flagCoronal: 1 if coronal; 0 if axial
    %varargin: directory where the .dcm image lives to first convert this image to
    %.nii and extract header to subsequently use to write the new nifti image

    
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

%when the header to write the .nii comes from an image that is still in
%.dcm format and first needs to be converted to .nii to use niftiinfo
if (~isempty(varargin))
    dcmFolderHeader = varargin{1,1};
    convertDicom2Nii(dcmFolderHeader,niiHeaderFolder,niiImgHeaderName)
end
header_nii = niftiinfo([niiHeaderFolder,niiImgHeaderName,'.nii']);
%--------write to .nii the masked magnitude TE 1
img2Save2NiiFlipped = (flip(flip(permute(img2Save2Nii,[2 1 3]),2),3));
if flagCoronal == 1
    img2Save2NiiFlipped = flip(flip(flip(permute(img2Save2Nii,[2 1 3]),1),2),3);
end
header_nii.BitsPerPixel = 64;
header_nii.Datatype = 'double';
niftiwrite(img2Save2NiiFlipped,[niiFolder,niiImgName],header_nii);

end

