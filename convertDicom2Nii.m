function convertDicom2Nii(dcmFolder,niiFolder,niiImgName)
%convertDicom2Nii 

%     Inputs:
%         dcmFolder - folder with dicom data
%         niiFolder - folder with nifti data
%         niiImgName - image name for nifti data



% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

fmt = 0; %nii uncompressed, 1 is compressed .nii.gz

cd(niiFolder)
delete *.nii
dicm2nii(dcmFolder, niiFolder, fmt)
S = dir(fullfile('*.nii')); %get the current name
movefile(S.name,[niiImgName,'.nii']) % rename the file
delete *.mat %delete the .mat dile
cd '/Users/gabrielabelsley/OneDrive - Perspectum Ltd/DPhil_SecondProject/DataScriptsResults/3DT1MappingScriptsFunctions/3DT1MapSiemensGit'


end

