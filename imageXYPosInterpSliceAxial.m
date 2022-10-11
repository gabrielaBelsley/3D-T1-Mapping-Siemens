function [imgMatchResCoordPosPatient] = imageXYPosInterpSliceAxial(img1,imgHeaderLowRes,imgHeaderHighRes,mapTitleLowRes,mapTitleHighRes)
    
%imageXYPosInterpSliceAxial Interpolation in 2 dimensions (in-plane) for matching
    %axial slice locations
    
    %     Input: 
%         Img1: normally low resolution image (but the function can also downsample from a high to a low resolution)
%         image we want to interpolate respecting the coordinates from
%         ImagePostion Patient. Interpolate to the resolution of image 2 (equal number of rows, columns and pixel spacing);
%         imgHeaderLowRes: Dicom Header of low resolution image
%         imgHeaderHighRes: Dicom Header of high resolution image
%         mapTitleLowRes: Name of acquisition protocol for image 1
%         mapTitleHighRes: Name of acquisition protocol for image 2

%     Output:
%            imgMatchResCoordPosPatient: interpolated in-plane image
    
    % Example: image 1 = B0 Map, image 2 = B1 Map
    %or % Image 1 = Siemens B1 & image 2 = 3D VIBE
    
  % Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021
  
  % DOUBLE CHECK that SliceLocation between the 2 images is ~ equal:
    fprintf(1, '\n'); %Blank Line
    fprintf(['IMPORTANT: Check the [',imgHeaderLowRes.ProtocolName,'] and the [',imgHeaderHighRes.ProtocolName,'] seq are at interp similar SliceLocations: \n'])
    fprintf([mapTitleLowRes,': Slice Location: %0.4f \n'], imgHeaderLowRes.ImagePositionPatient(3)) 
    fprintf([mapTitleHighRes,': Slice Location: %0.4f \n'], imgHeaderHighRes.ImagePositionPatient(3)) 
    fprintf(1, '\n'); %Blank Line
    
    
    %Original image before interpolation through matching image position
    firstPixelCenterCoord = imgHeaderLowRes.ImagePositionPatient;
    pixelSpacing = double(imgHeaderLowRes.PixelSpacing);
    nRows = double(imgHeaderLowRes.Rows);
    nColumns = double(imgHeaderLowRes.Columns);
    
    %Center coordinates of each pixel according to the scanner coordinate system
    xOriginal = (firstPixelCenterCoord(1):pixelSpacing(1):(nColumns-1)*pixelSpacing(1)+firstPixelCenterCoord(1));
    yOriginal = (firstPixelCenterCoord(2):pixelSpacing(2):(nRows-1)*pixelSpacing(2)+firstPixelCenterCoord(2));
    [xOriginalGrid,yOriginalGrid] = meshgrid(xOriginal,yOriginal);
    
    
    firstPixelCenterCoordImg2 = imgHeaderHighRes.ImagePositionPatient;
    pixelSpacingImg2 = double(imgHeaderHighRes.PixelSpacing);
    nRowsImg2 = double(imgHeaderHighRes.Rows);
    nColumnsImg2 = double(imgHeaderHighRes.Columns);
    
    %Center coordinates of each pixel according to the scanner coordinate system
    %x = columns
    %y = rows
    %Matlab convention: (r,c) = (y,x)
    xImg2 = (firstPixelCenterCoordImg2(1):pixelSpacingImg2(1):(nColumnsImg2-1)*pixelSpacingImg2(1)+firstPixelCenterCoordImg2(1));
    yImg2 = (firstPixelCenterCoordImg2(2):pixelSpacingImg2(2):(nRowsImg2-1)*pixelSpacingImg2(2)+firstPixelCenterCoordImg2(2));
    [xImg2Grid,yImg2Grid] = meshgrid(xImg2,yImg2);
    
    imgMatchResCoordPosPatient = interp2(xOriginalGrid,yOriginalGrid,img1,xImg2Grid,yImg2Grid,'linear',NaN);
    
end

