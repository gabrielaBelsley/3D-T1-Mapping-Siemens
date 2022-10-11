function [imgMatchResCoordPosPatient] = imageXYZPosInterpSliceAxial(img1,imgHeaderLowRes,imgHeaderHighRes,nSlicesImg2,mapTitleLowRes,mapTitleHighRes,varargin)

%imageXYZPosInterpSliceAxial Interpolation in 3 dimensions for axial slices

    % img1 = Low Resolution Image
    % img2 = High Resolution Image
%     Input: 
%         Img1: normally low resolution image (but the function can also downsample from a high to a low resolution)
%         imgHeaderLowRes: Dicom Header of low resolution image
%         imgHeaderHighRes: Dicom Header of high resolution image
%         nSlicesImg2: number of slices for high resolution image
%         mapTitleLowRes: Name of acquisition protocol for image 1
%         mapTitleHighRes: Name of acquisition protocol for image 2
%         varargin: 
%         if nargin == 6 : direction of acquisition foot-head:
%         flagFootHead = 1 then acquired from foot-head - this is the case for
%         Philips - but not for Siemens OCMR
%         if nargin == 7 : change dicom slice location of image 1 to be
%         equal to the slice location of image 2, used when interpolating
%         GRE-EPI to VFA as the slice-locations are very close but not
%         exactly the same. See consequences of not using this varargin in DPhil_SecondProject/CodePipeline/InterpolationB12VFARes.docx
%         if nargin == 8 : specify image slice location before interpolation

%     Output:
%            imgMatchResCoordPosPatient: interpolated image

 % Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021
   

    fprintf(1, '\n'); %Blank Line
    fprintf (['Interpolating Low Resolution ',mapTitleLowRes,'to high resolution',mapTitleHighRes])
    fprintf(1, '\n'); %Blank Line
    
    % Only Inteprolates Coordinates that are NOT NaN: commented out because
    % Linear can handle NaN: it just places those pixels to NaN
    %CoordDiffNaN     = ~isnan(img1);
    %     valid     = ~isnan(ImgBefInterp);
    % ImgBefInterpValid = ImgBefInterp(valid);
 
    flagFootHead = 0; %Acquired from head>foot direction (default OCMR-Prisma)
    if nargin>6
        flagFootHead = varargin{1,1}; %Acquired from foot>head direction  (Philips)
    end
    if nargin>7
        PlaceEPIZEqual2VFAZ= varargin{1,2}; 
    end
    if nargin>8
        ZLocImgBefInterp= varargin{1,3}; 
    end
    
    %LowRes/Original image before interpolation 
    firstPixelCenterCoord = imgHeaderLowRes.ImagePositionPatient;
    pixelSpacing = double(imgHeaderLowRes.PixelSpacing);
    if isfield(imgHeaderLowRes,'SpacingBetweenSlices')
        % 2D MultiSlice acquisitions have SpacingBetweenSlices = slice thickness+gap
        pixelSpacing(3,1) = double(imgHeaderLowRes.SpacingBetweenSlices);
    else
        % 3D have no SpacingBetweenSlices, thus we just use SliceThickness
        pixelSpacing(3,1) = double(imgHeaderLowRes.SliceThickness);
    end
    nRows = double(imgHeaderLowRes.Rows);
    nColumns = double(imgHeaderLowRes.Columns);
    nSlices = size(img1,3);
    
    %Center coordinates of each pixel according to the scanner coordinate system
    xOriginal = (firstPixelCenterCoord(1):pixelSpacing(1):(nColumns-1)*pixelSpacing(1)+firstPixelCenterCoord(1));
    yOriginal = (firstPixelCenterCoord(2):pixelSpacing(2):(nRows-1)*pixelSpacing(2)+firstPixelCenterCoord(2));
    if flagFootHead == 0
        zOriginal = (firstPixelCenterCoord(3):-pixelSpacing(3):-(nSlices-1)*pixelSpacing(3)+firstPixelCenterCoord(3));
    elseif flagFootHead == 1
        zOriginal = (firstPixelCenterCoord(3):pixelSpacing(3):(nSlices-1)*pixelSpacing(3)+firstPixelCenterCoord(3));
    end
    
    
    %High Resolution image
    firstPixelCenterCoordImg2 = imgHeaderHighRes.ImagePositionPatient;
    pixelSpacingImg2 = double(imgHeaderHighRes.PixelSpacing);
    
    if isfield(imgHeaderHighRes,'SpacingBetweenSlices')
    pixelSpacingImg2(3,1) = double(imgHeaderHighRes.SpacingBetweenSlices);
    else
        pixelSpacingImg2(3,1) = double(imgHeaderHighRes.SliceThickness);
    end
    nRowsImg2 = double(imgHeaderHighRes.Rows);
    nColumnsImg2 = double(imgHeaderHighRes.Columns);
    
    %Center coordinates of each pixel according to the scanner coordinate system
    %x = columns
    %y = rows
    %Matlab convention: (r,c) = (y,x)
    xImg2 = (firstPixelCenterCoordImg2(1):pixelSpacingImg2(1):(nColumnsImg2-1)*pixelSpacingImg2(1)+firstPixelCenterCoordImg2(1));
    yImg2 = (firstPixelCenterCoordImg2(2):pixelSpacingImg2(2):(nRowsImg2-1)*pixelSpacingImg2(2)+firstPixelCenterCoordImg2(2));
    if flagFootHead == 0
        zImg2 = (firstPixelCenterCoordImg2(3):-pixelSpacingImg2(3):-(nSlicesImg2-1)*pixelSpacingImg2(3)+firstPixelCenterCoordImg2(3));
    elseif flagFootHead == 1
        zImg2 = (firstPixelCenterCoordImg2(3):pixelSpacingImg2(3):(nSlicesImg2-1)*pixelSpacingImg2(3)+firstPixelCenterCoordImg2(3));      
    end
    [xImg2Grid,yImg2Grid,zImg2Grid] = meshgrid(xImg2,yImg2,zImg2);
    
    if exist('PlaceEPIZEqual2VFAZ','var') && ~isempty(PlaceEPIZEqual2VFAZ)
       zOriginal = zImg2(PlaceEPIZEqual2VFAZ(:,1));
    end
    
    if exist('ZLocImgBefInterp','var')%B0 with more slices than greEPI
        zOriginal = ZLocImgBefInterp.';
    end

    [xOriginalGrid,yOriginalGrid,zOriginalGrid] = meshgrid(xOriginal,yOriginal,zOriginal);
    
    
    imgMatchResCoordPosPatient = interp3(xOriginalGrid,yOriginalGrid,zOriginalGrid,img1,xImg2Grid,yImg2Grid,zImg2Grid,'linear',NaN);
% use linear and not cubic or makima:
% Reason: in locations where the B1+ was set to 1e6 the makima and cubic
% give wrong interpolation values (they use the 1e6 in the calculation of the interpolation)
% Linear also has the advantage of working with NaN, thus we do not need to
% set 1e6 in B1MapScript, makima and cubic do not do interpolation when
% ExtrapValue in B1MapScript.m is NaN, thus the reason we had to place a
% large number so it would still do the interpolation and then we just
% ignored regions in the interpolated region with very high values. The
% problem is it was corrupting the B1+ correction values inside the phantom
% when there was a region inside the phantom placed to 1e6. Thus the
% neighbouring pixels where the B1+Corr was done correctly were corrupted
% with the 1e6 during the interpolation and giving wrong values in the
% interpolated image. Linear simply places the interpolated value to NaN
% when the original value was NaN

% LR = [5 10 12.5 NaN]
% XLR = 1:4;
% valid = ~isnan(LR);
% XHR = 1:0.5:4;
% XLR(valid)
% LR(valid)
% HR=interp1(XLR(valid),LR(valid),XHR,'linear')
% Linear gives the right values, but makima not
    
end

