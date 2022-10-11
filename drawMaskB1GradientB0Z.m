function [maskForB1Map_gradientB0Z,Coord_maskForB1Map_gradientB0Z,maskForB1Map_gradientB0Z_VFARes,maskForB0Map] = drawMaskB1GradientB0Z(flagCoronal,gre2DEPIHeader,nslices_greEPI,nr_greEPI,nc_greEPI,faHeader,faData,nslices_VFA,nFA_VFA,nr_VFA,nc_VFA,B0Header,nr_B0,nc_B0,nslices_B0)

%DRAWMASKB1GRADIENTB0Z draw a mask on VFA to use in calculateB1MapB0GradZCorr.m
%  draw a mask on VFA and interpolate to:
    %B0 map resolution to use in the B0 map calculation
    %greEPI Res to use in calculateB1MapB0GradZCorr.m as this is very time
    %consuming without a mask (uses Bloch simulator)

%     Inputs:
%       flagCoronal - image orientation Coronal or Axial
%       gre2DEPIHeader - GRE-EPI data Dicom Header 
%       nslices_greEPI - number of slices for GRE-EPI data
%       nr_greEPI - number of rows for GRE-EPI data
%       nc_greEPI - number of columns for GRE-EPI data
%       faHeader - vfa data Dicom Header
%       faData - vfa data
%       nslices_VFA - number of VFA slices 
%       nFA_VFA - number of vfa FAs collected
%       nr_VFA - number of rows for VFA data
%       nc_VFA - number of columns for VFA data
%       B0Header - Header GRE Magnitude data collected in the B0 protocol
%       nr_B0,nc_B0 - number of rows and columns in the B0 map
%       nslices_B0 - number of slices in the B0 map 


%     Outputs:
%       maskForB1Map_gradientB0Z: Mask at the GREEPI resolution to calculate the B1+ map with slice
%                                 profile and B0 gradient through slice (Z) corrections
%       Coord_maskForB1Map_gradientB0Z: x,y coordinates of the points forming the mask
%       maskForB1Map_gradientB0Z_VFARes: Mask to calculate the B1+ map (maskForB1Map_gradientB0Z), interpolated to the VFA resolution
%       maskForB0Map: B1 map Mask interpolated to the B0 map resolution

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

%============ Step 1: Draw a Liver mask on VIBE FA15 to apply to the 15 slices of the undistorted greEPI images=========

%---------Find the VFA slice matching B1+ map slice----------------------%
[~,VFAMatchingB1Slice] = t1VFAMatchingB1Slice(gre2DEPIHeader,faHeader);

%----------0. Draw Liver mask on the vfa at greEPI reoslution----------
vfaFA15_AllSlices = zeros(nr_VFA,nc_VFA,nslices_greEPI);
for islice = 1:nslices_VFA
    vfaFA15_AllSlices(:,:,islice) = faData{islice,nFA_VFA}(:,:);
end

maskForB1Map_gradientB0Z_VFARes = zeros(nr_VFA,nc_VFA,nslices_greEPI);
Coord_maskForB1Map_gradientB0Z = cell(nslices_greEPI,1);
maskForB1Map_gradientB0Z = zeros(nr_greEPI,nc_greEPI,nslices_greEPI);
maskForB0Map = zeros(nr_B0,nc_B0,nslices_B0);
figure('Position', [350 200 1300 1000])
for  islice= 1:15%nslices_greEPI
    imagesc(squeeze(vfaFA15_AllSlices(:,:,VFAMatchingB1Slice(islice,1))),[0 200])
    colormap gray
    colorbar
    title(['Slice B1: ',num2str(islice),'; Slice VFA: ',num2str(VFAMatchingB1Slice(islice,1))])
    % Draw ROI in the liver
    ROI = drawpolygon(gca,'Color','r');
    Coord_maskForB1Map_gradientB0Z{islice,1} = ROI.Position;
    maskForB1Map_gradientB0Z_VFARes(:,:,islice)  = ROI.createMask();
    
    %----------1. Interpolate mask drawn on VFA to greEPI resolution -----------
    ImgBefInterp = maskForB1Map_gradientB0Z_VFARes(:,:,islice);
    HeaderImgBefInterp = faHeader{VFAMatchingB1Slice(islice,1),1}; %VFA SPGR
    HeaderImgAfterInterp = gre2DEPIHeader{islice,1}; %2D GRE
    
    if flagCoronal == 0
        [VFAMask_greEPIRes] = imageXYPosInterpSliceAxial(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,'VFA Mask','VFA Mask greEPI Res');
    elseif flagCoronal == 1
        [VFAMask_greEPIRes] = imageXYPosInterpSliceCoronal(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,'VFA Mask','VFA Mask greEPI Res');
        
    end
    VFAMask_greEPIRes(VFAMask_greEPIRes>0)=1;
    %interpolation will place pixels at the boundaries with values between
    %0 and 1, to get a mask we need to have logical values. Decided to
    %place values larger than 0 to 1 to avoid loosing liver when
    %downsampling the mask
    maskForB1Map_gradientB0Z(:,:,islice) = VFAMask_greEPIRes;
    
        %----------2. Interpolate mask drawn on VFA to B0 acquisition resolution -----------
    ImgBefInterp = maskForB1Map_gradientB0Z_VFARes(:,:,islice);
    HeaderImgBefInterp = faHeader{VFAMatchingB1Slice(islice,1),1}; %VFA SPGR
    HeaderImgAfterInterp = B0Header{islice,1}; %B0 Header Mag TE1
    
    if flagCoronal == 0
        [VFAMask_B0Res] = imageXYPosInterpSliceAxial(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,'VFA Mask','VFA Mask B0 Res');
    elseif flagCoronal == 1
        [VFAMask_B0Res] = imageXYPosInterpSliceCoronal(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,'VFA Mask','VFA Mask B0 Res');
        
    end
    
    VFAMask_B0Res(isnan(VFAMask_B0Res))=0;
    VFAMask_B0Res(VFAMask_B0Res>0)=1;
    %interpolation will place pixels at the boundaries with values between
    %0 and 1, to get a mask we need to have logical values. Decided to
    %place values larger than 0 to 1 to avoid loosing liver when
    %downsampling the mask. If they are placed to zero then the B0=0Hz and
    %B1+ corrected for B0gradientZ is not calculated for that slice.
    maskForB0Map(:,:,islice) = VFAMask_B0Res;


end


% Visualize VFA greEPI Resolution Mask
% figure()
% nCnt = 0;
% for isubplot = 1:nslices_greEPI
%     nCnt = nCnt+1;
%     s = subplot(3,5,nCnt);
%     plotColormap(maskForB1Map_gradientB0Z(:,:,isubplot),['Slice: ',num2str(isubplot)],[0 1],'gray','',1.5,1,s)
% end
% sgtitle('Mask for B1+ Map calculation')


end

