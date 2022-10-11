function [distCorrFAAlpha,distCorrFA2Alpha,imgB0MatchB1_Hz,imgB0MatchB1_NoMask_Hz,varargout]=correctDistortionEPI(flagCoronal,DataDir,gre2DEPIHeader,FatSup,faArray_greEPI,nFA_greEPI,nr_greEPI,nc_greEPI,nslices_greEPI,B0DataMagTE1,B0HeaderMagTE1,B0DataPhaseTE1,B0DataPhaseTE2,B0HeaderPhaseTE1,B0HeaderPhaseTE2,nr_B0,nc_B0,nslices_B0,maskB0Map,protocolGREEPI,flagB1AvgFAs)

%CORRECTDISTORTIONEPI distortion correction of greEPI Data using B0 map and fsl fugue

%   Input:
%       flagCoronal - image orientation Coronal or Axial
%       DataDir - directory where data lives
%       gre2DEPIHeader - GRE-EPI data Dicom Header 
%       FatSup - Fat Suppression method: Fat Saturation (FS) or Water excitation (WE)
%       faArray_greEPI - FAs used to acquire GRE EPI data: 65 and 130 degrees 
%       nFA_greEPI - number of GRE-EPI FAs collected
%       nr_greEPI - number of rows for GRE-EPI data
%       nc_greEPI - number of columns for GRE-EPI data
%       nslices_greEPI - number of slices for GRE-EPI data
%       B0DataMagTE1 - GRE Magnitude data collected at the first TE in the B0 protocol
%       B0HeaderMagTE1 - Header of GRE Magnitude data collected at the first TE in the B0 protocol
%       B0DataPhaseTE1 - GRE Phase data collected at the first TE in the B0 protocol
%       B0DataPhaseTE2 - GRE Phase data collected at the second TE in the B0 protocol
%       B0HeaderPhaseTE1 - Header of GRE Phase data collected at the first TE in the B0 protocol
%       B0HeaderPhaseTE2 - Header of GRE Phase data collected at the second TE in the B0 protocol
%       nr_B0 - number of rows in the B0 map
%       nc_B0 - number of columns in the B0 map
%       nslices_B0 - number of slices in the B0 map 
%       maskB0Map - liver mask B0 map 
%       protocolGREEPI - name of GRE-EPI forlder '2DGRE_EPI/'
%       flagB1AvgFAs - 1 when averaging data from two acquisitions of FA65

%   Output:
%       distCorrFAAlpha: GRE EPI FA Alpha distortion corrected
%       distCorrFA2Alpha: GRE EPI FA 2Alpha distortion corrected
%       imgB0MatchB1_Hz: B0 Map in Hz matching the GREEPI resolution with Liver Mask
%       imgB0MatchB1_NoMask_Hz: B0 Map in Hz matching the GREEPI resolution without Liver Mask
%       varargout: in case repetitions of the GRE EPI FA Alpha and 2 Alpha
%       were acquired, it also ouputs the GRE EPI FA AlphaR/2AlphaR distortion corrected
%
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

%-------Prelude and Fugue User Guide--------------------
%Uses fsl prelude to unwrap the phase and calculate the B0 Maps 
%uses fsl fugue for distortion Correction of GRE-EPI data
%https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide
%------------------------------------------------------

%============ Step 1: Draw Liver mask for each of the B0DataMagTE1 15 slices,apply Liver Mask to B0DataMagTE1 and write to .nii  ===================

%---------Find the VFA slice matching B1+ map slice----------------------%
%should also work for the B0 map as the B0 and B1 have matching slice locations
% [~,VFAMatchingB1Slice] = t1VFAMatchingB1Slice(gre2DEPIHeader,vfaHeader);
% 
% maskedROIs = zeros(nr_B0,nc_B0,nslices_B0);
% CenterCoord = cell(nslices_B0,1);
% figure('Position', [150 200 1600 1000])
% for  islice= 8%%1:nslices_B0
%     subplot(1,2,2)
%     imagesc(squeeze(vfaData{VFAMatchingB1Slice(islice,1),nFA_VFA}),[0 300])
%     axis image
%     subplot(1,2,1)
%     imagesc(squeeze(B0DataMagTE1{islice,1}),[0 300])
%     axis image
%     zoom(1.3)
%     colormap gray
%     title(['Slice: ',num2str(islice)])
%     Select ROI with mouse
%     ROI = drawpolygon(gca,'Color','r');
%     CenterCoord{islice,1} = ROI.Position;
%     maskedROIs(:,:,islice)  = ROI.createMask();
% 
% 
% end

%Apply the mask to the magnitude data TE1
MagTE1Liver = zeros(nr_B0,nc_B0,nslices_B0);
for  islice= 1:nslices_B0
    MagTE1Liver(:,:,islice) = maskB0Map(:,:,islice).*B0DataMagTE1{islice,1};
end
% imtilePlot(maskB0Map,1:nslices_B0,'B_{0} Map Mask','B_{0} [Hz]','parula',[0 1],3,5)

%Check how the mask fits the 15 slices for both the MagTE1 and PhaseTE1
% for  islice= 1:nslices_B0
%     figure('Position', [350 200 1000 800])
%     subplot(1,2,1)
%     imagesc(B0DataMagTE1{islice,1})
%     hold on
%     drawpolygon('Position',CenterCoord{islice,1} ,'Color','red');
%     hold off
%     subplot(1,2,2)
%     imagesc(B0DataPhaseTE1{islice,1})
%     hold on
%     drawpolygon('Position',CenterCoord{islice,1} ,'Color','red');
%     hold off
% end

%Save the masked Magnitude TE1 to .nii for the B0 Map calculation
dcmFolder = [DataDir,'B0_2DGREMultiSlice/Magnitude/TE1/FA15/'];
niiFolderMag_TE1 = [DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/MagTE1/'];
writeNifti(MagTE1Liver,'MagTE1LiverMask',niiFolderMag_TE1,niiFolderMag_TE1,'MagTE1',flagCoronal,dcmFolder)

%============ Step 2:  Map Phase from [-4095 4096] to [0 2pi] and write to .nii ===================
Phase_TE1 = zeros(nr_B0,nc_B0,nslices_B0);
Phase_TE2 = zeros(nr_B0,nc_B0,nslices_B0);
phase02pi_TE1 = zeros(nr_B0,nc_B0,nslices_B0);
phase02pi_TE2 = zeros(nr_B0,nc_B0,nslices_B0);


for iSlice = 1: nslices_B0
    for r = 1:nr_B0
        for c = 1:nc_B0
            %[-4095, 4096]
            Phase_TE1(r,c,iSlice) = (B0HeaderPhaseTE1{iSlice,1}.RescaleSlope).*B0DataPhaseTE1{iSlice,1}(r,c)+(B0HeaderPhaseTE1{iSlice,1}.RescaleIntercept);
            Phase_TE2(r,c,iSlice) = (B0HeaderPhaseTE2{iSlice,1}.RescaleSlope).*B0DataPhaseTE2{iSlice,1}(r,c)+(B0HeaderPhaseTE2{iSlice,1}.RescaleIntercept);
            
            % between [0 2pi]
            phase02pi_TE1(r,c,iSlice) =  Phase_TE1(r,c,iSlice) .* (2*pi/(4096+4095))+ 4095*(2*pi/(4096+4095)); %get the phase evolution at each pixel
            phase02pi_TE2(r,c,iSlice) =  Phase_TE2(r,c,iSlice) .* (2*pi/(4096+4095))+ 4095*(2*pi/(4096+4095)); %get the phase evolution at each pixel
        end
    end
end

%Histogram of PhaseTE1 and PhaseTE2 to check the phase is in range [0 2pi]
% islice = 8;
% figure()
% subplot(1,2,1)
% histogram(phase02pi_TE1(:,:,islice))
% title('TE1: Check the phase data is between 0 and 2pi')
% subplot(1,2,2)
% histogram(phase02pi_TE2(:,:,islice))
% title('TE2: Check the phase data is between 0 and 2pi')
% min(phase02pi_TE1(:,:,islice),[],'all')
% max(phase02pi_TE1(:,:,islice),[],'all')
% min(phase02pi_TE2(:,:,islice),[],'all')
% max(phase02pi_TE2(:,:,islice),[],'all')


%----------Write the Phase Data in range [0 2pi] to .nii for both TE1 and TE2-------------------
% Write Phase02pi_TE1 to .nii 
dcmFolder = [DataDir,'B0_2DGREMultiSlice/Phase/TE1/FA15/'];
niiFolderPhase_TE1 = [DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/PhaseTE1/'];
writeNifti(phase02pi_TE1,'Phase02pi_TE1',niiFolderPhase_TE1,niiFolderPhase_TE1,'PhaseTE1',flagCoronal,dcmFolder)


% Write Phase02pi_TE2 to .nii 
dcmFolder = [DataDir,'B0_2DGREMultiSlice/Phase/TE2/FA15/'];
niiFolderPhase_TE2 = [DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/PhaseTE2/'];
writeNifti(phase02pi_TE2,'Phase02pi_TE2',niiFolderPhase_TE2,niiFolderPhase_TE2,'PhaseTE1',flagCoronal,dcmFolder)

%============Step 3: Convert greEPI FAAlpha, FA2Alpha .dcm to .nii ========================

%FA Alpha, 2Alpha (AlphaR, 2AlphaR)
niiFolder_greEPI = [DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/greEPI/',FatSup];
for iFA = 1:nFA_greEPI
dcmFolder_FAAlpha = [DataDir,protocolGREEPI,FatSup,'FA',faArray_greEPI{1,iFA},'/'];
convertDicom2Nii(dcmFolder_FAAlpha,[niiFolder_greEPI,'FA',faArray_greEPI{1,iFA},'/'],['FA',faArray_greEPI{1,iFA}])
end

setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI'); % this to tell what the output type would be

if flagB1AvgFAs ==1 %when averaging data from two acquisitions of FA65 and 2 acquisitions of FA130
    cd([DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/'])

    %Average FA 65 and FA 65R
    AvgFAs65 = '/usr/local/fsl/bin/fslmaths greEPI/FS/FA65/FA65 -add greEPI/FS/FA65R/FA65R -div 2  greEPI/FS/FA65Avg';
    system(AvgFAs65);
    %Average FA 130 and FA 130R
    AvgFAs130 = '/usr/local/fsl/bin/fslmaths greEPI/FS/FA130/FA130 -add greEPI/FS/FA130R/FA130R -div 2  greEPI/FS/FA130Avg';
    system(AvgFAs130);
end
%============Step 4: fsl unwrap phase using prelude ==================================

cd([DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/'])

%No Mask
UnwrapPhase1 = '/usr/local/fsl/bin/prelude -a MagTE1/MagTE1 -p PhaseTE1/Phase02pi_TE1 -o PhaseTE1/Phase02pi_TE1_Unwrap';
system(UnwrapPhase1);
UnwrapPhase2 = '/usr/local/fsl/bin/prelude -a MagTE1/MagTE1 -p PhaseTE2/Phase02pi_TE2 -o PhaseTE2/Phase02pi_TE2_Unwrap';
system(UnwrapPhase2);

%Using Mask defined in step 1
UnwrapPhase1_Masked = '/usr/local/fsl/bin/prelude -a MagTE1/MagTE1LiverMask -p PhaseTE1/Phase02pi_TE1 -o PhaseTE1/Phase02pi_TE1_Unwrap_Mask';
system(UnwrapPhase1_Masked);
UnwrapPhase2_Masked = '/usr/local/fsl/bin/prelude -a MagTE1/MagTE1LiverMask -p PhaseTE2/Phase02pi_TE2 -o PhaseTE2/Phase02pi_TE2_Unwrap_Mask';
system(UnwrapPhase2_Masked);

%============ Step 5: fsl calculate B0 Map ========================================
%DeltaTE = (teList(2,1) - teList(1,1)); %seconds
%To do include DeltaTE in command fslmaths instead of hard coding as it is
%different for Trio
calculateB0Map = '/usr/local/fsl/bin/fslmaths PhaseTE2/Phase02pi_TE2_Unwrap -sub PhaseTE1/Phase02pi_TE1_Unwrap -mul 1000 -div 2.39 fieldmap_rads -odt float';
system(calculateB0Map);
%Mask
calculateB0MapMask = '/usr/local/fsl/bin/fslmaths PhaseTE2/Phase02pi_TE2_Unwrap_Mask -sub PhaseTE1/Phase02pi_TE1_Unwrap_Mask -mul 1000 -div 2.39 fieldmap_rads_Mask -odt float';
system(calculateB0MapMask);

%=================== Step 6: Read the .nii fieldmap created from unwraped phase data ===================
cd('/Users/gabrielabelsley/OneDrive - Perspectum Ltd/DPhil_SecondProject/DataScriptsResults/3DT1MappingScriptsFunctions/3DT1MapSiemensGit/')
B0MapFSl = readNifti([DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/fieldmap_rads_Mask.nii'],flagCoronal);
B0MapFSl_NoMask = readNifti([DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/fieldmap_rads.nii'],flagCoronal);

%=================== Step 7: Interpolate B0 Map to gre 2D EPI Data resolution ===================
%----------CHECK whether we need interpolation in the Slice direction: ----------------------
% in principle not as the B1 and B0 acquisitions have been matched in terms of number of slices, slice gap and slice thickness. 
sliceLocationsB0 = zeros(nslices_B0,1);
sliceLocationsB1 = zeros(nslices_B0,1);
% Check that the slice locations between B0 and gre B1 Map match
for iSliceB0 = 1:nslices_B0
    sliceLocationsB0(iSliceB0,1) =  B0HeaderMagTE1{iSliceB0,1}.SliceLocation;
end

for iSliceB1 = 1:nslices_greEPI
    sliceLocationsB1(iSliceB1,1) =  gre2DEPIHeader{iSliceB1,1}.SliceLocation;
end

nslices_compare=min(nslices_B0,nslices_greEPI);

if any(sliceLocationsB0(1:nslices_compare,1)-sliceLocationsB1(1:nslices_compare,1)> 1e-3)
    fprintf ('B1 and B0 slice locations are NOT equal, NEED INTERPOLATION in the slice direction \n')
    flag2DInterpolation = 0;
else
    fprintf ('B1 and B0 slice locations are equal, do not need interpolation in the slice direction \n')
    flag2DInterpolation = 1;
end
format bank
fprintf ('Slice Locations B0, B1 \n')
for iSliceB1 = 1:nslices_compare
    disp(['B0: ',num2str(sliceLocationsB0(iSliceB1,1)),' B1: ',num2str(sliceLocationsB1(iSliceB1,1))])
end

imgB0MatchB1= zeros(nr_greEPI,nc_greEPI,nslices_greEPI);
imgB0MatchB1_NoMask= zeros(nr_greEPI,nc_greEPI,nslices_greEPI);
if flag2DInterpolation == 1
    for iSliceB0Map = 1:nslices_B0
        ImgBefInterp = B0MapFSl(:,:,iSliceB0Map);
        ImgBefInterp2 = B0MapFSl_NoMask(:,:,iSliceB0Map);
        HeaderImgBefInterp = B0HeaderMagTE1{iSliceB0Map,1};
        % given the slice locations match we can use the same index (iSliceB0Map) for the imgB0MatchB1
        % to iterate through the slices
        HeaderImgAfterInterp = gre2DEPIHeader{iSliceB0Map,1};
        if flagCoronal == 0
            [imgB0MatchB1(:,:,iSliceB0Map)] = imageXYPosInterpSliceAxial(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,'B0 Map','greEPI Alpha');
            [imgB0MatchB1_NoMask(:,:,iSliceB0Map)] = imageXYPosInterpSliceAxial(ImgBefInterp2,HeaderImgBefInterp,HeaderImgAfterInterp,'B0 Map No Mask','greEPI Alpha');
            
        elseif flagCoronal == 1
            [imgB0MatchB1(:,:,iSliceB0Map)] = imageXYPosInterpSliceCoronal(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,'B0 Map','greEPI Alpha');
            [imgB0MatchB1_NoMask(:,:,iSliceB0Map)] = imageXYPosInterpSliceCoronal(ImgBefInterp2,HeaderImgBefInterp,HeaderImgAfterInterp,'B0 Map No Mask','greEPI Alpha');
            
        end
        
    end
    %-----------------------------------------------------------------
    % when we need 3D interpolation along the slice direction as the locations between gre and B0 do not match
else
    
    ImgBefInterp = B0MapFSl;
    ImgBefInterp2 = B0MapFSl_NoMask;
    HeaderImgBefInterp = B0HeaderMagTE1{1,1};
    HeaderImgAfterInterp = gre2DEPIHeader{1,1};
    NslicesHighRes=nslices_greEPI;
    if flagCoronal == 0
        [imgB0MatchB1] = imageXYZPosInterpSliceAxial(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,NslicesHighRes,'B0 Map','GRE EPI 60');
        [imgB0MatchB1_NoMask] = imageXYZPosInterpSliceAxial(ImgBefInterp2,HeaderImgBefInterp,HeaderImgAfterInterp,NslicesHighRes,'B0 Map No Mask','GRE EPI 60');
        
    elseif flagCoronal == 1
        [imgB0MatchB1] = imageXYZPosInterpSliceCoronal(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,NslicesHighRes,'B0 Map','GRE EPI 60');
        [imgB0MatchB1_NoMask] = imageXYZPosInterpSliceCoronal(ImgBefInterp2,HeaderImgBefInterp,HeaderImgAfterInterp,NslicesHighRes,'B0 Map No Mask','GRE EPI 60');
        
        
    end
    
end

%=================== Step 8: Write the B0Map at greEPI resolution to .nii =================== 

niigreEPIFAAlphaFolder = [DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/greEPI/',FatSup,'/FA',faArray_greEPI{1,1},'/'];
niiFolder = [DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/'];
writeNifti(imgB0MatchB1,'B0MapgreEpiRes',niiFolder,niigreEPIFAAlphaFolder,['FA',faArray_greEPI{1,1}],flagCoronal)

imgB0MatchB1_Hz = imgB0MatchB1./(2*pi);%to get to Hz
imgB0MatchB1_NoMask_Hz = imgB0MatchB1_NoMask./(2*pi);%to get to Hz

%=================== Step 9: fsl distortion correction with fugue ===============================

%========== Calulating Echo spacing fugue ==========
BWPerPixelPhaseEncode = gre2DEPIHeader{1, 1}.Private_0019_1028;
NumbPELines = double(gre2DEPIHeader{1, 1}.Rows);
% % if Interpolation is On this has to be 104(ReconMatrix PE lines) and not 52 (AcqMatrix PE lines) 
% format long
EchoSpacing = 1/(BWPerPixelPhaseEncode*NumbPELines);
% TotalReadoutTime = 1./BWPerPixelPhaseEncode

cd([DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/'])

if flagB1AvgFAs ==0
    %FA 60
    %unwarpFAAlpha = ['/usr/local/fsl/bin/fugue -i greEPI/',FatSup,'/FA60/FA60 --dwell=',num2str(EchoSpacing),' --loadfmap=B0MapgreEpiRes --unwarpdir=y- -u greEPI/',FatSup,'/NoRegist/FA60_unwarped'];
    %FA 65
    unwarpFAAlpha = ['/usr/local/fsl/bin/fugue -i greEPI/',FatSup,'/FA65/FA65 --dwell=',num2str(EchoSpacing),' --loadfmap=B0MapgreEpiRes --unwarpdir=y- -u greEPI/',FatSup,'/NoRegist/FA65_unwarped'];
    system(unwarpFAAlpha);
    
    %FA 120
    %unwarpFA2Alpha = ['/usr/local/fsl/bin/fugue -i greEPI/',FatSup,'/FA120/FA120 --dwell=',num2str(EchoSpacing),' --loadfmap=B0MapgreEpiRes --unwarpdir=y- -u greEPI/',FatSup,'/NoRegist/FA120_unwarped'];
    %FA 130
    unwarpFA2Alpha = ['/usr/local/fsl/bin/fugue -i greEPI/',FatSup,'/FA130/FA130 --dwell=',num2str(EchoSpacing),' --loadfmap=B0MapgreEpiRes --unwarpdir=y- -u greEPI/',FatSup,'/NoRegist/FA130_unwarped'];
    system(unwarpFA2Alpha);
    
    
    if nFA_greEPI>2 %there is FA65R and FA130R to unwarp
        %unwarpFAAlphaR = ['/usr/local/fsl/bin/fugue -i greEPI/',FatSup,'/FA65R/FA65R --dwell=0.0003 --loadfmap=B0MapgreEpiRes --unwarpdir=y- -u greEPI/',FatSup,'/NoRegist/FA65R_unwarped'];
            unwarpFAAlphaR = ['/usr/local/fsl/bin/fugue -i greEPI/',FatSup,'/FA',faArray_greEPI{1,1},'/FA',faArray_greEPI{1,1} '--dwell=0.0003 --loadfmap=B0MapgreEpiRes --unwarpdir=y- -u greEPI/',FatSup,'/NoRegist/FA65R_unwarped'];
        system(unwarpFAAlphaR);
        
        unwarpFA2AlphaR = ['/usr/local/fsl/bin/fugue -i greEPI/',FatSup,'/FA130R/FA130R --dwell=0.0003 --loadfmap=B0MapgreEpiRes --unwarpdir=y- -u greEPI/',FatSup,'/NoRegist/FA130R_unwarped'];
        system(unwarpFA2AlphaR);
   
    end
else
    unwarpFAAlphaAvg = ['/usr/local/fsl/bin/fugue -i greEPI/',FatSup,'/FA65Avg --dwell=0.0003 --loadfmap=B0MapgreEpiRes --unwarpdir=y- -u greEPI/',FatSup,'/NoRegist/FA65Avg_unwarped'];
    system(unwarpFAAlphaAvg);
    
    unwarpFA2AlphaAvg = ['/usr/local/fsl/bin/fugue -i greEPI/',FatSup,'/FA130Avg --dwell=0.0003 --loadfmap=B0MapgreEpiRes --unwarpdir=y- -u greEPI/',FatSup,'/NoRegist/FA130Avg_unwarped'];
    system(unwarpFA2AlphaAvg);
end

%Return to script directory
cd('/Users/gabrielabelsley/OneDrive - Perspectum Ltd/DPhil_SecondProject/DataScriptsResults/3DT1MappingScriptsFunctions/3DT1MapSiemensGit/')

niigreEPIundistortedFld = [DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/greEPI/',FatSup,'NoRegist/'];
if flagB1AvgFAs ==0
distCorrFAAlpha = readNifti([niigreEPIundistortedFld,'FA',faArray_greEPI{1,1},'_unwarped.nii'],flagCoronal);
distCorrFA2Alpha = readNifti([niigreEPIundistortedFld,'FA',faArray_greEPI{1,2},'_unwarped.nii'],flagCoronal);
if nFA_greEPI>2 %there is FA65R and FA130R to unwarp
    distCorrFAAlphaR = readNifti([niigreEPIundistortedFld,'FA',faArray_greEPI{1,3},'_unwarped.nii'],flagCoronal);
    varargout{1,1} = distCorrFAAlphaR;
    distCorrFA2AlphaR = readNifti([niigreEPIundistortedFld,'FA',faArray_greEPI{1,4},'_unwarped.nii'],flagCoronal);
    varargout{1,2} = distCorrFA2AlphaR;
            
end
else
    distCorrFAAlpha = readNifti([niigreEPIundistortedFld,'FA',faArray_greEPI{1,1},'Avg_unwarped.nii'],flagCoronal);
    distCorrFA2Alpha = readNifti([niigreEPIundistortedFld,'FA',faArray_greEPI{1,2},'Avg_unwarped.nii'],flagCoronal);
end


%========== End of script %========== 
end

