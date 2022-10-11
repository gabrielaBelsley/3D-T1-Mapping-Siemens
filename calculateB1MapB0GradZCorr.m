function [B1Map_B0GradZCorr,ratioIntAcq] = calculateB1MapB0GradZCorr(flagCoronal,DataDir,gre2DEPIHeader,FatSup,greEPI_FAAlpha,greEPI_FA2Alpha,nslices_greEPI,B0MatchB1_NoMask_Hz,B1Map_SliceProfCorr,ThGhostsEPI,maskForB1Map_gradientB0Z,flagPhantom)

%calculateB1MapB0GradZCorr B1+ Map with through slice off-resonance correction
%
%     Inputs:
%       flagCoronal - boolean flag specifying  the image orientation
%       necessary for readNifti function
%       DataDir - directory where GRE EPI lives
%       gre2DEPIHeader - Dicom Header of GRE EPI data
%       FatSup - Fat Suppression method: Fat Sat (FS) or Water excitation (WE)
%       greEPI_FAAlpha,greEPI_FA2Alpha: 'FA65', 'FA130'
%       nslices_greEPI: number of GRE EPI slices 
%       B0MatchB1_NoMask_Hz: B0 Map at the B1+ map resolution
%       B1Map_SliceProfCorr: B1+ Map with slice profile correction to use as initial guess for the B1+ Map with B0-Z correction
%       ThGhostsEPI: Th on the signal used to eliminate the EPI background ghosts
%       maskForB1Map_gradientB0Z: B1+ Map mask 
%       flagPhantom: boolean, for Phantom Data place to 1
%
%     Outputs:
%         B1Map_B0GradZCorr: B1+ Map after slice profile correction and gradient B0-Z correction
%         ratioIntAcq: Ratio between signal intensity at FA 130 and 65 degrees
%
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


%============ Step 1:Read the .nii undistorted images ===================


niigreEPIundistortedFld = [DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/greEPI/',FatSup,'NoRegist/'];
distCorrFAAlpha = readNifti([niigreEPIundistortedFld,'FA',greEPI_FAAlpha,'_unwarped.nii'],flagCoronal);
distCorrFA2Alpha = readNifti([niigreEPIundistortedFld,'FA',greEPI_FA2Alpha,'_unwarped.nii'],flagCoronal);
% distCorrFAAlpha = readNifti([niigreEPIundistortedFld,'FA',greEPI_FAAlpha,'InterpOff_unwarped.nii'],flagCoronal);
% distCorrFA2Alpha = readNifti([niigreEPIundistortedFld,'FA',greEPI_FA2Alpha,'InterpOff_unwarped.nii'],flagCoronal);

%place distCorr Data into a cell
distCorrgre2DEPIData = cell(nslices_greEPI,2);
for islice = 1:nslices_greEPI
    distCorrgre2DEPIData{islice,1} = distCorrFAAlpha(:,:,islice); %gre2DEPIData{islice,1};%No Dist Corr     
    distCorrgre2DEPIData{islice,2} = distCorrFA2Alpha(:,:,islice); %gre2DEPIData{islice,2};%No Dist Corr  
end


%============ Step 2: Calulate B1+ map taking into account B0 gradient across Z ========================

SliceThickness = gre2DEPIHeader{1,1}.SliceThickness;%mm
OriginalProtocolSlThick = 8;%mm
SliceThicknessFactor = OriginalProtocolSlThick/SliceThickness;

nominalAngleAlpha = deg2rad(gre2DEPIHeader{1,1}.FlipAngle);
nominalAngle2Alpha = deg2rad(gre2DEPIHeader{1,2}.FlipAngle);
kFactor = nominalAngle2Alpha/nominalAngleAlpha;%k*alpha
T1 = 900; %ms
T2 = 30; %ms

TE = gre2DEPIHeader{1,1}.EchoTime; %ms %time from the centre of the RF pulse until the echo

ExtrapValue = NaN;
[B1Map_B0GradZCorr,ratioIntAcq] = b1MapCorrection_gradientB0SliceProfile(distCorrgre2DEPIData,ExtrapValue,ThGhostsEPI,B0MatchB1_NoMask_Hz,B1Map_SliceProfCorr,T1,T2,TE,nominalAngleAlpha,kFactor,maskForB1Map_gradientB0Z,SliceThicknessFactor,FatSup(1:2),flagPhantom);


%========== End of script %========== 

end

