function [B1Map,maxB1bias,B1MapThGhostsEPI,ThGhostsEPI] = calculateB1MapSliceProfCorr(flagCoronal,DataDir,gre2DEPIHeader,FatSup,greEPI_FAAlpha,greEPI_FA2Alpha,nslices_greEPI,varargin)

%calculateB1MapSliceProfCorr B1+ Map with Slice Profile Correction
%
%     Inputs:
%       flagCoronal - boolean flag specifying  the image orientation
%       necessary for readNifti function
%       DataDir - directory where GRE EPI lives
%       gre2DEPIHeader - Dicom Header of GRE EPI data
%       FatSup - Fat Suppression method: Fat Sat (FS) or Water excitation (WE)
%       greEPI_FAAlpha,greEPI_FA2Alpha: 'FA65', 'FA130'
%       nslices_greEPI: number of GRE EPI slices 
%       varargin: B0 Map, for FS this is not necessary in practice as the
%       signal does not change with off-resonance, but for the WE it is
%       essential

%
%     Outputs:
%         B1Map: B1+ Map after slice profile correction
%         maxB1bias: maximum B1+ factor than can be calculated using only
%         magnitude data as the phase data is not available from the scanner.
%         B1MapThGhostsEPI: B1+ Map after slice profile correction with a threshold for the EPI ghosts
%         ThGhostsEPI: Th on the signal used to eliminate the EPI background ghosts
%
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


%============ Step 1:Read the .nii undistorted images ===================



niigreEPIundistortedFld = [DataDir,'B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/greEPI/',FatSup,'NoRegist/'];
distCorrFAAlpha = readNifti([niigreEPIundistortedFld,'FA',greEPI_FAAlpha,'_unwarped.nii'],flagCoronal);
distCorrFA2Alpha = readNifti([niigreEPIundistortedFld,'FA',greEPI_FA2Alpha,'_unwarped.nii'],flagCoronal);
% %Interpolation Off
%     distCorrFAAlpha = readNifti([niigreEPIundistortedFld,'FA',greEPI_FAAlpha,'InterpOff_unwarped.nii'],flagCoronal);
%     distCorrFA2Alpha = readNifti([niigreEPIundistortedFld,'FA',greEPI_FA2Alpha,'InterpOff_unwarped.nii'],flagCoronal);

%place distCorr Data into a cell
distCorrgre2DEPIData = cell(nslices_greEPI,2);
for islice = 1:nslices_greEPI
    distCorrgre2DEPIData{islice,1} =  distCorrFAAlpha(:,:,islice); %  gre2DEPIData{islice,1};%No Dist Corr 
    distCorrgre2DEPIData{islice,2} = distCorrFA2Alpha(:,:,islice);  % gre2DEPIData{islice,2};%No Dist Corr 
end


%============ Step 2: Calculate the B1+ Correction Factor for each pixel ===================

% input variables to the b1MapCorrection.m function
nominalAngleAlpha = (gre2DEPIHeader{1,1}.FlipAngle);
centreSlicegreEPI = round(nslices_greEPI/2);

%Eliminate EPI ghosts
figure('name','Check Pixel Value for Background: Eliminate fat Chemical Shift artefacts and EPI ghosts'); 
ax = axes;
plotColormap(distCorrgre2DEPIData{centreSlicegreEPI,1},{'What is the pixel value in background?','Eliminate EPI ghosts', 'Please click on one of the ghosts to get the Pixel intensity'},[0 800],'gray','[a.u.]',0,0,ax)
[xGhost,yGhost] = ginput(1);
%note: ginput and imagesc will permute the coordinates
ThGhostsEPI = (1/5).*(distCorrgre2DEPIData{centreSlicegreEPI,1}(round(yGhost),round(xGhost))); %filter around 1/10*signal inside Phantom

%pixels with a ratio outside of the look-up surface are not extrapolated and instead set to NaN
ExtrapValue = NaN; 

flagB1GS = 0;
flagDoubleAngle =0;


%-------Alpha,2Alpha LUS: Look-up-Surface------%
SliceThickness = gre2DEPIHeader{1,1}.SliceThickness;%mm
if strcmp(FatSup(1:2),'FS')
    
    %-------Alpha,2Alpha LUS: Look-up-Surface------%
    % Slice Profile and off resonance Corrections included in look up surface - simulated offline
    if SliceThickness == 8 %mm
        load('B1MapOptimAngles/2FASliceProfOffRes300Hz_RatioFSkFactor2_Angles[1_360]deg_Gzref_Position-1_0.1_1.mat','ratioIntenSliceProfOffres','AmbiguityAngleSS')%,'offResonance')
    elseif SliceThickness == 4 %mm
        load('B1MapOptimAngles/2FASliceProfOffRes_RatioFSkFactor2_Angles[1_360]deg_Gzref_Position-1_0.1_1_SliceThick4mm.mat','ratioIntenSliceProfOffres','AmbiguityAngleSS')%,'offResonance')
    end
    
    maxB1bias = AmbiguityAngleSS(101)/nominalAngleAlpha;
    %Negate B0 Map: check this is the correct way to use off-resonance
    % imgMatchB0Neg = -imgB0MatchB1;
    %for FS don't need to use an off-resonance vector as the ratio in the LUS
    %doesn't change as a function of off-resonance 
    %Without off-resonance
    BackgroundNoisegreEPI = 10; %Background signal is about 10-30 but we want to be conservative otherwise B1 is set to NaN
    [B1Map,~] = b1MapCorrection(distCorrgre2DEPIData,flagB1GS,flagDoubleAngle,ratioIntenSliceProfOffres,[],[],ExtrapValue,nominalAngleAlpha,AmbiguityAngleSS,BackgroundNoisegreEPI);
    %Pretty B1+ Map to plot, without EPI ghosts
    B1MapThGhostsEPI = b1MapCorrection(distCorrgre2DEPIData,flagB1GS,flagDoubleAngle,ratioIntenSliceProfOffres,[],[],ExtrapValue,nominalAngleAlpha,AmbiguityAngleSS,ThGhostsEPI);
    
elseif strcmp(FatSup(1:2),'WE') %With off-resonance: Water Excitation
    
    imgMatchB0Neg=varargin{1,1};
    
    %-------Alpha,2Alpha LUS: Look-up-Surface------%
    % Slice Profile and off resonance Corrections included in look up surface - simulated offline
    load('B1MapOptimAngles/2FASliceProfOffRes300Hz_RatioWEkFactor2_Angles[1_360]deg_Position-1_0.1_1.mat','ratioIntenSliceProfOffres','AmbiguityAngleWE','offResonance')

    
    [B1Map,~] = b1MapCorrection(distCorrgre2DEPIData,flagB1GS,flagDoubleAngle,ratioIntenSliceProfOffres,imgMatchB0Neg,offResonance,ExtrapValue,nominalAngleAlpha,AmbiguityAngleWE,ThGhostsEPI);
end

%========== End of script %========== 

end

