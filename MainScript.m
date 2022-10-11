% Master Script for B1+ Mapping GRE-EPI and 3D VFA-SPGR T1 Mapping for
% SIEMENS data acquired on a 3T Prisma scanner OCMR

% also calculates B1+ Map PreRF (default Siemens B1+ Map) and uses the B1+
% to correct the VFA for the T1 mapping(t1Map_preRFB1)

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022

clc;clearvars;close all

%Define Data Directories
MainFolder = '/Volumes/Seagate Expansion Drive/1.DPhil_DATA/DeLiver/ValidationHealthyVolunteers/Siemens/'; %Location of Data: Change to your directory
Subject = 'Vol5/'; 
Scanner = 'Prisma/';
Run = 'Run1/';
DataDir = [MainFolder,Subject,Scanner,Run];

%FS GRE-EPI: Slice thickness 8mm
protocolGREEPI = '2DGRE_EPI/';
FatSup = 'FS/';%'WE/GrappaOff_IntOff_PreScanOff/';
if Run(4)=='1' || Run(4)=='2'
    faArray_greEPI = {'65','130'};
elseif Run(4)=='3'
    faArray_greEPI = {'65','130','65R','130R'};
end
nslices_greEPI = 15;

%FS GRE-EPI Slice thickness 4mm
FatSup_SlTh4mm = 'FS_SliceThick4mm/';
faArray_greEPI_SlTh4mm = {'65','130'};


%WE GRE-EPI
faArrayWE = [];%WE

%B0 Map
protocolB0MagTE1 = '/B0_2DGREMultiSlice/Magnitude/TE1/';
protocolB0MagTE2 = '/B0_2DGREMultiSlice/Magnitude/TE2/';
protocolB0PhaseTE1 = '/B0_2DGREMultiSlice/Phase/TE1/';
protocolB0PhaseTE2 = '/B0_2DGREMultiSlice/Phase/TE2/';
fa_B0 = {'15'};
nslices_B0 = 15;


%3D VFA SPGR 
protocolVFA = '3DVIBE_VFA_Dixon/';
faArray_VFA = {'2','2R','15','15R'};
if Run(4)=='3'
    faArray_VFA_ControlFAs = {'3','6','9','12','15'};
end
nslices_VFA = 48;

%preRF B1+ Map: Siemens B1+ Map
protocolPreRFB1 = 'Siemens2DMultiSlice_B1Map/';
faPreRFB1 = {'80'};
nslices_PreRFB1 = 9;
SliceOrderImagPosZ_PreRFB1 = [5,9,4,8,3,7,2,6,1];

%External Functions
FolderCode_dcm2nii = genpath('xiangruili-dicm2nii-6b0c640');
addpath(FolderCode_dcm2nii)

%% Read Data acquired

%FS GRE-EPI
acqDirectory = [MainFolder,Subject,Scanner,Run,protocolGREEPI,FatSup];
[gre2DEPIData,gre2DEPIHeader] = readImage(acqDirectory,faArray_greEPI,nslices_greEPI);
[nslices_greEPI,nFA_greEPI] = size(gre2DEPIData);
[nr_greEPI,nc_greEPI] = size(gre2DEPIData{1,1});
[gre2DEPIDataMatrix] = cell2matrix(gre2DEPIData,nr_greEPI,nc_greEPI,nslices_greEPI,nFA_greEPI);

%B0 Map
%TE1
acqDirectory = [MainFolder,Subject,Scanner,Run,protocolB0MagTE1];
[B0DataMagTE1,B0HeaderMagTE1] = readImage(acqDirectory,fa_B0,nslices_B0);
acqDirectory = [MainFolder,Subject,Scanner,Run,protocolB0PhaseTE1];
[B0DataPhaseTE1,B0HeaderPhaseTE1] = readImage(acqDirectory,fa_B0,nslices_B0);
%TE2
acqDirectory = [MainFolder,Subject,Scanner,Run,protocolB0MagTE2];
[B0DataMagTE2,B0HeaderMagTE2] = readImage(acqDirectory,fa_B0,nslices_B0);
acqDirectory = [MainFolder,Subject,Scanner,Run,protocolB0PhaseTE2];
[B0DataPhaseTE2,B0HeaderPhaseTE2] = readImage(acqDirectory,fa_B0,nslices_B0);
[nr_B0,nc_B0] = size(B0DataMagTE1{1,1});
nslices_B0 = size(B0DataMagTE1,1);

%3D VFA SPGR: optimal FAs
acqDirectory = [MainFolder,Subject,Scanner,Run,protocolVFA];
[vfaData,vfaHeader] = readImage(acqDirectory,faArray_VFA,nslices_VFA);
[nslices_VFA,nFA_VFA] = size(vfaData);
[nr_VFA,nc_VFA] = size(vfaData{1,1});

%3D VFA SPGR: control FAs
if Run(4)=='3'
acqDirectory = [MainFolder,Subject,Scanner,Run,protocolVFA];
[vfaData_ControlFAs,vfaHeader_ControlFAs] = readImage(acqDirectory,faArray_VFA_ControlFAs,nslices_VFA);
[nslices_VFA,nFA_VFA_ControlFAs] = size(vfaData_ControlFAs);
end

%preRF B1+ Map: Siemens B1+ Map
acqDirectory = [MainFolder,Subject,Scanner,Run,protocolPreRFB1];
[preRFB1Data,preRFB1Header] = readImage(acqDirectory,faPreRFB1,nslices_PreRFB1);
[nslices_preRFB1,nFA_preRFB1] = size(preRFB1Data);
[nr_preRFB1,nc_preRFB1] = size(preRFB1Data{1,1});


%Orientation
if B0HeaderPhaseTE1{1,1}.ImageOrientationPatient == [1;0;0;0;0;-1]
    flagCoronal = 1;
    flagAxial = 0;
elseif B0HeaderPhaseTE1{1,1}.ImageOrientationPatient == [1;0;0;0;1;0]
    flagAxial = 1;
    flagCoronal = 0;
end


%_____________Visualize Raw greEPI + Raw VFA________________
%head,middle,foot slice 
ThreeSlicesgreEPI = [4,8,12];
ThreeSlicesVFA = [11,24,37];

RawData2Plot = cell(nFA_greEPI+nFA_VFA,3);
for islice = 1:3
    RawData2Plot{1,islice} = gre2DEPIData{ThreeSlicesgreEPI(1,islice),1};
    RawData2Plot{2,islice} = gre2DEPIData{ThreeSlicesgreEPI(1,islice),2};
    for iFA = 1:nFA_VFA
        RawData2Plot{iFA+2,islice} = vfaData{ThreeSlicesVFA(1,islice),iFA};
    end
end



close all
%Plot raw greEPI FA 65,130 degrees for 3 slices: 4,8,12
imtilePlot(RawData2Plot(1:2,:),1:6,'Slices(Row) 4,8,12 greEPI Raw Data: FA(Col) 65,130 degrees','Raw Signal [a.u.]','gray',[0 800],3,2,1,12)
%Plot raw VFA FA 2,2,15,15R degrees for 3 slices: 11,24,37
imtilePlot(RawData2Plot(3:6,:),1:12,'Slices(Row) 11,24,37 VFA Raw Data: FA(Col) 2,2R,15,15R degrees','Raw Signal [a.u.]','gray',[0 500],3,4,1,12)


%% Validate MRI Acquisition Parameters: frequencies, Reference Amplitude, Shim Currents, Adjustment Volume, ImagePosition Patient X,Y,Z, Voltages B1 and VFA
close all;clc
infoB1VFASiemens = validateAcqParameters_GREEPIB0VFA(gre2DEPIHeader,[],B0HeaderPhaseTE1,vfaHeader,faArray_greEPI,faArrayWE,faArray_VFA,fa_B0,Subject,protocolGREEPI,[],protocolVFA,protocolB0PhaseTE1,nFA_greEPI,nFA_VFA);

open infoB1VFASiemens

%% Liver Mask VFA FA15 degrees

% Draw a mask on VFA FA 15 + interpolate to:
%1. B0 resolution to use in correctDistortionEPI.m for unwrapping phase data
%2. greEPI resolution to use in calculateB1MapB0GradZCorr.m
% use the same mask for B1 Map greEPI Sl Thick 4mm or greEPI FA65R,130R
   
[maskForB1Map_gradientB0Z,Coord_maskForB1Map_gradientB0Z,maskForB1Map_gradientB0Z_VFARes,maskB0Map] = drawMaskB1GradientB0Z(flagCoronal,gre2DEPIHeader,nslices_greEPI,nr_greEPI,nc_greEPI,vfaHeader,vfaData,nslices_VFA,nFA_VFA,nr_VFA,nc_VFA,B0HeaderMagTE1,nr_B0,nc_B0,nslices_B0);

flagColorbar=0;
imtilePlot(maskForB1Map_gradientB0Z,1:nslices_greEPI,'B_{0} Map Mask B1 Resolution','Mask','parula',[0 1],3,5,flagColorbar,12)
imtilePlot(maskB0Map,1:nslices_B0,'B_{0} Map Mask B0 Resolution','Mask','parula',[0 1],3,5,flagColorbar,12)

%% Distortion correction greEPI data
flagB1AvgFAs = 0;
if nFA_greEPI ==2
    %Applying Mask to Mag TE 1 B0 Data
    [distCorrFAAlpha,distCorrFA2Alpha,B0MatchB1_Hz,B0MatchB1_NoMask_Hz] = correctDistortionEPI(flagCoronal,DataDir,gre2DEPIHeader,FatSup,faArray_greEPI,nFA_greEPI,nr_greEPI,nc_greEPI,nslices_greEPI,B0DataMagTE1,B0HeaderMagTE1,B0DataPhaseTE1,B0DataPhaseTE2,B0HeaderPhaseTE1,B0HeaderPhaseTE2,nr_B0,nc_B0,nslices_B0,maskB0Map,protocolGREEPI,flagB1AvgFAs);

elseif nFA_greEPI>2 %greEPI: Sl Thick 8mm FA65R130R
    %Applying Mask to Mag TE 1 B0 Data
    [distCorrFAAlpha,distCorrFA2Alpha,B0MatchB1_Hz,B0MatchB1_NoMask_Hz,distCorrFAAlphaR,distCorrFA2AlphaR] = correctDistortionEPI(flagCoronal,DataDir,vfaHeader,vfaData,nFA_VFA,gre2DEPIHeader,FatSup,faArray_greEPI,nFA_greEPI,nr_greEPI,nc_greEPI,nslices_greEPI,B0DataMagTE1,B0HeaderMagTE1,B0DataPhaseTE1,B0DataPhaseTE2,B0HeaderPhaseTE1,B0HeaderPhaseTE2,nr_B0,nc_B0,nslices_B0,maskB0Map,protocolGREEPI,flagB1AvgFAs);
    
    
    %greEPI: Sl Thick 4mm
    if ~isempty(gre2DEPIHeader_SlTh4mm{1,1})
    [distCorrFAAlpha_SlTh4mm,distCorrFA2Alpha_SlTh4mm] = correctDistortionEPI(flagCoronal,DataDir,vfaHeader,vfaData,nFA_VFA,gre2DEPIHeader_SlTh4mm,FatSup_SlTh4mm,faArray_greEPI_SlTh4mm,nFA_greEPI_SlTh4mm,nr_greEPI,nc_greEPI,nslices_greEPI,B0DataMagTE1,B0HeaderMagTE1,B0DataPhaseTE1,B0DataPhaseTE2,B0HeaderPhaseTE1,B0HeaderPhaseTE2,nr_B0,nc_B0,nslices_B0,maskB0Map,protocolGREEPI,flagB1AvgFAs);
    end

end


%_____________Visualize greEPI distortion corrected Sagittal and Coronal views + B0 Map________________

% Select point in the axial in the Liver: [r,c]=[SagittalSl,CoronalSl]
centreSlicegreEPI = 8;
figure()
ax = axes;
plotColormap(distCorrFAAlpha(:,:,centreSlicegreEPI),{'greEPI Undistorted','Choose col in liver (Cut Sagittal)','Choose row in liver (Cut Coronal)'},[0 800],'gray','[a.u.]',0,0,ax)

[x,y] = ginput(1);
% note: ginput as imagesc will permute the coordinates
sliceCoronalgreEPI = round(y); %row
sliceSagittalgreEPI = round(x); %col

figure('name','Sagittal View: Distorted vs Undistorted greEPI')
%FA Alpha
s1 = subplot(2,2,1);
plotColormap(fliplr(imrotate(squeeze(gre2DEPIDataMatrix(:,sliceSagittalgreEPI,:,1)),-90)),['FA ',faArray_greEPI{1,1},' DISTORTED'],[0 800],'gray','[a.u.]',0,0,s1)
s2 = subplot(2,2,3);
plotColormap(fliplr(imrotate(squeeze(distCorrFAAlpha(:,sliceSagittalgreEPI,:)),-90)),['FA ',faArray_greEPI{1,1},' UNDISTORTED'],[0 800],'gray','[a.u.]',0,0,s2)
%FA 2Alpha
s3 = subplot(2,2,2);
plotColormap(fliplr(imrotate(squeeze(gre2DEPIDataMatrix(:,sliceSagittalgreEPI,:,2)),-90)),['FA ',faArray_greEPI{1,2},' DISTORTED'],[0 800],'gray','[a.u.]',0,0,s3)
s4 = subplot(2,2,4);
plotColormap(fliplr(imrotate(squeeze(distCorrFA2Alpha(:,sliceSagittalgreEPI,:)),-90)),['FA ',faArray_greEPI{1,2},' UNDISTORTED'],[0 800],'gray','[a.u.]',0,0,s4)

figure('name','Coronal View: Distorted vs Undistorted greEPI')
%FA Alpha
s1 = subplot(2,2,1);
plotColormap(fliplr(imrotate(squeeze(gre2DEPIDataMatrix(sliceCoronalgreEPI,:,:,1)),-90)),['FA ',faArray_greEPI{1,1},' DISTORTED'],[0 800],'gray','[a.u.]',0,0,s1)
s2 = subplot(2,2,3);
plotColormap(fliplr(imrotate(squeeze(distCorrFAAlpha(sliceCoronalgreEPI,:,:,1)),-90)),['FA ',faArray_greEPI{1,1},' UNDISTORTED'],[0 800],'gray','[a.u.]',0,0,s2)
%FA 2Alpha
s3 = subplot(2,2,2);
plotColormap(fliplr(imrotate(squeeze(gre2DEPIDataMatrix(sliceCoronalgreEPI,:,:,2)),-90)),['FA ',faArray_greEPI{1,2},' DISTORTED'],[0 800],'gray','[a.u.]',0,0,s3)
s4 = subplot(2,2,4);
plotColormap(fliplr(imrotate(squeeze(distCorrFA2Alpha(sliceCoronalgreEPI,:,:,1)),-90)),['FA ',faArray_greEPI{1,2},' UNDISTORTED'],[0 800],'gray','[a.u.]',0,0,s4)

figure('name','B0 Map: Axial, Sagittal, Coronal')
for islice = 1:length(ThreeSlicesgreEPI)
s1 = subplot(2,3,islice);
plotColormap(squeeze(B0MatchB1_Hz(:,:,ThreeSlicesgreEPI(islice))),['B0 Map Axial Slice ',num2str(ThreeSlicesgreEPI(islice))],[-200 200],'parula','[Hz]',2,1,s1)
end
s2 = subplot(2,3,4);
plotColormap(fliplr(imrotate(squeeze(B0MatchB1_Hz(:,sliceSagittalgreEPI,:,1)),-90)),'B0 Map Sagittal (rows=slice direction)',[-200 200],'parula','[Hz]',1,1,s2)
s3 = subplot(2,3,5);
plotColormap(fliplr(imrotate(squeeze(B0MatchB1_Hz(sliceCoronalgreEPI,:,:,1)),-90)),'B0 Map Coronal (rows=slice direction)',[-200 200],'parula','[Hz]',1,1,s3)

%B0 Map
imtilePlot(B0MatchB1_Hz,1:nslices_greEPI,'B_{0} Map Mask','B_{0} [Hz]','parula',[-300 300],3,5,1,12)
imtilePlot(B0MatchB1_NoMask_Hz,1:nslices_greEPI,'B_{0} Map No Mask','B_{0} [Hz]','parula',[-200 200],3,5,1,12)


%% Calculate B1+ Map with slice profile correction using distortion corrected greEPI Alpha and 2Alpha
clc

greEPI_FAAlpha = faArray_greEPI{1,1}; %FA65
greEPI_FA2Alpha = faArray_greEPI{1,2}; %FA130

%greEPI: Fat Saturation Sl Thick 8mm
if strcmp(FatSup(1:2),'FS')
[B1Map_SliceProfCorr,maxB1bias,B1MapThGhostsEPI,ThGhostsEPI] = calculateB1MapSliceProfCorr(flagCoronal,DataDir,gre2DEPIHeader,FatSup,greEPI_FAAlpha,greEPI_FA2Alpha,nslices_greEPI);
flagGreEPI = 1;
[B1MapMatchVFA_SliceProf] = Interpolate2VFA(B1Map_SliceProfCorr,gre2DEPIHeader,vfaHeader,nslices_VFA,flagCoronal,flagGreEPI);
end

%greEPI: Water Excitation Sl Thick 8mm
if strcmp(FatSup(1:2),'WE')

    B1Map_SliceProfCorr_WE = calculateB1MapSliceProfCorr(flagCoronal,DataDir,gre2DEPIHeader,FatSup,greEPI_FAAlpha,greEPI_FA2Alpha,nslices_greEPI,B0MatchB1_Hz);
    flagGreEPI = 1;
    [B1MapMatchVFA_SliceProf_WE] = Interpolate2VFA(B1Map_SliceProfCorr_WE,gre2DEPIHeader,vfaHeader,nslices_VFA,flagCoronal,flagGreEPI);
end
%NOTE: to calculate B1 without slice profile correction: inside the
%function calculateB1MapSliceProfCorr place flagB1GS and flagDoubleAngle to 1



if nFA_greEPI >2
    greEPI_FAAlphaR = faArray_greEPI{1,3};
    greEPI_FA2AlphaR = faArray_greEPI{1,4};
    %greEPI: Sl Thick 8mm FA65R130R
    [B1Map_SliceProfCorr_FA65R130R] = calculateB1MapSliceProfCorr(flagCoronal,DataDir,gre2DEPIHeader,FatSup,greEPI_FAAlphaR,greEPI_FA2AlphaR,nslices_greEPI);
end


%_____________Visualize B1+ Map and ratio gre2Alpha/greAlpha____________
imtilePlot(B1Map_SliceProfCorr(1:nslices_greEPI,:),1:nslices_greEPI,'B_{1}+ Map without ThGhostsEPI','B_{1}+ Factor','parula',[0.5 1.5],3,5,1,12)
imtilePlot(B1MapThGhostsEPI(1:nslices_greEPI,:),1:nslices_greEPI,'B_{1}+ Map with ThGhostsEPI','B_{1}+ Factor','parula',[0.5 1.5],3,5,1,12)

%% Calculate B1+ Map with slice profile and B0-gradient through slice correction using the B1+ from the previous section as an initial guess

clc; close all

%place pixels with B0 less than -300 Hz to 0Hz
B0MatchB1_Hz(B0MatchB1_Hz<-300)=0;
B0MatchB1_Hz(B0MatchB1_Hz>300)=0;
BackgroundNoisegreEPI = 10;
%greEPI: Sl Thick 8mm
flagPhantom=1;
if strcmp(FatSup(1:2),'FS')
    [B1Map_B0GradZCorr,ratioIntAcq] = calculateB1MapB0GradZCorr(flagCoronal,DataDir,gre2DEPIHeader,FatSup,greEPI_FAAlpha,greEPI_FA2Alpha,nslices_greEPI,B0MatchB1_Hz,B1Map_SliceProfCorr,BackgroundNoisegreEPI,maskForB1Map_gradientB0Z,flagPhantom);

%ignore error when calculating only 1 slice in calculateB1MapB0GradZCorr
    B1Map_B0GradZCorr_AllSlices = zeros(nr_greEPI,nc_greEPI,nslices_greEPI);
    for islice = 2:nslices_greEPI-1
        B1Map_B0GradZCorr_AllSlices(:,:,islice)=B1Map_B0GradZCorr{islice,1};
    end
    flagGreEPI=1;
    [B1MapMatchVFA_NoFatCorrection] = Interpolate2VFA(B1Map_B0GradZCorr_AllSlices,gre2DEPIHeader,vfaHeader,nslices_VFA,flagCoronal,flagGreEPI);
    imtilePlot(B1MapMatchVFA_NoFatCorrection,1:nslices_VFA,'B_{1}+ Map VFA resolution','B_{1}+ Factor','parula',[0.5 1.5],6,8,1,12)
end

%greEPI: WE Sl Thick 8mm
if strcmp(FatSup(1:2),'WE')
    
    [B1Map_B0GradZCorr_WE,ratioIntAcq_WE] = calculateB1MapB0GradZCorr(flagCoronal,DataDir,gre2DEPIHeader,FatSup,greEPI_FAAlpha,greEPI_FA2Alpha,nslices_greEPI,B0MatchB1_Hz,B1Map_SliceProfCorr_WE,BackgroundNoisegreEPI,maskForB1Map_gradientB0Z,flagPhantom);
    B1Map_B0GradZCorr_WE_AllSlices = zeros(nr_greEPI,nc_greEPI,nslices_greEPI);
    for islice = 2:nslices_greEPI-1
        B1Map_B0GradZCorr_WE_AllSlices(:,:,islice)=B1Map_B0GradZCorr_WE{islice,1};
    end
    flagGreEPI=1;
    [B1MapMatchVFA_WE_NoFatCorrection] = Interpolate2VFA(B1Map_B0GradZCorr_WE_AllSlices,gre2DEPIHeader,vfaHeader,nslices_VFA,flagCoronal,flagGreEPI);
    imtilePlot(B1MapMatchVFA_WE_NoFatCorrection,1:nslices_VFA,'B_{1}+ Map VFA resolution','B_{1}+ Factor','parula',[0.5 1.5],6,8)
end


if nFA_greEPI>2
    [B1Map_B0GradZCorr_FA65R130R,ratioIntAcq_FA65R130R] = calculateB1MapB0GradZCorr(flagCoronal,DataDir,gre2DEPIHeader,FatSup,greEPI_FAAlphaR,greEPI_FA2AlphaR,nslices_greEPI,B0MatchB1_Hz,B1Map_SliceProfCorr_FA65R130R,BackgroundNoisegreEPI,maskForB1Map_gradientB0Z);
end


%_____________Visualize B1+ Map and ratio gre2Alpha/greAlpha____________
imtilePlot(B1Map_B0GradZCorr,2:nslices_greEPI-1,'B_{1}+ Map Before Removal Fat','B_{1}+ Factor','parula',[0.2 1.5],3,5,1,11)
imtilePlot(ratioIntAcq(1:(nslices_greEPI),:),2:nslices_greEPI-1,'Mxy_{2*Alpha}/Mxy_{Alpha}','[a.u.] Factor','parula',[0.5 1.5],3,5,1,12)


%% Remove fat band from B1+ Map due to chemical shift in EPI + Interpolate B1+ Map to VFA resolution (to use in T1 map calculation)

LiverMask_NoHoles = removeHolesLiverMaskFromB1Map(B1Map_B0GradZCorr,nr_greEPI,nc_greEPI,nslices_greEPI);
flagGreEPI = 1;
[B1AllCorrections_NoFat] = removeFatBandB1(B1Map_B0GradZCorr,nslices_greEPI,nr_greEPI,nc_greEPI,LiverMask_NoHoles);
[B1MapMatchVFA] = Interpolate2VFA(B1AllCorrections_NoFat,gre2DEPIHeader,vfaHeader,nslices_VFA,flagCoronal,flagGreEPI);


if nFA_greEPI >2
[B1AllCorrections_NoFat_FA65R130R] = removeFatBandB1(B1Map_B0GradZCorr_FA65R130R,nslices_greEPI,nr_greEPI,nc_greEPI,LiverMask_NoHoles);
[B1MapMatchVFA_FA65R130R] = Interpolate2VFA(B1AllCorrections_NoFat_FA65R130R,gre2DEPIHeader,vfaHeader,nslices_VFA,flagCoronal,flagGreEPI);
end


%_____________Visualize B1+ Map After Fat removal ____________
imtilePlot(B1AllCorrections_NoFat,1:nslices_greEPI,'B_{1}+ Map Chemical Shift Fat Band removed','B_{1}+ Factor','parula',[0.2 1.5],3,5,1,11)

%_____________Visualize Liver Mask derived from B1+ Map without holes ____________
imtilePlot(LiverMask_NoHoles,1:nslices_greEPI,'Liver Mask: No holes','Logical','parula',[0 1],3,5,1,11)

%_____________Visualize B1+ Map at high resolution VFA____________
imtilePlot(B1MapMatchVFA,1:nslices_VFA,'B_{1}+ Map VFA resolution','B_{1}+ Factor','parula',[0.2 1.5],8,8,1,11)

[B0MapMatchVFA] = Interpolate2VFA(B0MatchB1_Hz,gre2DEPIHeader,vfaHeader,nslices_VFA,flagCoronal,flagGreEPI);
%_____________Visualize B0 Map at high resolution VFA____________
imtilePlot(B0MapMatchVFA,1:nslices_VFA,'B_{0} Map VFA resolution Masked','B_{0} [Hz]','parula',[-100 200],6,8,1,11)


%% preRF B1+ Map: Siemens default B1+ map

preRFB1_Ordered = ones(nr_preRFB1,nc_preRFB1,nslices_preRFB1);
preRFB1Header_Ordered = cell(nslices_preRFB1,1);
for s = 1:nslices_preRFB1
    preRFB1_Ordered(:,:,s) = preRFB1Data{SliceOrderImagPosZ_PreRFB1(s),1}(:,:)./(preRFB1Header{SliceOrderImagPosZ_PreRFB1(s),1}.FlipAngle.*10);
    preRFB1Header_Ordered{s,1} = preRFB1Header{SliceOrderImagPosZ_PreRFB1(s),1};
end
%_____________Visualize B1+ Map preRF____________
imtilePlot(preRFB1_Ordered,1:nslices_preRFB1,'preRF B_{1}+ Map','B_{1}+ Factor','parula',[0.5 1.5],2,5,1,12)

flagGreEPI = 0;
[B1MapMatchVFA_preRFB1] = Interpolate2VFA(preRFB1_Ordered,preRFB1Header_Ordered,vfaHeader,nslices_VFA,flagCoronal,flagGreEPI);

flagGreEPI = 1;
[maskForB1Map_VFARes] = Interpolate2VFA(maskForB1Map_gradientB0Z,gre2DEPIHeader,vfaHeader,nslices_VFA,flagCoronal,flagGreEPI);
maskForB1Map_VFARes(maskForB1Map_VFARes~=1)=0;
imtilePlot(maskForB1Map_VFARes,1:nslices_VFA,'B_{1}+ Map mask at VFA resolution','Logical','parula',[0 1],6,8,1,12)

B1MapMatchVFA_preRFB1_masked = ones(nr_VFA,nc_VFA,nslices_VFA);
for islice = 1:nslices_VFA
B1MapMatchVFA_preRFB1_masked(:,:,islice) = B1MapMatchVFA_preRFB1(:,:,islice).*maskForB1Map_VFARes(:,:,islice);
end

%_____________Visualize B1+ Map preRF at high resolution VFA____________
imtilePlot(B1MapMatchVFA_preRFB1_masked,1:nslices_VFA,'preRF B_{1}+ Map VFA resolution','B_{1}+ Factor','parula',[0.5 1.5],6,8,1,12)



%% Calculate T1 using B1+ map and incomplete spoiling correction


clc;close all

% Code to interpolate just a single slice: used when calculating only 1 B1+
% map slice in function calculateB1MapSliceProfCorr due to time constraints. 
SliceB1=10;
SliceVFA=31;%VFA slice matching Z location of SliceB1
ImgBefInterp = B1AllCorrections_NoFat(:,:,SliceB1);
HeaderImgBefInterp = gre2DEPIHeader{SliceB1,1};
HeaderImgAfterInterp = vfaHeader{SliceVFA,1};
[B1MapMatchVFAsingleSlice] = imageXYPosInterpSliceAxial(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,'B1 Map','vfa'); 
B1MapMatchVFA = zeros(nr_VFA,nc_VFA,nslices_VFA);
B1MapMatchVFA(:,:,SliceVFA)=B1MapMatchVFAsingleSlice;
figure()
ax=axes;
plotColormap(B1MapMatchVFA(:,:,SliceVFA),'B1+ Map VFA Res',[0.5 1.5],'parula','B1+ Factor',0,0,ax)


%-----greEPI: FATSAT Sl Thick 8mm---------
%***B1 Map SLice Profile Correction + B0 Gradient Z Correction + NO INCOMP SPOIL Correction***
flagInCompSpoilCorr=1;
slicesVFA_T1Map = SliceVFA;%1:nslices_VFA; for all slices
[m0Map,b1Map,t1Map,rsqMap,ciT1Map,NormSSE,deltaT1,Norm2Residuals,ResidualsFit,CorrFactorSpoil,lbT1,ubT1] = calculateT1Map(vfaHeader,vfaData,nFA_VFA,maxB1bias,B1MapMatchVFA,flagInCompSpoilCorr,slicesVFA_T1Map);


%--------preRF B1+--------
[m0Map_preRFB1,b1Map_preRFB1,t1Map_preRFB1,rsqMap_preRFB1,ciT1Map_preRFB1,NormSSE_preRFB1,deltaT1_preRFB1,Norm2Residuals_preRFB1,ResidualsFit_preRFB1,CorrFactorSpoil_preRFB1,lbT1_preRFB1,ubT1_preRFB1] = calculateT1Map(vfaHeader,vfaData,nFA_VFA,maxB1bias,B1MapMatchVFA_preRFB1_masked,flagInCompSpoilCorr,slicesVFA_T1Map);

%_____________Visualize T1 Map Single slice____________
%B1+ GEEPI
figure()
ax=axes;
plotColormap(t1Map(:,:,SliceVFA),'T_1 Map w/ GE-EPI B_1^+',[500 1500],'parula','T_1 (ms)',1.2,0,ax,1,'',200,0)
figure()
ax=axes;
plotColormap(B1MapMatchVFA(:,:,SliceVFA),'B_1^+ Map GE-EPI',[0.5 1.5],'parula','B_1^+ factor',1.2,0,ax,1,'',0.2,0)

%B1+ PreRF
figure()
ax=axes;
plotColormap(t1Map_preRFB1(:,:,SliceVFA),'T_1 Map w/ PreRF B_1^+',[500 1500],'parula','T_1 (ms)',1.2,0,ax,1,'',200,0)
figure()
ax=axes;
plotColormap(B1MapMatchVFA_preRFB1_masked(:,:,SliceVFA),'B_1^+ Map PreRF',[0.5 1.5],'parula','B_1^+ factor',1.2,0,ax,1,'',0.2,0)
%%
%_____________Visualize T1 Map All slices____________
imtilePlot(t1Map,1:nslices_VFA,'T_{1} Map','T_{1} [ms]','parula',[lbT1 ubT1],6,8,1,12)

imtilePlot(t1Map_preRFB1,1:nslices_VFA,'T_{1} Map','T_{1} [ms]','parula',[lbT1 ubT1],6,8,1,12)
