function [infoB1VFASiemens]=validateAcqParameters_GREEPIB0VFA(gre2DEPIHeaderFS,gre2DEPIHeaderWE,B0Header,vfaHeader,faArraygreEPIFS,faArraygreEPIWE,faArrayVFA,faB0,subject,protocolB1FS,protocolB1WE,protocolVFA,protocolB0,nFA_greEPI,nFA_VFA)

%VALIDATEACQPARAMETERS_GREEPIB0VFA: conduct acquisition patrameter checks to the greEPI Fat Sat, VIBE and B0 Map
% important before proceeding with image processing

%     Inputs:
%         gre2DEPIHeaderFS - GRE EPI Header Fat Saturation
%         gre2DEPIHeaderWE - GRE EPI Header Water Excitation
%         B0Header - B0 Map Header
%         vfaHeader - VFA SPGR Header
%         faArraygreEPIFS - Cell array with FAs acquired in GRE EPI FS: e.g. {'65','130'}
%         faArraygreEPIWE - Cell array with FAs acquired in GRE EPI WE: e.g. {'65','130'}
%         faArrayVFA - Cell array with FAs acquired in VFA SPGR: e.g. {'2','2R','15','15R'}
%         faB0 - Cell array with FAs acquired in B0: e.g. {'15'}
%         subject - data subject acquisition
%         protocolB1FS,protocolB1WE,protocolVFA,protocolB0: names of folders where the B1, B0 and VFA acquisitions were saved
%         nFA_greEPI,nFA_VFA: number of FAs in the GREEPI and VFA acquisition

%     Outputs:
%         infoB1VFASiemens: table with information extracted from headers
%                           and ASCONV of different acquisitions
%
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

%------------------------------------------------------------------------
%Check the FAs in the header match the protocol FAs 
%------------------------------------------------------------------------
faHeader_VIBE=zeros(1,nFA_VFA);
faProtocol_VIBE=zeros(1,nFA_VFA);
for iFA = 1:nFA_VFA
    faHeader_VIBE(1,iFA) = vfaHeader{1,iFA}.FlipAngle;
    if ~isnan(str2double(faArrayVFA{1,iFA}))
        faProtocol_VIBE(1,iFA) = str2double(faArrayVFA{1,iFA});
    elseif isnan(str2double(faArrayVFA{1,iFA}))
        faProtocol_VIBE(1,iFA) = str2double(faArrayVFA{1,iFA}(1:end-1));
    end
end

if any(faHeader_VIBE(1,:)~=faProtocol_VIBE(1,:))
    fprintf ('VIBE FAs acquired DO NOT match protocol \n')
end

faHeader_greEPI=zeros(1,nFA_greEPI);
faProtocol_greEPI=zeros(1,nFA_greEPI);
for iFA = 1:nFA_greEPI
    faHeader_greEPI(1,iFA) = gre2DEPIHeaderFS{1,iFA}.FlipAngle;
    if ~isnan(str2double(faArraygreEPIFS{1,iFA}))
        faProtocol_greEPI(1,iFA) = str2double(faArraygreEPIFS{1,iFA});
    elseif isnan(str2double(faArraygreEPIFS{1,iFA}))
        faProtocol_greEPI(1,iFA) = str2double(faArraygreEPIFS{1,iFA}(1:end-1));
    end
end

if any(faHeader_greEPI(1,:)~=faProtocol_greEPI(1,:))
    fprintf ('greEPI FatSat FAs acquired DO NOT match protocol \n')
end

%------------------------------------------------------------------------
%Check the information in Headers matches within each acquisition (e.g. B1+ FA60 and 120) and between acquisitions (e.g. B0 imgFreq = B1+ imgFreq)
%------------------------------------------------------------------------
% Checks only for middle slice for each RF pulse played. In principle all
% the slices played at FA = 60 should have the same imaging frequency
%Matching slices you wish to investigate the headers
% Note: I don't think there is a need to check whether the imaging
% frequency is equal for all 15 slices in B1, it should be sufficient to
% only check for middle slices of each acquisition, i.e. sliceB1 = 8 and slicevFA = 24.
sliceB1 = 8; 
sliceVFA = 24;
sliceB0 = 8;
[infoB1VFASiemens] = validateImagesHeader(subject,nFA_greEPI,nFA_VFA,gre2DEPIHeaderFS,gre2DEPIHeaderWE,B0Header,vfaHeader,faArraygreEPIFS,faArraygreEPIWE,faArrayVFA,faB0,sliceB1,sliceVFA,sliceB0,protocolB1FS,protocolB1WE,protocolVFA,protocolB0);
[nFields,nFAs] = size(infoB1VFASiemens);

%------------------------------------------------------------------------
% Flag FREQUENCY: Check that the frequency is equal between gre-epi FS, WE and B0
%------------------------------------------------------------------------
clc
FreqFS60 = infoB1VFASiemens{3,2}; 
% if ~isempty(faArraygreEPIWE) % when faArrayWE is not empty
%     FreqFSWEB0 = [infoB1VFASiemens{3,3:2+nFA_greEPI-1},infoB1VFASiemens{3,2+nFA_greEPI}];
%     FreqVIBE = infoB1VFASiemens{3,2+nFA_greEPI+1:nFAs};
%     % end corresponds to B0 protocol, only important for WE

    FreqFSWEB0 = [infoB1VFASiemens{3,3:2+nFA_greEPI-1},infoB1VFASiemens{3,2+nFA_greEPI}];
    FreqVIBE = infoB1VFASiemens{3,2+nFA_greEPI+1:nFAs};

DiffFreq = FreqFSWEB0 - FreqFS60; 
if all(DiffFreq(:) == 0) 
    fprintf('\n')
    fprintf ('FS, WE and gre B0 have matching frequencies \n')
else
    fprintf('\n')
    fprintf ('ALERT: FS, WE and gre B0 frequencies do NOT match \n')
    fprintf('FS120-FS60 = %i Hz',DiffFreq(1))
    if ~isempty(faArraygreEPIWE)
        fprintf('WE60-FS60 = %i Hz',DiffFreq(2))
        fprintf('WE120-FS60 = %i Hz',DiffFreq(3))
        fprintf('B0-FS60 = %i Hz',DiffFreq(4))
    end
end
VibeDiff = FreqVIBE -FreqFS60;
fprintf('\n')
fprintf ('Difference between VIBE frequency and gre-epi (Freq_FS60 - Freq_VFA) [Hz]: %2g %2g %2g %2g %2g \n',VibeDiff)
fprintf('\n')

%------------------------------------------------------------------------
%   Flag Reference Amplitude: Check that the Reference amplitude is equal between gre-epi and VFA
%------------------------------------------------------------------------
RefAmplitudeFS60 = infoB1VFASiemens{13,2}; 
RefAmplitude = [infoB1VFASiemens{13,3:nFAs}];
DiffRefAmp = RefAmplitude - RefAmplitudeFS60; 
if all(DiffRefAmp(:) == 0) 
    fprintf('\n')
    fprintf ('FS, WE, vfa and gre B0 have equal Reference Amplitudes \n')
else
    fprintf('\n')
    fprintf ('ALERT: FS, WE, VFA and gre B0 Reference Amplitude do NOT match \n')
     fprintf ('ALERT: Will need to divide the Nominal Angle 60 and 120 by Ratio \n')
    fprintf('FS60/120 = %i/%i [V] \n',RefAmplitudeFS60(1),RefAmplitude(1))
    if ~isempty(faArraygreEPIWE) % when faArrayWE is not empty
        fprintf('VFA = %i [V] \n',RefAmplitude(4:end-1))
        fprintf('WE60/120 = %i/%i [V] \n',RefAmplitude(2),RefAmplitude(3))
        fprintf('B0 = %i [V] \n',RefAmplitude(end))
    else % No faArrayWE 
        fprintf('VFA Ref Amplitude [V] = ')
        fprintf('%i ',RefAmplitude(2:end-1))
        fprintf('\n')
    end
end


%------------------------------------------------------------------------
% Dicom + ASCCONV info: check for differences between gre-epi alpha and 2 alpha
%------------------------------------------------------------------------
indexFa60 = 1; % corresponds to 60 degree FA, saved in first column of gre2DEPIHeader
[dicomInfoGRE60] = infoDicomHeader(gre2DEPIHeaderFS,sliceB1,indexFa60);% Extract most important dicom Header Information
[headerInfo60] = getAscconv(gre2DEPIHeaderFS{sliceB1,indexFa60},dicomInfoGRE60);
indexFa120 = 2; % corresponds to 120 degree FA, saved in second column of gre2DEPIHeader
[dicomInfoGRE120] = infoDicomHeader(gre2DEPIHeaderFS,sliceB1,indexFa120);% Extract most important dicom Header Information
[headerInfo120] = getAscconv(gre2DEPIHeaderFS{sliceB1,indexFa120},dicomInfoGRE120);


% Check the differences in the ASCCONV fields between GRE EPI 60 and 120 degrees
fprintf(1, '\n'); %Blank Line
fprintf(['Differences in ASCCONV field between GRE EPI 60',char(176), ' and ',faArraygreEPIFS{2},char(176),': \n'])
diffAscconv(headerInfo60,headerInfo120,'GRE 60','GRE 120');

%------------------------------------------------------------------------
%   Flag Shimming: Check that Shim Currents and adjustment volume are equal
%   between greEPI alpha and 2 alpha, B0 and VFA (all FAs)
%------------------------------------------------------------------------
RefShimCurrents_GREEPIFAAlpha = infoB1VFASiemens{14,2}; 
ShimCurrents = zeros(nFAs-2,length(RefShimCurrents_GREEPIFAAlpha));
DiffShimCurrents = zeros(nFAs-2,length(RefShimCurrents_GREEPIFAAlpha));
for iFA = 3:nFAs
    if isempty([infoB1VFASiemens{14,iFA}])
    else
        ShimCurrents(iFA-2,:) = [infoB1VFASiemens{14,iFA}];
        DiffShimCurrents(iFA-2,:) = RefShimCurrents_GREEPIFAAlpha-ShimCurrents(iFA-2,:);
    end
end

if all(DiffShimCurrents(:) == 0) 
    fprintf('\n')
    fprintf ('FS, WE, vfa and gre B0 have equal Shim Currents \n')
else
    fprintf('\n')
    fprintf ('ALERT: FS, WE, VFA and gre B0 Shim Currents do NOT match \n')
end

RefAdjVolume_GREEPIFAAlpha = [infoB1VFASiemens{15:nFields,2}];
AdjVolume = zeros(nFAs-2,nFields-15+1);
DiffAdjVolume = zeros(nFAs-2,nFields-15+1);
for iFA = 3:nFAs
    if isempty([infoB1VFASiemens{15:nFields,iFA}])
    else
        AdjVolume(iFA-2,:) = [infoB1VFASiemens{15:nFields,iFA}];
        DiffAdjVolume(iFA-2,:) = RefAdjVolume_GREEPIFAAlpha-AdjVolume(iFA-2,:);
    end
end

if all(DiffAdjVolume(:) == 0) 
    fprintf('\n')
    fprintf ('FS, WE, vfa and gre B0 have equal Green Adjustment Volume \n')
else
    fprintf('\n')
    fprintf ('ALERT: FS, WE, VFA and gre B0 Green Adjustment Volume CHANGED \n')
end

%------------------------------------------------------------------------
%   Flag ImagePosition: Check that Z Position is equal between greEPI alpha and 2 alpha, B0 and VFA (all FAs)
%       %Check image position X,Y is equal between greEPI alpha and 2 alpha
%       %Check image position X,Y is equal between VFA (all FAs)
%------------------------------------------------------------------------
RefImagePosPatientZ_GREEPIFAAlpha  = [infoB1VFASiemens{10,2}];
ImagePosPatientZ_GREEPIFAAlpha  = [infoB1VFASiemens{10,3:nFAs}];
DiffAdjImagePosPatientZ = RefImagePosPatientZ_GREEPIFAAlpha-ImagePosPatientZ_GREEPIFAAlpha;


if all(DiffAdjImagePosPatientZ(:) < 1e-3) 
    fprintf('\n')
    fprintf ('FS, WE, vfa and gre B0 have equal Image Position Patient Z \n')
else
    fprintf('\n')
    fprintf ('ALERT: FS, WE, VFA and gre B0 do NOT have Image Position Patient Z \n')
end

ImagePosPatientXY_GREEPIFAAlpha  = [infoB1VFASiemens{8:9,2}];
ImagePosPatientXY_GREEPIFA2Alpha  = [infoB1VFASiemens{8:9,3}];
if isequal(ImagePosPatientXY_GREEPIFAAlpha,ImagePosPatientXY_GREEPIFA2Alpha)
    fprintf('\n')
    fprintf ('FS Alpha and 2Alpha have equal Image Position Patient X,Y \n')
else
    fprintf('\n')
    fprintf ('FS Alpha and 2Alpha do NOT have equal Image Position Patient X,Y  \n')
end
if nFA_greEPI>2
    ImagePosPatientXY_GREEPIFAAlphaR  = [infoB1VFASiemens{8:9,4}];
    ImagePosPatientXY_GREEPIFA2AlphaR  = [infoB1VFASiemens{8:9,5}];
    if isequal(ImagePosPatientXY_GREEPIFAAlphaR,ImagePosPatientXY_GREEPIFA2AlphaR)
        fprintf('\n')
        fprintf ('FS AlphaR and 2AlphaR have equal Image Position Patient X,Y \n')
    else
        fprintf('\n')
        fprintf ('FS AlphaR and 2AlphaR do NOT have equal Image Position Patient X,Y  \n')
    end
end
ImagePosPatientX_VFA  = [infoB1VFASiemens{8,2+nFA_greEPI+1:nFAs}];
ImagePosPatientY_VFA  = [infoB1VFASiemens{9,2+nFA_greEPI+1:nFAs}];

if all(ImagePosPatientX_VFA==ImagePosPatientX_VFA(1)) && all(ImagePosPatientY_VFA==ImagePosPatientY_VFA(1)) 
    fprintf('\n')
    fprintf ('VFA all FAs have equal Image Position Patient X,Y \n')
else
    fprintf('\n')
    fprintf ('VFA all FAs do NOT have equal Image Position Patient X,Y  \n')
end
%------------------------------------------------------------------------
% Flag VOLTAGES: Check voltages increase linearly with FA
%------------------------------------------------------------------------
%uncomment when validating WE
% protocolB1 = protocolB1FS;
% if isequal(protocolB1,protocolB1FS)
%     flagWE = 0;
%     faArray = faArrayFS;
% elseif isequal(protocolB1,protocolB1WE)
%         flagWE = 1;
%         faArray = faArrayWE;
%         
% end
% if flagWE ==0
    VoltageFS = [infoB1VFASiemens{4,2:3}];
    FAFactor = gre2DEPIHeaderFS{1,2}.FlipAngle/gre2DEPIHeaderFS{1,1}.FlipAngle;
    if (abs(VoltageFS(1,2) - FAFactor*VoltageFS(1,1)))<0.1
        fprintf ('\n FS voltages linearly increase with FA by a factor of: %i \n',FAFactor)
    else
        fprintf('\n ALERT: FS voltages do NOT increase linearly \n')
    end
    
    if nFA_greEPI>2
        VoltageFS = [infoB1VFASiemens{4,4:5}];
        FAFactor = gre2DEPIHeaderFS{1,4}.FlipAngle/gre2DEPIHeaderFS{1,3}.FlipAngle;
        if (abs(VoltageFS(1,2) - FAFactor*VoltageFS(1,1)))<0.1
            fprintf ('\n FS voltages alphaR/2alphaR linearly increase with FA by a factor of: %i \n',FAFactor)
        else
            fprintf('\n ALERT: FS voltages do NOT increase linearly \n')
        end
    end
% else
%     VoltageWE1 = [infoB1VFASiemens{4,4:5}];
%     VoltageWE2 = [infoB1VFASiemens{5,4:5}];
%     
%     FAFactor = gre2DEPIHeader{1,2}.FlipAngle/gre2DEPIHeader{1,1}.FlipAngle;
%     if (abs(VoltageWE1(1,2) - FAFactor*VoltageWE1(1,1)))<0.1 && (abs(VoltageWE2(1,2) - FAFactor*VoltageWE2(1,1)))<0.1
%         fprintf ('\n WE voltages linearly increase with FA by: %i \n',FAFactor)
%     else
%         fprintf('\n ALERT: WE voltages do NOT increase linearly \n')
%     end
% end

%--------------------------------------------------------------------
%       Check the VFA Header information against the gre-epi
%--------------------------------------------------------------------
%Commented out as there are lots of differences as we expect from 2
%different sequences

% Dicom + ASCCONV info
% sliceVFA  = 24; % only necessary to check imaging frequencies and voltages for 1 slice of VFA, the other 48 slices should have same parameters, only the slice location should change 
% [dicomInfoVFA] = infoDicomHeader(faHeader,sliceVFA); % Most important info from DICOM Header
% [headerInfoVFA] = getAscconv(faHeader{sliceVFA,1},dicomInfoVFA);


% Check the differences in the ASCCONV fields between VFA and GRE EPI Alpha
% fprintf(1, '\n'); %Blank Line
% fprintf(['Differences in ASCCONV field between [VFA] and [GRE EPI 60]',char(176),': \n'])
% diffAscconv(headerInfo60,headerInfoVFA,'GRE 60','VFA');

%------------------------------------------------------------------------
%   Flag Voltages VFA: Check if VFA voltages increase linearly with FA
%------------------------------------------------------------------------
faArrayVFA = zeros(1,size(vfaHeader,2));
for iFA = 1:size(vfaHeader,2)
faArrayVFA(1,iFA) = vfaHeader{1,iFA}.FlipAngle;
end

%Peak Voltage B1+:
VoltageVFA = [infoB1VFASiemens{4,2+nFA_greEPI+1:nFAs}];
% The B1+ peak voltage does not increase linearly with the FA on the Prisma,
% i.e. it is stretching the pulse duration instead of increasing the
% peak B1+ voltage. It is not a problem to multiply all FAs by the same B1+ field
% factor as this is a physical property that should be independent of the
% pulse and affect equally a 3 or 15 degree FA. 

% if all(diff(VoltageVFA)==0)
%     VoltageVFA = [infoB1VFASiemens{4,6:10}];
% end
LinearVoltageIncrease = (faArrayVFA./faArrayVFA(1,1)).*VoltageVFA(1,1);
VFAVoltIncrease = VoltageVFA - LinearVoltageIncrease;
if all(abs(VFAVoltIncrease)<0.1)
    fprintf ('\n OK: VFA voltages linearly increase with FA value. \n')
else
    fprintf('\n ALERT: VFA voltages do NOT increase linearly with FA value. \n')
end
figure()
plot(faArrayVFA,VoltageVFA,'--*b')
hold on
plot(faArrayVFA,LinearVoltageIncrease,'--sm')
hold off
ylabel('Voltage [V]')
xlabel('FA [deg]')
legend('Acquisition voltage','Linear increase voltage','location','Northwest')
set(gca,'FontSize',14,'FontName','Arial')
%========== End of script %========== 
end

