function [infoAscconv] = getAscconv(dicomHeader, dicomInfo)

%getAscconv: extract ASCONV (has more information than Dicom Headers such as B1 voltages)

%     Inputs:
%       dicomHeader - Header of dicom data
%       dicomInfo - Info extracted from dicom data using function infoDicomHeader

%     Outputs:
%       infoAscconv: asconv information

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

% Acess important inormation in ASCCONV of dicom images through Regular Expressions

% Matlab explanation: https://uk.mathworks.com/help/matlab/matlab_prog/regular-expressions.html
% A regular expression is a sequence of characters that defines a certain
% pattern. You normally use a regular expression to search text for a group
% of words that matches the pattern, for example, while parsing program
% input or while processing a block of text.

% Websites for regular expressions:
% http://regexlib.com/Search.aspx?k=decimal&AspxAutoDetectCookieSupport=1
% https://medium.com/factory-mind/regex-tutorial-a-simple-cheatsheet-by-examples-649dc1c3f285
% https://www.regular-expressions.info/tutorial.html
% https://developer.mozilla.org/en-US/docs/Web/JavaScript/Guide/Regular_Expressions

% Usage:
% extract ASCCONV
% in here, info is the output from dicominfo(dicomfilename).

%extractStringFromTok(regexp(XProt,'"%CustomerSeq%\\([A-Za-z0-9_]+)"','tokens'))
%extractFromTok(regexp(XProt,'sTXSPEC\.asNucleusInfo\[0\]\.lFrequency\s+=\s+([0-9]+)','tokens'))

lengthDicomInfo = length(dicomInfo);
infoAscconv = dicomInfo;
XProt = char(dicomHeader.Private_0029_1020)';
inds = regexp(XProt,'ASCCONV');
ASCCONV = XProt(inds(1):inds(2));

%Sequence Name
infoAscconv{1+lengthDicomInfo,1} = 'Sequence Name';

SeqName = regexp(ASCCONV,'tSequenceFileName\s+=\s+""%SiemensSeq%\\(\w+)""','tokens');

if isempty(SeqName)
    SeqName = regexp(ASCCONV,'tProtocolName\s+=\s+""(\w+)""','tokens');
end
if isempty(SeqName)
SeqName = regexp(ASCCONV,'tSequenceFileName\s+=\s+""%CustomerSeq%\\(\w+)""','tokens');
end
infoAscconv{1+lengthDicomInfo,2} = extractStringFromTok(SeqName);
  %  infoAscconv{1+lengthDicomInfo,2} = extractStringFromTok(regexp(ASCCONV,'tSequenceFileName\s+=\s+""%SiemensSeq%\\(\w+[0-9]\w+_\w+)""','tokens'));

%Ref SNR
infoAscconv{2+lengthDicomInfo,1} = 'Reference SNR';
infoAscconv{2+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'dRefSNR\s+=\s+([0-9]+.[0-9]+)','tokens'));

% Flip Angle
infoAscconv{3+lengthDicomInfo,1} = ['Flip Angle [',char(176),']'];
infoAscconv{3+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'adFlipAngleDegree\[0\]\s+=\s+([0-9]+.[0-9]+)','tokens'));

%B0 Nominal
infoAscconv{4+lengthDicomInfo,1} = 'B0 Nominal';
infoAscconv{4+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'\w+\.\w+B0\s+=\s+([0-9]+.[0-9]+)','tokens'));

%Shim Currents: 4 shim currents
infoAscconv{5+lengthDicomInfo,1} = 'Shim Currents [V]';
infoAscconv{5+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sGRADSPEC.alShimCurrent\[\d]\s+=\s+([-+]?[0-9]+)','tokens'));

% Number of pulses: parallel imaging with ref lines at start have one more
% pulse than ref lines integrated or no reference lines
infoAscconv{6+lengthDicomInfo,1} = 'Number of Pulses';
infoAscconv{6+lengthDicomInfo,2} =  extractFromTok(regexp(ASCCONV,'TXSPEC.lNoOfTraPulses\s+=\s+([0-9]+)','tokens'));

%B1 Shim mode
infoAscconv{7+lengthDicomInfo,1} = 'B1 Shim Mode';
infoAscconv{7+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'TXSPEC.lB1ShimMode\s+=\s+([0-9]+)','tokens'));

% frequency
infoAscconv{8+lengthDicomInfo,1} = 'Imaging Frequency [Hz]';
infoAscconv{8+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sTXSPEC.asNucleusInfo\[0\].lFrequency\s+=\s+([-+]?[0-9]+)','tokens'));

% Reference Voltage
infoAscconv{9+lengthDicomInfo,1} = 'Reference Voltage [V]';
infoAscconv{9+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sTXSPEC.asNucleusInfo\[0\].flReferenceAmplitude\s+=\s+([-+]?\d+(\.flAmplitude)?)','tokens'));

% B1 Reference Voltage
infoAscconv{10+lengthDicomInfo,1} = 'B1 Reference Voltage [V]';
infoAscconv{10+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sTXSPEC.asNucleusInfo\[0\].flCompProtectionB1PlusRefAmpl\s+=\s+([-+]?\d+(\.\d+)?)','tokens'));



% B1 Shim Type + B0 Correction Type
infoAscconv{11+lengthDicomInfo,1} = 'B1 Shim Type';
infoAscconv{12+lengthDicomInfo,1} = 'B0 Correction Type';

if extractFromTok(regexp(ASCCONV,'sTXSPEC.aPTXRFPulse.__attribute__.size\s+=\s+([0-9]+)','tokens')) == 0
    infoAscconv{11+lengthDicomInfo,2} = [];
    infoAscconv{12+lengthDicomInfo,2} = [];
else
    infoAscconv{11+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sTXSPEC.aPTXRFPulse\[0\].lB1ShimType\s+=\s+([0-9]+)','tokens'));
    infoAscconv{12+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sTXSPEC.aPTXRFPulse\[0\].lB0CorrectionType\s+=\s+([0-9]+)','tokens'));
end



% Volume Information: Thickness, FOV, Position (Sag,Cor,Tra)  
infoAscconv{13+lengthDicomInfo,1} = 'Adj Volume Thickness [mm]';
infoAscconv{13+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.dThickness\s+=\s+(\d+(.\d]+)?)','tokens'));
infoAscconv{14+lengthDicomInfo,1} = 'Adj Volume Phase FOV';
infoAscconv{14+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.dPhaseFOV\s+=\s+(\d+(.\d+)?)','tokens'));
infoAscconv{15+lengthDicomInfo,1} = 'Adj Volume Readout FOV';
infoAscconv{15+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.dReadoutFOV\s+=\s+(\d+(.\d+)?)','tokens'));
infoAscconv{16+lengthDicomInfo,1} = 'Adj Volume Sagittal Pos [mm]';
infoAscconv{16+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.sPosition.dSag\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{17+lengthDicomInfo,1} = 'Adj Volume Coronal Pos [mm]';
infoAscconv{17+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.sPosition.dCor\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{18+lengthDicomInfo,1} = 'Adj Volume Axial Pos [mm]';
infoAscconv{18+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.sPosition.dTra\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{19+lengthDicomInfo,1} = 'Axial Adj Volume?';
infoAscconv{19+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.sNormal.dTra\s+=\s+([-+]?\d+(.\d+)?)','tokens'));



% Slice Information: Thickness, FOV, Position (Sag,Cor,Tra)  
infoAscconv{20+lengthDicomInfo,1} = 'Slice Thickness [mm]';
infoAscconv{20+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sSliceArray.asSlice\[0\].dThickness\s+=\s+(\d+(.\d]+)?)','tokens'));
infoAscconv{21+lengthDicomInfo,1} = 'Slice Phase FOV';
infoAscconv{21+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sSliceArray.asSlice\[0\].dPhaseFOV\s+=\s+(\d+(.\d+)?)','tokens'));
infoAscconv{22+lengthDicomInfo,1} = 'Slice Readout FOV';
infoAscconv{22+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sSliceArray.asSlice\[0\].dReadoutFOV\s+=\s+(\d+(.\d+)?)','tokens'));
infoAscconv{23+lengthDicomInfo,1} = 'Slice Sagittal Pos [mm]';
infoAscconv{23+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sSliceArray.asSlice\[0\].sPosition.dSag\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{24+lengthDicomInfo,1} = 'Slice Coronal Pos [mm]';
infoAscconv{24+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sSliceArray.asSlice\[0\].sPosition.dCor\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{25+lengthDicomInfo,1} = 'Slice Axial Pos [mm]';
infoAscconv{25+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sSliceArray.asSlice\[0\].sPosition.dTra\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{26+lengthDicomInfo,1} = 'Axial Slice?';
infoAscconv{26+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sSliceArray.asSlice\[0\].sNormal.dTra\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
% Distance between slices
infoAscconv{27+lengthDicomInfo,1} = 'Distance between Slices [%]';
infoAscconv{27+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sGroupArray.asGroup\[0\].dDistFact\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
% Gap between slices in mm 
infoAscconv{28+lengthDicomInfo,1} = 'Gap between Slices [mm]';
infoAscconv{28+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sGroupArray.sPSat.dGap\s+=\s+([-+]?\d+(.\d+)?)','tokens'));




%K-space Lines
infoAscconv{29+lengthDicomInfo,1} = 'k-Space Phase Partial Fourier?';
infoAscconv{29+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sKSpace.dSeqPhasePartialFourierForSNR\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{30+lengthDicomInfo,1} = 'k-Space Base Res';
infoAscconv{30+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sKSpace.lBaseResolution\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{31+lengthDicomInfo,1} = 'k-Space Phase Encoding Lines';
infoAscconv{31+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sKSpace.lPhaseEncodingLines\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{32+lengthDicomInfo,1} = 'k-Space Partitions';
infoAscconv{32+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sKSpace.lPartitions\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{33+lengthDicomInfo,1} = 'k-Space Images Per Slab';
infoAscconv{33+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sKSpace.lImagesPerSlab\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{34+lengthDicomInfo,1} = 'k-Space Radial Views';
infoAscconv{34+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sKSpace.lRadialViews	\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{35+lengthDicomInfo,1} = 'k-Space Radial Interleaves per Image';
infoAscconv{35+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sKSpace.lRadialInterleavesPerImage\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{36+lengthDicomInfo,1} = 'k-Space Lines Per Shot:';
infoAscconv{36+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sKSpace.lLinesPerShot\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{37+lengthDicomInfo,1} = 'k-Space Phase Partial Fourier';
infoAscconv{37+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sKSpace.ucPhasePartialFourier\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{38+lengthDicomInfo,1} = 'k-Space Slice Partial Fourier';
infoAscconv{38+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sKSpace.ucSlicePartialFourier\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{39+lengthDicomInfo,1} = 'k-Space Averaging Mode';
infoAscconv{39+lengthDicomInfo,2} =extractFromTok(regexp(ASCCONV,'sKSpace.ucAveragingMode\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{40+lengthDicomInfo,1} = 'k-Space EPI Factor';
infoAscconv{40+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sFastImaging.lEPIFactor\s+=\s+([-+]?\d+(.\d+)?)','tokens'));




%Parallel Imaging
infoAscconv{41+lengthDicomInfo,1} = 'Accelleration Factor PE';
infoAscconv{41+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sPat.lAccelFactPE\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{42+lengthDicomInfo,1} = 'Accelleration Factor 3D';
infoAscconv{42+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sPat.lAccelFact3D\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{43+lengthDicomInfo,1} = 'Reference Lines PE';
infoAscconv{43+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sPat.lRefLinesPE\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{44+lengthDicomInfo,1} = 'PAT Mode';
infoAscconv{44+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sPat.ucPATMode\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{45+lengthDicomInfo,1} = 'Reference Scan Mode';
infoAscconv{45+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sPat.ucRefScanMode\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoAscconv{46+lengthDicomInfo,1} = 'Average all frames?';
infoAscconv{46+lengthDicomInfo,2} = extractStringFromTok(regexp(ASCCONV,'sPat.ucTPatAverageAllFrames\s+=\s+([-+]?\d+\w\d+)','tokens'));
infoAscconv{47+lengthDicomInfo,1} = 'Caipirinha Mode';
infoAscconv{47+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sPat.ulCaipirinhaMode\s+=\s+([-+]?\d+(.\d+)?)','tokens'));

% Overall Image Scale Factor
infoAscconv{48+lengthDicomInfo,1} = 'Overall Image Scale Factor';
infoAscconv{48+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,'sCoilSelectMeas.dOverallImageScaleFactor\s+=\s+(\d+\.[0-9]+)','tokens'));

 % RF Pulse Name and Voltages
 for iPulse = 1:infoAscconv{6+lengthDicomInfo,2} 
 % Name
infoAscconv{48+iPulse+lengthDicomInfo,1} =  [extractStringFromTok(regexp(ASCCONV,['sTXSPEC.aRFPULSE\[',num2str(iPulse-1),'\].tName\s+=\s+""(\w+(\s+)?\w+)""'],'tokens')),' [V]'];
 % Voltage
infoAscconv{48+iPulse+lengthDicomInfo,2} = extractFromTok(regexp(ASCCONV,['sTXSPEC.aRFPULSE\[',num2str(iPulse-1),'\].flAmplitude\s+=\s+(\d+(\.\d+)?)'],'tokens'));
 end
 
 
 
end
