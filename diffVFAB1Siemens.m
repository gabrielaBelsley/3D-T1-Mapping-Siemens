function [infoFreqVoltImgPos] = diffVFAB1Siemens(Header)
%diffVFAB1Siemens 
%   extract ImagePositionPatient, ImagingFrequency, Voltages for VFA, B1+ FS,
% B1+ WE and B1+ Siemens, Shim Currents and Adjustment Volume

%   Input:
%       Header - header to extract info

%   Output:
%       infoFreqVoltImgPos - Cel with Info about Frequency, Voltages and Image
%       position locations, Shim Currents and Adjustment Volume


% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021
        
        
% Get the ASCONV Header
XProt = char(Header.Private_0029_1020)';
inds = regexp(XProt,'ASCCONV');
ASCCONV = XProt(inds(1):inds(2));



infoFreqVoltImgPos{1,1} = 'Imaging Frequency [Hz]';
infoFreqVoltImgPos{1,2} =extractFromTok(regexp(ASCCONV,'sTXSPEC.asNucleusInfo\[0\].lFrequency\s+=\s+([-+]?[0-9]+)','tokens'));


NPulses=  extractFromTok(regexp(ASCCONV,'TXSPEC.lNoOfTraPulses\s+=\s+([0-9]+)','tokens'));

% RF Pulse Name and Voltages
 for iPulse = 1:NPulses 
 % Name
infoFreqVoltImgPos{1+iPulse,1} =  [extractStringFromTok(regexp(ASCCONV,['sTXSPEC.aRFPULSE\[',num2str(iPulse-1),'\].tName\s+=\s+""(\w+(_\w+\s+\d+)?)""'],'tokens')),' [V]'];
 % Voltage
infoFreqVoltImgPos{1+iPulse,2} = extractFromTok(regexp(ASCCONV,['sTXSPEC.aRFPULSE\[',num2str(iPulse-1),'\].flAmplitude\s+=\s+(\d+(\.\d+)?)'],'tokens'));
 end
 if NPulses < 4
     blankPulses = 4-NPulses; % when the RF pulses is less than 4 (e.g. only 2 pulses), we want to leave Voltage 3 and voltage 4 with []
     NPulses = NPulses+blankPulses;
 end
 ImgPosPatient = Header.ImagePositionPatient;
infoFreqVoltImgPos{2+NPulses,1} = 'Slice Sagittal Pos [mm]';
infoFreqVoltImgPos{2+NPulses,2} = ImgPosPatient(1);
infoFreqVoltImgPos{3+NPulses,1} = 'Slice Coronal Pos [mm]';
infoFreqVoltImgPos{3+NPulses,2} = ImgPosPatient(2);
infoFreqVoltImgPos{4+NPulses,1} = 'Slice Axial Pos [mm]';
infoFreqVoltImgPos{4+NPulses,2} = ImgPosPatient(3);
infoFreqVoltImgPos{5+NPulses,1} = 'Slice Location [mm]';
infoFreqVoltImgPos{5+NPulses,2} = Header.SliceLocation;
infoFreqVoltImgPos{6+NPulses,1} = 'Overall Image Scale Factor';
infoFreqVoltImgPos{6+NPulses,2} = extractFromTok(regexp(ASCCONV,'sCoilSelectMeas.dOverallImageScaleFactor\s+=\s+(\d+\.[0-9]+)','tokens'));
%(\.\d+) or  [0-9]+, will work
infoFreqVoltImgPos{7+NPulses,1} = 'Reference Amplitude [V]';
infoFreqVoltImgPos{7+NPulses,2} =extractFromTok(regexp(ASCCONV,'sTXSPEC.asNucleusInfo\[0\].flReferenceAmplitude\s+=\s+([-+]?\d+(\.flAmplitude)?)','tokens'));

%Shim Currents: 4 shim currents
infoFreqVoltImgPos{8+NPulses,1} = 'Shim Currents [V]';
infoFreqVoltImgPos{8+NPulses,2} = extractFromTok(regexp(ASCCONV,'sGRADSPEC.alShimCurrent\[\d]\s+=\s+([-+]?[0-9]+)','tokens'));

% Volume Information: Thickness, FOV, Position (Sag,Cor,Tra)  
infoFreqVoltImgPos{9+NPulses,1} = 'Adj Volume Thickness [mm]';
infoFreqVoltImgPos{9+NPulses,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.dThickness\s+=\s+(\d+(.\d]+)?)','tokens'));
infoFreqVoltImgPos{10+NPulses,1} = 'Adj Volume Phase FOV';
infoFreqVoltImgPos{10+NPulses,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.dPhaseFOV\s+=\s+(\d+(.\d+)?)','tokens'));
infoFreqVoltImgPos{11+NPulses,1} = 'Adj Volume Readout FOV';
infoFreqVoltImgPos{11+NPulses,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.dReadoutFOV\s+=\s+(\d+(.\d+)?)','tokens'));
infoFreqVoltImgPos{12+NPulses,1} = 'Adj Volume Sagittal Pos [mm]';
infoFreqVoltImgPos{12+NPulses,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.sPosition.dSag\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoFreqVoltImgPos{13+NPulses,1} = 'Adj Volume Coronal Pos [mm]';
infoFreqVoltImgPos{13+NPulses,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.sPosition.dCor\s+=\s+([-+]?\d+(.\d+)?)','tokens'));
infoFreqVoltImgPos{14+NPulses,1} = 'c';
infoFreqVoltImgPos{14+NPulses,2} = extractFromTok(regexp(ASCCONV,'sAdjData.sAdjVolume.sPosition.dTra\s+=\s+([-+]?\d+(.\d+)?)','tokens'));



end

