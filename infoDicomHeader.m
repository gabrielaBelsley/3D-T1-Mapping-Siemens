function [dicomInfo] = infoDicomHeader(dicomHeader,sliceB1Match,varargin)

%infoDicomHeader: extract most important information from the Dicom Header


% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

    %Most Important Dicom Headers to understand signal  
    
    if nargin == 3
        faHeader = varargin{1};
    else
        faHeader = 1;
    end
    
    dicomInfo{1,1} = 'Echo Time [ms]';
    dicomInfo{1,2} = dicomHeader{sliceB1Match,faHeader}.EchoTime;
    dicomInfo{2,1} = 'Repetition Time [ms]';
    dicomInfo{2,2} = dicomHeader{sliceB1Match,faHeader}.RepetitionTime;
    dicomInfo{3,1} = ['Flip Angle [',char(176),']'];
    dicomInfo{3,2} = dicomHeader{sliceB1Match,faHeader}.FlipAngle;
    dicomInfo{4,1} = 'FOV';
    dicomInfo{4,2} = dicomHeader{sliceB1Match,faHeader}.Private_0051_100c; %FOV
    dicomInfo{5,1} = 'Slice Thickness';
    dicomInfo{5,2} = dicomHeader{sliceB1Match,faHeader}.SliceThickness;
    dicomInfo{6,1} = 'Pixel Spacing';
    dicomInfo{6,2} = dicomHeader{sliceB1Match,faHeader}.PixelSpacing;
    dicomInfo{7,1} = 'Number of Rows';
    dicomInfo{7,2} = dicomHeader{sliceB1Match,faHeader}.Rows;
    dicomInfo{8,1} = 'Number of Columns';
    dicomInfo{8,2} = dicomHeader{sliceB1Match,faHeader}.Columns;
    dicomInfo{9,1} = 'Image Position Patient';
    dicomInfo{9,2} = dicomHeader{sliceB1Match,faHeader}.ImagePositionPatient;
    dicomInfo{10,1} = 'FatSat or WE';
    dicomInfo{10,2} = dicomHeader{sliceB1Match,faHeader}.ScanOptions;%WE;
    dicomInfo{11,1} = 'MR Acquisition Type';
    dicomInfo{11,2} = dicomHeader{sliceB1Match,faHeader}.MRAcquisitionType;
    dicomInfo{12,1} = 'Sequence';
    dicomInfo{12,2} = dicomHeader{sliceB1Match,faHeader}.SequenceName;
    dicomInfo{13,1} = 'Frequency [MHz]';
    dicomInfo{13,2} = dicomHeader{sliceB1Match,faHeader}.ImagingFrequency;
    
%     sliceLoc = [];
%     [nSlices,~] = size(dicomHeader);
%     for iSlice = 1:nSlices
%            sliceLoc = cat(1,sliceLoc, dicomHeader{iSlice,faHeader}.SliceLocation);
%     end
    dicomInfo{14,1} = 'Slice Location';
    dicomInfo{14,2} = dicomHeader{sliceB1Match,faHeader}.SliceLocation;
end

