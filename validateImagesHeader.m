function [infoB1VFASiemens] = validateImagesHeader(subject,nFA_greEPI,nFA_VFA,gre2DEPIHeaderFS,gre2DEPIHeaderWE,B0Header,faHeader,faArrayFS,faArrayWE,faArrayVFA,faB0,sliceB1,sliceVFA,sliceB0,protocolB1FS,protocolB1WE,protocolVFA,protocolB0)

%validateImagesHeader: extract important information from the Header of the
%                       B1, B0 and VFA SPGR acquisition

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

    %intialize structure fields with Header names we want to check
    infoB1VFASiemens{1,1} = subject;
    infoB1VFASiemens{3,1} = 'Imaging Frequency [Hz]';
    infoB1VFASiemens{4,1} = 'Voltage 1 [V]';
    infoB1VFASiemens{5,1} = 'Voltage 2 [V]';
    infoB1VFASiemens{6,1} = 'Voltage 3 [V]';
    infoB1VFASiemens{7,1} = 'Voltage 4 [V]';
    infoB1VFASiemens{8,1} = 'Slice Sagittal Pos [mm]';
    infoB1VFASiemens{9,1} = 'Slice Coronal Pos [mm]';
    infoB1VFASiemens{10,1} = 'Slice Axial Pos [mm]';
    infoB1VFASiemens{11,1} = 'Slice Location [mm]';
    infoB1VFASiemens{12,1} = 'Overall Image Scale Factor';
    infoB1VFASiemens{13,1} = 'Reference Amplitude [V]';
    infoB1VFASiemens{14,1} = 'Shim Currents [V]';
    infoB1VFASiemens{15,1} = 'Adj Volume Thickness [mm]';
    infoB1VFASiemens{16,1} = 'Adj Volume Phase FOV';
    infoB1VFASiemens{17,1} = 'Adj Volume Readout FOV';
    infoB1VFASiemens{18,1} = 'Adj Volume Sagittal Pos [mm]';
    infoB1VFASiemens{19,1} = 'Adj Volume Coronal Pos [mm]';
    infoB1VFASiemens{20,1} = 'Adj Volume Coronal Pos [mm]';

    % B1+ FS
    for iFA = 1:nFA_greEPI
        infoB1VFASiemens{1,iFA+1} = strcat(protocolB1FS,string(faArrayFS(iFA)));
        [FSB1] = diffVFAB1Siemens(gre2DEPIHeaderFS{sliceB1,iFA});
        infoB1VFASiemens{2,iFA+1}=FSB1; 
        for ientries = 3:length(FSB1)+2 % the first two are titles so we only populate at field 3,fields extracted and printed in rows 3:12
            infoB1VFASiemens{ientries,iFA+1}=FSB1{ientries-2,2};
        end
    end
    
    
    % B1+ WE
    if ~isempty(faArrayWE)
    

        for iFA = 1:length(faArrayWE)
            infoB1VFASiemens{1,3+iFA} = strcat(protocolB1WE,string(faArrayWE(iFA)));
            [WEB1] = diffVFAB1Siemens(gre2DEPIHeaderWE{sliceB1,iFA});
            infoB1VFASiemens{2,3+iFA}=WEB1;
            for ientries = 3:length(WEB1)+2
                infoB1VFASiemens{ientries,2+iFA+1}=WEB1{ientries-2,2};
            end
        end
    end

    % VFA Protocol
    for iFA = 1:nFA_VFA
        infoB1VFASiemens{1,nFA_greEPI+2+iFA} = strcat(protocolVFA,string(faArrayVFA(iFA)));
        [VFAB1] = diffVFAB1Siemens(faHeader{sliceVFA,iFA});
        infoB1VFASiemens{2,nFA_greEPI+2+iFA}=VFAB1;
        for ientries = 3:length(VFAB1)+2
            infoB1VFASiemens{ientries,nFA_greEPI+2+iFA}=VFAB1{ientries-2,2};
        end
    end
    
    
    % B0 maps: LMS T2* Protocol, LMS IDEAL or 3DgreB0

    if ~isempty(protocolB0)

        infoB1VFASiemens{1,nFA_greEPI+2} = strcat(protocolB0,string(faB0(1)));
        [B0] = diffVFAB1Siemens(B0Header{sliceB0,1});
        infoB1VFASiemens{2,nFA_greEPI+2}=B0;
        for ientries = 3:length(B0)+2
            infoB1VFASiemens{ientries,nFA_greEPI+2}=B0{ientries-2,2};
        end

    end
end
