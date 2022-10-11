function [imgMatchVFA] = Interpolate2VFA(B1Map,gre2DEPIHeader,vfaHeader,nslices_VFA,flagCoronal,flagGreEPI)

%INTERP2VFA  3D interpolation GREEPI to VFA resolution

%    Input:
%       B1Map: B1+ map to interpolate to VFA resolution
%       gre2DEPIHeader: Dicom Header of GRE-EPI data used to calculate the B1+ map
%       vfaHeader: vfa data header
%       nslices_VFA: number of vfa slices
%       flagCoronal: orientation of image: 1-coronal, 0-axial
%       flagGreEPI: 1 to place slice location of B1+ map equal to VFA slice
%       location
%    Output:
%       imgMatchVFA: B1+ map interpolated to VFA resolution        

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


%============ Step 1: Interpolate 2D Mulitlsice B1+/B0 Map at greEPI Resolution to VFA Resolution ===================

if flagGreEPI==1
    %---------Find the VFA slice matching B1+ map slice----------------------%
    [~,VFAMatchingB1Slice] = t1VFAMatchingB1Slice(gre2DEPIHeader,vfaHeader);
    
    if iscell(B1Map)
        [nr,nc] = size(B1Map{1,1});
        nslices = length(B1Map);
        B1Map3DMatrix = zeros(nr,nc,nslices);
        for islice = 1:nslices
            B1Map3DMatrix(:,:,islice)=B1Map{islice,1};
        end
    else
        B1Map3DMatrix = B1Map;
    end
    
    ImgBefInterp = B1Map3DMatrix;
    HeaderImgBefInterp = gre2DEPIHeader{1,1};
    HeaderImgAfterInterp = vfaHeader{1,1}; %VFA SPGR
    NslicesHighRes = nslices_VFA;
    if flagCoronal == 0
        %Important Note: added VFAMatchingB1Slice as input in
        %imageXYZPosInterpSliceAxial so that the z positions of the greEPI are
        %equal to the z positions of the vfa. 
        [imgMatchVFA] = imageXYZPosInterpSliceAxial(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,NslicesHighRes,'B1 Map','VFA',0,VFAMatchingB1Slice);
    elseif flagCoronal == 1
        [imgMatchVFA] = imageXYZPosInterpSliceCoronal(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,NslicesHighRes,'B1 Map','VFA');
        
    end
    
else
    
    if iscell(B1Map)
        [nr,nc] = size(B1Map{1,1});
        nslices = length(B1Map);
        B1Map3DMatrix = zeros(nr,nc,nslices);
        for islice = 1:nslices
            B1Map3DMatrix(:,:,islice)=B1Map{islice,1};
        end
    else
        B1Map3DMatrix = B1Map;
    end
    
    ImgBefInterp = B1Map3DMatrix;
    HeaderImgBefInterp = gre2DEPIHeader{1,1};
    HeaderImgAfterInterp = vfaHeader{1,1}; %VFA SPGR
    NslicesHighRes = nslices_VFA;
    if flagCoronal == 0
        [imgMatchVFA] = imageXYZPosInterpSliceAxial(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,NslicesHighRes,'B1 Map','VFA',0);
    elseif flagCoronal == 1
        [imgMatchVFA] = imageXYZPosInterpSliceCoronal(ImgBefInterp,HeaderImgBefInterp,HeaderImgAfterInterp,NslicesHighRes,'B1 Map','VFA');
        
    end
end
%========== End of script %========== 
end

