function [B1SliceMatchingVFA,VFAMatchingB1Slice,sliceLocationsB1,sliceLocationsVFA] = t1VFAMatchingB1Slice(headerGRE,headerVFA)

%Input: 
    %VFA Header
    %greEPI Header
%Output: 
    %   B1SliceMatchingVFA: find the B1 matching slice to the VFA slice acquired: B1SliceMatchingVFA
    %   VFAMatchingB1Slice: find the vfa matching slice to the B1 slice acquired: VFAMatchingB1Slice
    %  sliceLocationsB1: vector with slice locations for B1+ map acquisition
    %  sliceLocationsVFA: vector with slice locations for VFA map acquisition

        
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

    

    [nSlicesT1,~] = size(headerVFA);
    [nSlicesB1,~] = size(headerGRE);
    sliceLocationsB1 = zeros(nSlicesB1,1);
    sliceLocationsVFA = zeros(nSlicesT1,1);
    VFAMatchingB1Slice = zeros(nSlicesB1,1);
    B1SliceMatchingVFA = zeros(nSlicesT1,1);

    
    for iSliceB1 = 1:nSlicesB1
        sliceLocationsB1(iSliceB1,1) =  headerGRE{iSliceB1,1}.ImagePositionPatient(3);
    end
    
    for iSliceT1 = 1:nSlicesT1
           sliceLocationsVFA(iSliceT1,1) = headerVFA{iSliceT1,1}.ImagePositionPatient(3);
    end
    
    
    % for each B1+ slice go through all VFA slices and find the closest VFA slice
    for iSliceB1 = 1:nSlicesB1 
        % subtract VFA slice locations from gre single slice location
        diffSliceLoc_VFA = sliceLocationsVFA-ones(nSlicesT1,1).*headerGRE{iSliceB1,1}.ImagePositionPatient(3); % subtract vfa locations from B1 single slice location
        [~,VFAMatchingB1Slice(iSliceB1,1)] = min(abs(diffSliceLoc_VFA));
        VFAMatchingB1Slice(iSliceB1,2) = headerVFA{VFAMatchingB1Slice(iSliceB1,1),1}.ImagePositionPatient(3); %Z Poisition VFA
        VFAMatchingB1Slice(iSliceB1,3) = headerGRE{iSliceB1,1}.ImagePositionPatient(3); %Z Poisition greEPI
    end
        
    % for each T1 slice go through all B1 slices and find the closest B1 slice
    for iSliceT1 = 1:nSlicesT1 
        % subtract B1 (gre) slice locations from VFA single slice location
        diffSliceLoc_GRE = sliceLocationsB1-ones(nSlicesB1,1).*headerVFA{iSliceT1,1}.ImagePositionPatient(3); 
        [~,B1SliceMatchingVFA(iSliceT1,1)] = min(abs(diffSliceLoc_GRE)); 
        %28082019: instead of less than 2 (abs(diffSliceLoc) < 2) change to minimum
        %check where the difference between locations is less than 2, it
        % should be less than the slice thickness of the vfa which is 3 mm,
        % because if the found vfa matching slice is more than 3 mm
        % distance relative to the gre location, it means there is a closer
        % vfa slice that we should be using instead.
        
        B1SliceMatchingVFA(iSliceT1,2) = headerGRE{B1SliceMatchingVFA(iSliceT1,1),1}.ImagePositionPatient(3); %Z Poisition greEPI
        B1SliceMatchingVFA(iSliceT1,3) = headerVFA{iSliceT1,1}.ImagePositionPatient(3); %Z Poisition VFA
    end
    %fprintf('GRE EPI 60 matching VFA slice # is: %d \n', vfaMatchingB1Slice)

end

