function [B1_est,B1Guess,df] = B1Factor_gradientB0SliceProfile(B0Map,B1Map,ratioIntAcq,SliceB1,Pixel,T1,T2,TE,nominalAngle,kFactor,SliceThicknessFactor,FatSupress,maskForB1Map_gradientB0Z,flagPhantom)

%B1Factor_gradientB0SliceProfile:
%apply NLLS fit between the ratio of signals acquired at FA130 and 65 and the Bloch simulated
%ratio, simulated across z (slice direction) between -1cm and 1cm including the B0 value interpolated at steps of 0.1cm.
%Summary: Outputs B1+ corrected for a varying B0 across the slice profile

%   Input:
%         B0Map - B0 Map interpolated to greEPI resolution
%         B1Map - B1+ Map calculated with code assuming constant B0 across the slice profile
%         ratioIntAcq - ratio Signal130/Signal65
%         SliceB1 and Pixel - Slice and Pixel to recalculate the B1+ correcting for a gradient B0 across the slice profile
%         T1 - calculate T1 recovery from end of RF pulse until time TE
%         T2 - calculate T2 decay from end of RF pulse until time TE
%         TE - 11 ms
%         nominalAngle: 65 degrees
%         kFactor: 2
%         SliceThicknessFactor: original protocol had 8mm->SliceThicknessFactor=1, but we later also experimented with 4mm->SliceThicknessFactor=2
%         FatSupress: FS or WE
%         maskForB1Map_gradientB0Z: mask for B1+ Map
%         flagPhantom: 1 when dealing with phantom data

%   Output:
        %B1_est: B1+ correcting for B0 variations through slice
        %B1Guess: original B1+ calculated with code (calculateB1MapSliceProfCorr.m) assuming constant B0 across the slice profile
        %df: off-resonance gradient across the slice position

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

addpath('B0 Maps/')
addpath('B0 Maps/SLMtools/')



sliceposition = -1:0.1:1; %cm

%=============== Assign an off-resonance frequency to each slice position ============
% we only have off-resonance values at slice centre, i.e. at position -2,
% -1, 0,1,2 cm. Each slice centre is separated from the neighbour slice centre
% by 1cm = 0.4(slicethickness/2)+0.2(gap)+0.4(slicethickness/2)
% we need the offresonance values at points -1:0.1:1 cm as this is the vector
% we use in the Bloch simulations Wind_blochMxyRFInterp_gradientB0SliceProfile_Gzref.m. 

%13042021: B0Interp - new function to extrapolate B0 and enable calculation of B1+
%in places where B0 is set to 0Hz for the neighbour slices: it uses a
%global fit across all slices with B0 different from 0, and then extrapolates for slices with B0=0Hz;
B015Slices_TargetPixel = squeeze(B0Map(Pixel(1,1),Pixel(1,2),:)); %B0 for one pixel, extract values across all 15 slices
Mask15Slices_TargetPixel = squeeze(maskForB1Map_gradientB0Z(Pixel(1,1),Pixel(1,2),:)); %Mask for one pixel across all 15 slices
if sum(Mask15Slices_TargetPixel(:))>3%For a given pixel there must be at least 4 slices which contain a B0 value, otherwise B1+ is not calculated and set to NaN

    % uses a global fit across all slices with B0 different from 0, and then
    %interpolates B0 for slice positions -1:0.1:1 cm used in the Bloch simulator
    B0Fit = B0Interp(B015Slices_TargetPixel,Mask15Slices_TargetPixel,SliceB1,flagPhantom);

    
    %off-resonance gradient across the slice position, i.e. varying off-resonance with slice position
    df = B0Fit;
    %constant off-resonance for each slice position
    % df_constant = linspace(OffResSliceCentre,OffResSliceCentre,length(slicepos));% Hz
    
    
    %============ Calculate the B1+ bias whose ratio matches the acquisition ratio ========
    Ratio_Acq = ratioIntAcq{SliceB1,1}(Pixel(1,1),Pixel(1,2)); %Ratio Signal2Alpha/Alpha Acquired
    opts = optimset('Display','off');
    %lsqnonlin will find the B1Bias that introduced in the function Mxy_te_gradientB0 will output Ratio_Acq
    B1Guess = B1Map{SliceB1,1}(Pixel(1,1),Pixel(1,2)); %use as initial point the B1+ calculated with code (calculateB1MapSliceProfCorr.m) assuming constant B0 across the slice profile
    if ~isnan(B1Guess) && B1Guess>0
        B1_est = lsqnonlin(@(B1Bias)  Mxy_te_gradientB0SliceProfile(B1Bias,nominalAngle,kFactor,df,T1,T2,TE,sliceposition,SliceThicknessFactor,FatSupress)-Ratio_Acq,B1Guess,[],[],opts);
        if B1_est<0.3 || B1_est>1.7 %if B1+ is larger than ±70% of nominal FA set B1_est to NaN as this is not expected in the liver at 3T
            B1_est = NaN;
        end
    else
        B1_est = NaN;
    end
else
    B1_est = NaN;
end

end

