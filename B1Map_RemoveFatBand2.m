function [B1Final,B1MapFatBand] = B1Map_RemoveFatBand2(B1MapFat,LiverMask)

%B1Map_RemoveFatBand2 To eliminate/mitigate the Fat band corrupting the liver B1+ factors:
% ﻿identify and remove the erroneous B1+ factor values calculated due to
% the addition of fat signal to the water signal

% Input:
    % B1MapFat: B1+ map with fat band
    % LiverMask: Mask of the liver 
% Output:
    % B1Final: B1+ map without the fat band
    % B1MapFatBand: pixels corresponding to the fat band   
    
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

addpath('auxillaryFunctionsRemoveFatBand2')

% clean up first row and save a LiverMask
B1MapFat(1,:) = NaN;
% instead of using the mask generated from ~isnan(B1MapFat) set to 1, we
% are now using the mask defined in the VFA data used in the B1B0gradientZ
% calculation as this avoids the holes in the B1+ map (with B0Z gradient
% correction) where the B1+ was not calculated and set to NaN. Using
% ~isnan(B1MapFat) the holes remain
%Pixels outside the liver were masked to NaN
% LiverMask = ~isnan(B1MapFat);
% LiverMask = imfill(LiverMask,'holes');
% figure()
% imagesc(LiverMask)


%First clean up the B1+ Map by identifying isolated extreme points 
Threshold = 0.075;
verbose = 0;
[B1MapFatNoExtrema] = EliminateIsolatedExtrema(B1MapFat, LiverMask, Threshold, verbose);


% Implementation of finding Fat Band in FatBand2 function
Thres = 2;
[B1MapFatBand] = FatBand2(B1MapFatNoExtrema, Thres, verbose);

B1FatBandNaN = B1MapFatNoExtrema;
B1FatBandNaN(B1MapFatBand) = NaN;

% Extrapolate to fill in the NaN values
B1MapNoFat = LiverMask.*inpaint_nans(B1FatBandNaN,2);
%solves a system of 2nd order differential equation=Laplacian to
%interpolate the points identified as outliers and fat band. Key idea:
%place laplacian to zero (=minimum curvature) to ensure the B1+ varies
%smoothly in the FE and PE directions.

B1MapNoFat(B1MapNoFat == 0) = NaN;



% Final Clean up Eliminate isolated extrema

Threshold = 0.05;
verbose = 0;

[B1Final] = EliminateIsolatedExtrema(B1MapNoFat, LiverMask, Threshold, verbose);



end
