function [ratioIntenSliceProfOffres,totalSignalAlpha,totalSignal2Alpha] = ratio_kAlpha_Alpha(mxy_alpha,mxy_2alpha)

%ratio_kAlpha_Alpha:
%integrate the signal across the slice profile for alpha and 2alpha
%calculate the Ratio of the integrated signals Signal2alpha/Sginalalpha 

%Input: mxy signal at each slice position for alpha and 2alpha
%Output: 
    %ratioIntenSliceProfOffres: ratio of integrated signals Signal2alpha/Sginalalpha
    %totalSignalAlpha,totalSignal2Alpha: Integrated signal at alpha and 2alpha
    
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


totalSignalAlpha= (trapz(mxy_alpha(:,1))); % complex integration across the z spatial position
totalSignal2Alpha= (trapz(mxy_2alpha(:,1))); % complex integration across the z spatial position
ratioIntenSliceProfOffres = abs(totalSignal2Alpha)./abs(totalSignalAlpha);
%Note: only take absolute value afterwards otherwise you cancel the phase of the spins
end

