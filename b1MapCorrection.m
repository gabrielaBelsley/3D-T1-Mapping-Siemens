function [B1CorrectionFactor,ratioIntAcq] = b1MapCorrection(greData,flagB1GS,flagDoubleAngle,ratioIntenSliceProfOffres,matchingB02DEPIRes,offResonance,ExtrapValue,nominalAngle,angleAmbiguityArray,ThGhostsEPI)

% b1MapCorrection Calculate the B1+ Correction Factor for each pixel with a
% correction for slice profile effects

%   Input:
%       gre2DEPIData - gre 2D signal
%       flagB1GS - set to 1 for analysis of GS B1+ data acquired with a 3D
%       acquisition, i.e. no slice profile
%       flagDoubleAngle - voltage for 2alpha is double the voltage for
%       alpha, sometimes this did not happen and the voltage ratio between
%       2alpha and alpha had to be extracted from FileMerge
%       ratioIntenSliceProfOffres - B1 look up surface (LUS) generated in script
%       b1LookupSurfaceSimulation.m
%       matchingB02DEPIRes - B0 Map scaled to the same dimensions of the gre2DEPIData
%       offResonance - range of offResonances from -300 Hz to 300 Hz in
%       steps of 100 Hz
%       ExtrapValue - NaN: pixels with a ratio larger than k or smaller
%       than 0 are placed to NaN. 
%       nominalAngle - nominal angle in degrees prescribed at the scanner
%       for alpha, used to calculate the B1+ factor=angleLUS/nominalAngle
%       angleAmbiguityArray - angle where the ratioIntenSliceProfOffres no
%       longer is injective, for the Hamming windowed it's at 95 degrees
%       (instead of 90 degrees for a sin function) due to slice profile effects
%       ThGhostsEPI - Th on the signal used to eliminate the EPI background ghosts

% For ExtrapValue started by using a large value defined by the user then cahnged to NaN.
% The large number helped to visualize pixels where the ratio was outside of the limits:
% ratio=2(small angles) and ~0 (large angles) - figure 5.2 ratio as a function of nominal FA
% 2 is the limit for small angles approximation (sinalpa~alpha) when we use k=2 (alpha,2alpha)
% Later, changed to NaN as the large number resulted in erroneous B1+ when interpolating the B1+ map to a higher resolution.
%   Output:
%       B1CorrectionFactor - B1+ Map corrected for slice profiles
%       ratioIntAcq - ratio between signals at 2alpha and alpha
%Date: 14/01/2019
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2019


[nSlices,~] =  size(greData);
[nRows,nCols] = size(greData{1,1});

% initialize variables
ratioIntAcq = cell(nSlices,1);
B1CorrectionFactor = cell(nSlices,1);
alpha = zeros(nSlices,nRows,nCols);
FAGS = cell(nSlices,1);
%B1CorrectionFactorExtrapolate = cell(nSlices,1);

for islice = 1:nSlices
    for row = 1:nRows
        for col = 1:nCols
            
            if (greData{islice,1}(row,col) <ThGhostsEPI ||  greData{islice,2}(row,col) <ThGhostsEPI )
                % Note: 24082019 eliminate ghosts from EPI images: restrict the calculation to values above a Th (ghosts in EPI have a smaller intensity)
                % Note the signal in the phantom is 500 or higher so instead of 10 we can
                % use a th of 100 and that will eliminate B1 map calculation from the background
                % background noise = low signal intensity below 10 is not considered
                B1CorrectionFactor{islice,1}(row,col) = ExtrapValue;
                ratioIntAcq{islice,1}(row,col) = greData{islice,1}(row,col);
                
            else
                numerator = double(greData{islice,2}(row,col)); %2alpha signal
                denominator = double(greData{islice,1}(row,col)); %alpha signal
                
                if flagB1GS == 1 && flagDoubleAngle == 1 % alpha2 = 2*alpha1 GS 3D GRE (Perfect Rectangular slice  OR NOT correcting for slice profile effects)
                    ratioIntAcq{islice,1}(row,col) = abs(numerator)./(2.*abs(denominator)); % when alpha2 = 2*alpha1, we can take the inverse cosine
                    
                    %Cos is limited between 0 and 1
                    if ratioIntAcq{islice,1}(row,col) > 1
                        B1CorrectionFactor{islice,1}(row,col) = ExtrapValue;
                    else
                        
                        alpha(islice,row,col) = acosd(ratioIntAcq{islice,1}(row,col));%in radians
                        FAGS{islice,1}(row,col) = alpha(islice,row,col);
                        B1CorrectionFactor{islice,1}(row,col) = alpha(islice,row,col)/(nominalAngle);
                    end
                    
                    
                elseif flagB1GS ==1 && flagDoubleAngle == 0 % GS 3D GRE but the reference voltages were not adjusted correctly at the scanner console
                    
                    ratioIntAcq{islice,1}(row,col) = abs(numerator)./(abs(denominator));
                    rfPulseAmplitude = 589.75/294.096405029; %740.055480957/425.734985352; %see from fileMerge the voltage amplitudes
                    nominal_angle = 30;
                    [ratioInten3DGRE] = lookupTable3DGRE(nominal_angle,rfPulseAmplitude);
                    B1CorrectionFactor{islice,1}(row,col)  = interp1(ratioInten3DGRE,0.01:0.01:2,ratioIntAcq{islice,1}(row,col),'linear',ExtrapValue);
                    
                    
                elseif flagB1GS == 0 % 2D GRE EPI sequence (correcting for slice profile effects)
                    
                    ratioIntAcq{islice,1}(row,col) = abs(numerator)./abs(denominator);
                    
                    if isempty(offResonance) % B1 map WITHOUT B0 off resonance correction
                        angleAmbiguity = angleAmbiguityArray(101,1); % using the Angle Ambiguity On-resonance because we assume on-resonance when we do not correct for off-resonance
                        % ratioIntenSliceProfOffres(1:angleAmbiguity,101): 101 represents the index of frequency 0Hz = on-resonance
                        angleEstimated = interp1(ratioIntenSliceProfOffres(1:angleAmbiguity,101),1:angleAmbiguity,ratioIntAcq{islice,1}(row,col),'linear',ExtrapValue);
                        B1CorrectionFactor{islice,1}(row,col) = angleEstimated./nominalAngle;
                        
                    else % WITH B0 off resonance correction
                        offresPixel = matchingB02DEPIRes(row,col,islice);
                        
                        if offresPixel < offResonance(1) || offresPixel > offResonance(end) || isnan(offresPixel)
                            %TopupB0
                            %if abs(offresPixel) < offResonance(1) || abs(offresPixel) > offResonance(end) || isnan(offresPixel)
                            B1CorrectionFactor{islice,1}(row,col) = ExtrapValue;
                        else
                            iOffres = find(offResonance == round(offresPixel));
                            angleAmbiguity = angleAmbiguityArray(iOffres);

%                            %Topup B0Map
%                            iOffres = find(offResonance == round(abs(offresPixel))); 
%                              angleAmbiguity = 120;
                            
                            % because not every signal ratio will be in the table we interpolate
                            angleEstimated = interp1(ratioIntenSliceProfOffres(1:angleAmbiguity,iOffres),1:angleAmbiguity,ratioIntAcq{islice,1}(row,col),'linear',ExtrapValue);
                            % Reason we changed from 'makima' interpolation to 'linear': see imageXYZPosInterpSliceAxial.m
                            B1CorrectionFactor{islice,1}(row,col) = angleEstimated./nominalAngle;
                            
                            
                        end
                    end
                end

            end
        end
    end
end

