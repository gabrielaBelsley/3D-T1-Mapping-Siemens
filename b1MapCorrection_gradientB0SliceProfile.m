function [B1CorrectionFactor,ratioIntAcq] = b1MapCorrection_gradientB0SliceProfile(greData,ExtrapValue,ThGhostsEPI,B0Map_Hz,B1Map_ConstantB0SliceProfile,T1,T2,TE,nominalAngle,kFactor,maskForB1Map_gradientB0Z,SliceThicknessFactor,FatSupress,flagPhantom)

% b1MapCorrection_gradientB0SliceProfile Take the ratio between the
% signals acquired at FA130 and 65 degrees and calculate the B1+ Factor
% corrected for slice profile effects and through slice B0 gradients
% function called by calculateB1MapB0GradZCorr.m

%     Inputs:
%       greData - GREEPI Data
%       ExtrapValue - NaN: pixels with a ratio larger than k (alpha,kalpha) or smaller
%       than 0 are placed to NaN. small angle approximation (sin(alpha)~alpha), thus sin(kalpha)/sin(alpha)~kalpha/alpha=k
%       ThGhostsEPI: Th on the signal used to eliminate the EPI background ghosts
%       B0Map_Hz: B0 Map in Hz
%       B1Map_ConstantB0SliceProfile: B1Map using only slice profile
%       correction assuming constant B0 across all the z positions within the slice thickness
%       T1,T2,TE: values for the Bloch simulation
%       nominalAngle: 65 degrees
%       kFactor: 2
%       maskForB1Map_gradientB0Z: mask for B1+ Map calculation restricted
%       to the liver - speed up B1+ factor calculation
%       SliceThicknessFactor: original protocol had 8mm->SliceThicknessFactor=1, but we later also
%       experimented with 4mm->SliceThicknessFactor=2
%       FatSupress: FS or WE
%       flagPhantom: 1 when dealing with phantom data
%
%     Outputs:
%         B1Map_B0GradZCorr: B1+ Map after slice profile correction and gradient B0-Z correction
%         ratioIntAcq: Ratio between signal intensity at FA 130 and 65 degrees
%
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

[nSlices,~] =  size(greData);
[nRows,nCols] = size(greData{1,1});

% initialize variables
ratioIntAcq = cell(nSlices,1);
B1CorrectionFactor = cell(nSlices,1);


for islice = 10%2:nSlices-1
    for row = 1:nRows
        for col = 1:nCols
            
            %Calculate B1+ only for regions inside the liver mask 
            if maskForB1Map_gradientB0Z(row,col,islice) == 0 || isnan(maskForB1Map_gradientB0Z(row,col,islice))%outside of the liver set B1+ to NaN
                B1CorrectionFactor{islice,1}(row,col) = ExtrapValue;
                ratioIntAcq{islice,1}(row,col) = greData{islice,1}(row,col);

            elseif (greData{islice,1}(row,col) <ThGhostsEPI ||  greData{islice,2}(row,col) <ThGhostsEPI || isnan(greData{islice,1}(row,col)) )
                % Note: 24082019 eliminate ghosts from EPI images: restrict the calculation to values above a Th (ghosts in EPI have a smaller intensity)
                % background noise, low signal intensity below Threshold (e.g. 10) is not considered
                % Note the signal in the phantom is 500 or higher so instead of a threshold=10 we can
                % use a threshold=100 and that will avoid B1 map calculation on the background
                B1CorrectionFactor{islice,1}(row,col) = ExtrapValue;
                ratioIntAcq{islice,1}(row,col) = greData{islice,1}(row,col);
                
            else
                numerator = double(greData{islice,2}(row,col)); %2alpha signal
                denominator = double(greData{islice,1}(row,col)); %alpha signal
                ratioIntAcq{islice,1}(row,col) = abs(numerator)./abs(denominator);
                
                
                [B1_est] = B1Factor_gradientB0SliceProfile(B0Map_Hz,B1Map_ConstantB0SliceProfile,ratioIntAcq,islice,[row,col],T1,T2,TE,nominalAngle,kFactor,SliceThicknessFactor,FatSupress,maskForB1Map_gradientB0Z,flagPhantom);
                B1CorrectionFactor{islice,1}(row,col) = B1_est;
                

                
                
            end
            
        end
    end
end


