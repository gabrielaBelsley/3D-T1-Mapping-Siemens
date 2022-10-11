function [B1MapFatBand] = FatBand2(B1Map, Thres, verbose)

%FatBand2 Function to identify the region of a possible Fat Band in the lower part
% of the liver due to chemical shift of Fat from the back. 
%   This function looks at the B1+ map column by column (PE direction)
%   It uses a "shooter" algorithm to predict the B1+ values based on the  
%   previous 2 points for a FPred (forward prediction) and the next 2 points for RPred (reverse prediction).
%   The rows where these two predictions are very different are indications
%   of a hard edge. When the absolute difference is greater than Thres*mean
%   over the search area these pixels are flagged as possible Fat band
%   pixels
%
%   inputs:  B1Map - the B1+ Map to be analyzed
%            Thres - Threshold in DeltaB1+ for identifying a Fat Band
%            verbose - logical flag if true then output images
%
%   outputs: B1MapFatBand - Logical Matrix with 1s covering the Fat Band

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


% Since we know the fat band is in the lower section, 
% begin the search half way down the liver  
% First find the midpoint row: rMid
NumDataPtsCol = sum(~isnan(B1Map),1);
% Identify the column with the largest number of B1+ numerical values
ThickestCol = find(NumDataPtsCol == max(NumDataPtsCol),1);
rStart = find(~isnan(B1Map(:,ThickestCol)),1);
rEnd = find(~isnan(B1Map(:,ThickestCol)),1,'last');
rMid = floor((rStart +rEnd)/2);

% Also identify columns with at least 5 B1+ data pixels below the mid row
% this will define the search area along the columns
cStart = find(sum(~isnan(B1Map(rMid:end,:))) >=5,1);
cEnd = find(sum(~isnan(B1Map(rMid:end,:))) >=5,1,'last');


% Initialize FatBand Logical Matrix
[Nr,Nc] = size(B1Map);
B1MapFatBand = false(Nr, Nc);

% search for FatBand in columns between cStart and cEnd 
for c= cStart:cEnd
    NoFat = 0;
    xBeg = rMid-2;
    xEnd = find(~isnan(B1Map(:,c)),1,'last')-3;% avoid calculating the prediction at the edges which don't have a correct B1+ value and the fat shouldn't be at the edge
    rowVect = xBeg:xEnd;
    Npts = length(rowVect);
    % Check there are at least 5 pixels between xBeg and XEnd
    if (Npts <= 5)
        continue  
    end
    FPred = zeros(Npts,1);
    RPred = zeros(Npts,1);
    
    %Calulations for each pixel in each column: 
    %Calculate predictions for the present point using the slope and
    %curvatures of the preceeding 2 points for forward prediction (FPred)
    %and the subsquent two points for the reverse prediction (RPred)
    for j = 1:Npts
        yMinus2 = B1Map(xBeg-3+j,c);
        yMinus1 = B1Map(xBeg-2+j,c);
        yPresent = B1Map(xBeg-1+j,c);
        yPlus1 = B1Map(xBeg+j,c);
        yPlus2 = B1Map(xBeg+1+j,c);
        FVel = yMinus1-yMinus2; %first derivative calculated though the backward difference
        FAcc = yPresent+yMinus2-2*yMinus1; % acceleration calculated by a central difference
        FPred(j) = yMinus1+FVel+FAcc/2;
        BVel = yPlus1-yPlus2;
        BAcc = yPresent+yPlus2-2*yPlus1;
        RPred(j) = yPlus1+BVel+BAcc/2;
    end
    
    % The metric is the contradiction: the absolute difference between these
    % two predictions FPred and RPred. If the B1+ map is homogeneous, i.e.
    % no fat band the forward and the reverse prediction should be very
    % similar
    Contradiction = abs(FPred-RPred);
    MeanContra = mean(Contradiction, 'omitnan');
    [MaxContra, MaxInd] = sort(Contradiction,'descend');
    
    % fat band was associated with pixels where the highest contradiction
    % is above twice the average absolute difference between the forward
    % and backwards predictions calculated along each column-Thres*Average
    % (Thres=2(large number) (or better could be Thres=2*std(Contradiction))) 
    % AND (&&) that the second highest difference is within 2 pixels of the
    % first highest difference. This limits the fat band extent along each
    % column to 3 pixels. The Fat band rarely spans more than 2 pixels in
    % the greEPI resolution (2*7mm)
    if (MaxContra(1)>Thres*MeanContra) && (abs(MaxInd(1)-MaxInd(2))<3)
        FatStart = min(rowVect(MaxInd(1)),rowVect(MaxInd(2)));
        FatEnd = max(rowVect(MaxInd(1)),rowVect(MaxInd(2)));
    else
        NoFat = 1; 
    end
     % Flag the appropriate pixels
        if (~NoFat)
            B1MapFatBand(FatStart:FatEnd,c) = 1;
        end

end
if verbose
    
    B1FatBandNaN = B1Map;
    B1FatBandNaN(B1MapFatBand) = NaN;
end
%     figure()
%     h = imagesc(B1FatBandNaN,[0.5 1.5]);
%     set(h, 'AlphaData', ~isnan(B1FatBandNaN))
%     set(gca, 'Color', [0, 0, 0])
%     colorbar
%     title('B1Map FatBand -> NaN')

%Rationale for this function 29052021
% calculating the gradient and the accelleration does not work because:
    %1. the fat band width varies across the columns 
    %2. fat band varies for different slices
    %3. sometimes the fat band has overestimated B1+ and sometimes
    %underestimated B1+
    %Conclusion: cannot have a general rule identifying a fat band
    %based on a global th on the gradient because this gradient will vary
    %with the width (denominator in gradient) and the intensity(numerator) of the B1+
%New approach:
    %1. calculate a forward prediciton of B1+ based on the vel and acc calculated
        %from the B1+ of the  2 rows previous to the target pixel
    %2. calculate the backward prediciton of B1+ based on the vel and acc calculated
        %from the B1+ of the 2 rows after the target pixel
    % 3. difference between the forward and the backward prediction for all pixels
        % in one column and then average. 
%    4. for each col, if the max difference is larger than 2*MeanDiff and
        % the row location of the second and the first max difference
        % is less than 3 pixels (otherwise could be an outlier or be at the
        % edge of the liver) apart that is where the fat band
        % starts. The fat band will never be larger than 2 pixels in the greEPI
        % resolution (2*7mm)

end

