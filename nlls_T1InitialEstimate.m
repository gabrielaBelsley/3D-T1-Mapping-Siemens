function [T1_InitialEst] = nlls_T1InitialEstimate(DataPerSlicePerPixel,faTrueRad,TR,lbT1,ubT1)


%nlls_T1InitialEstimate Initial T1 estimate from the ratio of two signals
%acquired at distinct FAs

%   Input:
        %DataPerSlicePerPixel: vector with data for a single pixel. Length=Number of FAs collected
        %faTrueRad: FA corrected for B1+ inhomogeneities
        %TR: repetition time in ms
        %lbT1,ubT1: lower and upper T1 bounds

%   Output:
        %T1_InitialEst: Initial T1 estimate
        
        
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


%identify the two lowest FAs in the acquisition
Data_FALowest = DataPerSlicePerPixel((1)); %Lowest FA
indexFALowest = 1;

index_FASecondLowest = find(diff(faTrueRad)>0,1);
index_FASecondLowest = index_FASecondLowest+1;
Data_FASecondLowest = DataPerSlicePerPixel(index_FASecondLowest); % Second Lowest FA

MeanT1UpperLowerBound = (ubT1+lbT1)/2; %average between ubT1 and lbT1

%----------Estimate T1 from the two lowest angles: less affected by spoiling-------------------------------
% with the baseline FA [3,6,9,12,15] degrees, it
% will be derived from the signal acquired at FA 3 and 6 degrees:
% Method to estimate T1 from 3,6 deg Signal:
%Option 1: Helms et al. 2008:?Quantitative FLASH MRI at 3T Using a Rational Approximation of the Ernst Equation
% T1_FA36 = 2*TR*((DataPerSlicePerPixel(1)/faTrueRad(1)-DataPerSlicePerPixel(2)/faTrueRad(2))/(DataPerSlicePerPixel(2)*faTrueRad(2)-DataPerSlicePerPixel(1)*faTrueRad(1)));
%if index_FASecondLowest is above the Ersnt Angle we can't calculate the T1
%estimate from this approximation as it's only valid for small angles, but we can use the ratio (Option 2)
%check the second FA is below the Ernst Angle
% MeanT1UpperLowerBound = (ubT1+lbT1)/2; %average between ubT1 and lbT1
% ErnstAngleMeanT1 = acosd(exp(-TR/MeanT1UpperLowerBound));  
% Option2: no approximations, worked out the
% formula from SPGR function by taking the ratio of two signals

NumE1 = (Data_FALowest*sin(faTrueRad(index_FASecondLowest))-Data_FASecondLowest*sin(faTrueRad(indexFALowest)));
DenE1 = (-Data_FASecondLowest*sin(faTrueRad(indexFALowest))*cos(faTrueRad(index_FASecondLowest))+Data_FALowest*sin(faTrueRad(index_FASecondLowest))*cos(faTrueRad(indexFALowest)));
E1 = NumE1/DenE1;
E1Max = exp(-TR/ubT1);
E1Min = exp(-TR/lbT1);

if (E1>E1Max) || (E1<E1Min)
    T1_InitialEst = MeanT1UpperLowerBound;
else
    T1_InitialEst = -TR/log(E1);
end

%---------------------------------------------------------------------------------------------

end

