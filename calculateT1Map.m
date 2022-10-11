function [m0MapNLLS,b1MapNLLS,t1MapNLLS,rsqMapNLLS,ciT1MapNLLS,NormSSE,deltaT1,Norm2Residuals,ResidualsFit,CorrFactorSpoil,lbT1,ubT1] = calculateT1Map(vfaHeader,vfaData,nFA_VFA,maxB1bias,B1MapMatchVFA,flagInCompSpoilCorr,slicesVFA_T1Map)

%CALCULATET1Map T1 map calculation using vfa data, B1+ map and incomplete
%spoiling correction

%     Inputs:
%       vfaHeader: VFA Header
%       vfaData: VFA Data
%       nFA_VFA: number of SPGR FAs
%       maxB1bias: maximum B1+ factor value
%       B1MapMatchVFA: B1+ map at the VFA resolution
%       nslices_VFA: number of slices in the VFA acquisition
%       flagInCompSpoilCorr: boolean flag, 1 to include incomplete spoiling correction

%     Outputs: All outputs at the VFA resolution
%       m0MapNLLS: M0 Map
%       b1MapNLLS: B1+ Map
%       t1MapNLLS: T1 Map
%       rsqMapNLLS: R^2 Map
%       ciT1MapNLLS: T1 confidence Interval
%       NormSSE: Normalized Sum Squared Errors: Sum Squared Error for each pixel normalized by the mean signal acquired for the nFA_VFA FAs
%       deltaT1: difference in the upper and lower bound of the T1 Confidence interval
%       Norm2Residuals: Squared norm of the residual
%       ResidualsFit: Residuals
%       CorrFactorSpoil: Incomplete Spoiling correction factor
%       lbT1,ubT1: lower and upper T1 bound

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


%============ Step 1: Read FAs, TR, number of slices ===================

faVFADeg = zeros(1,nFA_VFA);
for iFA = 1:nFA_VFA
    faVFADeg(1,iFA) = vfaHeader{1,iFA}.FlipAngle;
end
faVFARad=(faVFADeg.*pi)/180;

TR = vfaHeader{1,1}.RepetitionTime;

%============ Step 2: Incomplete Spoiling Correction of SPGR signal===================

%----072020: Spoiling Correction using the interpolation of surface Correction (which is a function of T1 and alpha)------------------
%tumours are inflamed and can have large T1s (>1500ms); for the healthy
%volunteer study I used smaller T1 bounds between 500ms and 1600ms
load('SpoilCorr/EPGSignalSS_FA1_120_T1500_3000_TR4.1_T230.mat','T1','alphaArray','sigSSNAlphaTheory','sigSSNAlphaEPG')
%Note: FDAPhantom Perspectum, T1mesPhantom, Calimetrix QC0025
%Need different Spoil Correction as T1 value range and T2 differs from in vivo
%T1mesPhantom 
%load('SpoilCorr/EPGSignalSS_FA1-20_T1100-2000_TR4.1_T246.mat','T1','alphaArray','sigSSNAlphaTheory','sigSSNAlphaEPG')

NT1 = length(T1);
Nangles = length(alphaArray); 
CorTheoryEPG50deg = zeros(NT1,Nangles);
for iT1 = 1:NT1
CorTheoryEPG50deg(iT1,:) = (sigSSNAlphaTheory(1,iT1,:)./sigSSNAlphaEPG(1,iT1,:));
end
[T1Grid,AlphaGrid] = ndgrid(T1,alphaArray);
Correction = CorTheoryEPG50deg;
%figure()
%surfc(T1Grid,AlphaGrid,Correction)
SpoilingCorrection = griddedInterpolant(T1Grid,AlphaGrid,Correction,'cubic');

lbT1 = min(T1);%ms
ubT1 = max(T1);%ms

if exist('flagInCompSpoilCorr','var')&& flagInCompSpoilCorr== 0 %No Incomplete Spoiling Correction
    SpoilingCorrection = [];
    lbT1 = 0;%min(T1);%ms
    ubT1 = 2000;%max(T1);%ms
end

%============ Step 3: Calculate T1 ===================
[m0MapNLLS,b1MapNLLS,t1MapNLLS,rsqMapNLLS,ciT1MapNLLS,NormSSE,deltaT1,Norm2Residuals,ResidualsFit,CorrFactorSpoil] = nllsfitM0B1T1Map_x0ublbSpoil_CorrSigAcq_GoodInitGuess(faVFARad,vfaData,slicesVFA_T1Map,TR,maxB1bias,B1MapMatchVFA,lbT1,ubT1,SpoilingCorrection);

%========== End of script %========== 
end

