function [ratioInten_te] = Mxy_te_gradientB0SliceProfile(B1Bias,nominalAngle,kFactor,df,T1,T2,TE,sliceposition,SliceThicknessFactor,FatSup)

%Mxy_te_gradientB0SliceProfile
%Summary: calculate the Signal130/Sginal65 assuming varying B0 across the slice profile
%       Input: 
%         B1Bias - unkown parameter to estimate
%         nominalAngle - alpha (65 degrees)
%         kfactor - k*alpha (2*alpha)
%         df - off-resonance for each z position used in the Bloch simulator
%         T1, T2 - relaxation times
%         TE: echo time
%         sliceposition - simulate mxy between -1 and 1 cm
%         SliceThicknessFactor: original protocol had 8mm-SliceThicknessFactor=1, but we later also experimented with 4mm-SliceThicknessFactor=2
%         FatSup - Fat Suppression method: Fat Sat (FS) or Water excitation (WE)

%       Output: 
%         ratioInten_te: ratio of integrated signals Signal2alpha/Sginalalpha

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


%path to Bloch simulation: Wind_blochMxyRFInterp_gradientB0SliceProfile
addpath('/Users/gabrielabelsley/OneDrive - Perspectum Ltd/DPhil_SecondProject/MRIPhysics/UCLACourses/AdvancedMRI/Lesson2_RFPulseDesign/matlab_demo/')



%============Mxy at end of the RF pulse============

nsamples_sliceProfile = length(sliceposition);
mx = ones(1,nsamples_sliceProfile);
mx_2alpha = ones(1,nsamples_sliceProfile);
my = ones(1,nsamples_sliceProfile);
my_2alpha = ones(1,nsamples_sliceProfile);
mz = ones(1,nsamples_sliceProfile);
mz_2alpha = ones(1,nsamples_sliceProfile);

%Wind_blochMxyRFInterp_gradientB0SliceProfile_Gzref uses new greEPI RF
%pulse sent by Liz in Nov 2020 using the new protocol with 11ms TE and
%grappa turned off: Liz used the parameters from 3rd sequence in protocol
%3DB1T1Protocol_inVivo_21102020.pdf; Fat Sat EPI pulse saved in 
%Meeting10112020/ConcomitantB0Field/greEPI_SimulationProtocol_RFD/SimulationProtocol_GRZ.dsv',tmax);
for ipos = 1:length(sliceposition)
    zpos = sliceposition(1,ipos);
    freq = df(1,ipos);
    if strcmp(FatSup(1:2),'FS')
        [~,mx(1,ipos),my(1,ipos),mz(1,ipos)] = Wind_blochMxyRFInterp_gradientB0SliceProfile_Gzref(B1Bias,nominalAngle,freq,zpos,SliceThicknessFactor,T1,T2);
        [~,mx_2alpha(1,ipos),my_2alpha(1,ipos),mz_2alpha(1,ipos)] = Wind_blochMxyRFInterp_gradientB0SliceProfile_Gzref(B1Bias,nominalAngle*kFactor,freq,zpos,SliceThicknessFactor,T1,T2);
    elseif strcmp(FatSup(1:2),'WE')
        [~,mx(1,ipos),my(1,ipos),mz(1,ipos)] = WE_blochMxyRFInterpGradShift(B1Bias,nominalAngle,freq,zpos,T1,T2);
        [~,mx_2alpha(1,ipos),my_2alpha(1,ipos),mz_2alpha(1,ipos)] = WE_blochMxyRFInterpGradShift(B1Bias,nominalAngle*kFactor,freq,zpos,T1,T2);

    end
end



%============ TE defined from the centre of the RF pulse and not the end of the Gz ref gradient==========
%get the time (microseconds) where the peak of the RF pulse occurs
if strcmp(FatSup(1:2),'FS')
    Time_CentreRFUntilGzRefEnd = (1820)*10^-3; %Gzref = 22.8mT/m
    %Time_CentreRFUntilGzRefEnd = (2443)*10^-3; %Gz ref = 5.96mT/m
    freePrecessionTime = TE -Time_CentreRFUntilGzRefEnd;% TE is defined from the centre of the RF pulse, but free precession only starts after the Gz Ref ends
elseif strcmp(FatSup(1:2),'WE')
    Time_CentreRFUntilGzRefEnd = (2887.5)*10^-3; %ms
    freePrecessionTime = TE -Time_CentreRFUntilGzRefEnd;% TE is defined from the centre of the RF pulse, but free precession only starts after the Gz Ref ends
end
%============ Calulate the Mxy for each slice position with T1 relaxation and T2 decay after time TE ============ 
%It uses the signal calculated at the end of the RF
%pulse. Each slice position has a certain off-resonance frequency saved in
%df which will result in a different propagation matrix

Mx_te = zeros(nsamples_sliceProfile,1);
My_te = zeros(nsamples_sliceProfile,1);
Mz_te = zeros(nsamples_sliceProfile,1);
Mx_te_2alpha = zeros(nsamples_sliceProfile,1);
My_te_2alpha = zeros(nsamples_sliceProfile,1);
Mz_te_2alpha = zeros(nsamples_sliceProfile,1);

for islicePos = 1:nsamples_sliceProfile %loop over the slice profile
    
   

M_alpha = [mx(1,islicePos);my(1,islicePos);mz(1,islicePos)]; %Magnetization immediately after you play excitation RF pulse alpha
M_2alpha = [mx_2alpha(1,islicePos);my_2alpha(1,islicePos);mz_2alpha(1,islicePos)]; %Magnetization immediately after you play excitation RF pulse 2*alpha



% ===== Get the Propagation Matrix ======
%propagates the signal from end of RF pulse until time TE
[Ate,Bte] = freeprecess(freePrecessionTime,T1,T2,df(1,islicePos));


Mte = Ate*M_alpha+Bte;	% Magnetization at t=TE for alpha
Mte_2alpha = Ate*M_2alpha+Bte;	% Magnetization at t=TE for 2*alpha


Mx_te(islicePos,1) = Mte(1,1);
My_te(islicePos,1) = Mte(2,1);
Mz_te(islicePos,1) = Mte(3,1);

Mx_te_2alpha(islicePos,1) = Mte_2alpha(1,1);%Mx is the first component
My_te_2alpha(islicePos,1) = Mte_2alpha(2,1);%My is the second component
Mz_te_2alpha(islicePos,1) = Mte_2alpha(3,1);%Mz is the third component
end

mxy_te = Mx_te+1i.*My_te; %mxy for alpha
mxy_te_2alpha = Mx_te_2alpha+1i.*My_te_2alpha; %mxy for 2*alpha

[ratioInten_te] = ratio_kAlpha_Alpha(mxy_te,mxy_te_2alpha);

end

