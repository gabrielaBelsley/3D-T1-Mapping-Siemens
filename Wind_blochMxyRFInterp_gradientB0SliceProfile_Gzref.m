function [mxy,mx,my,mz] = Wind_blochMxyRFInterp_gradientB0SliceProfile_Gzref(b1_bias,angle,f,sliceposition,SliceThicknessFactor,T1ms,T2ms)

% Wind_blochMxyRFInterp_gradientB0SliceProfile_Gzref
%Bloch simulator using the GRE EPI Hamming Windowed Sinc RF excitation pulse
%   Input:
%       b1_bias: B1+ Factor
%       angle: nominal excitation angle
%       f: off-resonance frequency [Hz]
%       SliceThicknessFactor: 1=slice thickness 8mm; 2=slice thickness 4mm
%       T1ms,T2ms: relxation times in ms
%   Output:
%       Transverse and Longitudinal Magnetization

%dependencies: uses bloch function written by Brian Hargreaves

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

%Gradient Z from simulations, with ramps
load('EPI_FID_SliceProfile/epi_fid_rf_grz.mat','rf_phase_interp','gz_interp','time_interp')
if ~exist('SliceThicknessFactor','var')
    SliceThicknessFactor = 1;
end
gz_SlThick = SliceThicknessFactor.*gz_interp;
g   = gz_SlThick./10; % mT/m to G/cm


%RF pulse
%time bandwidth: product of its temporal duration and spectral width (in frequency space)
%we normalize the amplitude from the simulator tool to be able to then
%multiply by the different nominal angles magnitude. 
%FA = integral (gamma*B1dt), sum(rf_phase_interp) is essentially this integral. 
rf_phase_interp = rf_phase_interp./sum(rf_phase_interp);    % normalized so that sum(h) = 1
%Because we normalized, if now we take the integral it will give the FA as FA = integral (gamma*B1dt)
rf_phase_interp = (b1_bias.*angle).*rf_phase_interp; % rf waveform, scaled so sum(rf) = flip angle
gamma = 2.*pi.*4257.6; % Hz/G
dt = (time_interp(2)-time_interp(1)).*1e-6; % in s
rf_phase_interp = rf_phase_interp./(gamma.*dt); % in Gauss
b1  = rf_phase_interp;% in Gauss

%time vector with same length as RF pulse
t = repmat(dt,size(b1));  % in s
%spatial position RF pulse
x = sliceposition;
T1 = T1ms*1e-3;%0.5; %10000%s 1
T2 = T2ms*1e-3;%0.04;%10000;%s 0.2


% Bloch Simulation
[mx,my,mz] = bloch(b1,g,t,T1,T2,f,x,0); 


% Transverse Magnetization
mxy = mx+1i.*my;
%mxy_onresonance = abs(mxy(:,round(end/2)));
end

