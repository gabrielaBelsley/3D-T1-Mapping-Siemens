function [mxy,mx,my,mz] = WE_blochMxyRFInterpGradShift(b1_bias,angle,f,sliceposition,T1ms,T2ms,varargin)

% WE_blochMxyRFInterpGradShift
%Bloch simulator using the GRE EPI Water Excitation (WE) pulse
% WE Pulse: 121 binomial RF pulse played with a spatially selective gradient

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


load('EPI_FID_SliceProfile/epi_WE_pulse.mat', 'rf', 'time')
load('EPI_FID_SliceProfile/WE_pulse_grad.mat', 'tGrad', 'yGrad')


%binominal 121
rf_phase = rf;
%samples = length(time);

%The pulse comes out of the simulator all positive, need to negate the side
%lobes for each of the 3 truncated sinc pulses forming the binomial WE pulse. 
%windowed sinc 1
rf_phase(248:289)=-rf(248:289); 
rf_phase(374:415)=-rf(374:415); 

%windowed sinc 2
rf_phase(494:535)=-rf(494:535); 
rf_phase(620:661)=-rf(620:661); 

%windowed sinc 3
rf_phase(740:781)=-rf(740:781); 
rf_phase(866:907)=-rf(866:907);

time_interp = 1:time(end);
rf_phase_interp = interp1(time, rf_phase,time_interp,'pchip');


% Intepolation of the gradient waveform so that step size in time is equal for the RF pulse and the gradient:
% Interpolate the gradient signal to have values in steps of 5, to match the
% steps in the rf_phase pulse, which is also 5. (the original step size of the gradient is 10). 
x = tGrad;  
v = yGrad;                                      % function to interpolate
xq = 1:(time_interp(2)-time_interp(1)):tGrad(525);   
yGradInterp = interp1(x,v,xq,'pchip');

%%%%%%%%%%%%%%%
% Calculating again the gradient zeroth moment, the peak of the RF did not
% align with the zeroth moment gradient, thus shifted by 3.
shift = -3;
yGradInterpShift = zeros(1,length(yGradInterp));
start = 1+abs(shift);
y=start:length(yGradInterp);
yGradInterpShift(y+shift) = yGradInterp(y);
%tested with shift of zero and it gives equal B1 Maps with B0-Z correction
%to using shift of -3

% if shift>0
% yGradInterpShift = [yGradInterpShift zeros(1,abs(shift))];
% else
%     yGradInterpShift = [zeros(1,abs(shift)) yGradInterpShift ];
% end

extrazeros = length(yGradInterpShift)-length(rf_phase_interp);
rf_phase_padded= [rf_phase_interp zeros(1,extrazeros)];
%%%%%%%%%%%


%time bandwidth: product of its temporal duration and spectral width (in frequency space)
rf_phase_padded = rf_phase_padded/sum(rf_phase_padded);    % normalized so that sum(h) = 1
rf_phase_padded = (b1_bias.*angle).*rf_phase_padded; % rf waveform, scaled so sum(rf) = flip angle
gamma = 2*pi*4257.6; % Hz/G
dt = (time_interp(2)-time_interp(1)).*1e-6; % in s
rf_phase_padded = rf_phase_padded./(gamma.*dt); % in Gauss
b1  = rf_phase_padded;              % in Gauss
g   = yGradInterpShift./10;    % zeros(1,samples*1.5);%mT/m to G/cm 10-3/(10-4*10^2)=10^-1


% gradient with flyback
% g(701) = -sum(g(1:455));
% g(456:700) = g(456:700)* -1;
% g(455) = -sum(g(1:455));


% Bloch Simulation
%Commented 03/11/2021
% if nargin == 5
%     f = varargin{1};
%     [mx,my,mz] = bloch(b1,g,t,T1,T2,f,x,0);
% elseif nargin ==6
%     f = varargin{1};
%     x = varargin{2};
%     [mx,my,mz] = bloch(b1,g,t,T1,T2,f,x,0);
% else
% 
%     [mx,my,mz] = bloch(b1,g,t,T1,T2,0,x,0);
% end

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

