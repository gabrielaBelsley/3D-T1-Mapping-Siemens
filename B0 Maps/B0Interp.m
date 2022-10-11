function B0Fit = B0Interp(B0Map,Mask, SliceNo,varargin)

%Function to calculate the interpolated B0 map along the slice direction
%  Uses John D'ericco's slm interpolator across the whole masked region
% uses a linear interpolation to extend the B0 map into a region outside the
% mask
%   Input:
%        B0Map: B0 map values (matched to B1 grid)
%        Mask: logical vector of zeros and ones that describe the mask along the slice direction for the chosen pixel
%        SliceNo: Slice number around which we want to interpolate the B0 map (output will extend from previous slice to next slice)

%   Output       
%       B0Fit: Interpolated B0 map from SliceNo-1 to SliceNo+1
%   
%   Internal parameters Step: spacing for points to be interpolated

%   Note: requires John D'ericco's slm program

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

Step = 0.1;%cm
flagPhantom = varargin{1,1};
NSlices = length(Mask);
Slices = 1:NSlices;
% find the start and end slices of the mask
StartSlice = find(Mask~=0,1);
EndSlice = find(Mask ~=0,1,'last');
if flagPhantom==1
    StartSlice = SliceNo-1;
    EndSlice = SliceNo+1;
end
%removing B0 pixels outside the mask
B0Map([1:StartSlice-1 EndSlice+1:NSlices])=[];
Slices([1:StartSlice-1 EndSlice+1:NSlices])=[];
TotalSlices = EndSlice-StartSlice+1;
Nknots = max(2,floor(TotalSlices/3));
%spline fit between knots: optimal at least 3 slices between the knots
%extrapolation of B0 with linear fit to the slices where B0 was not calculated
%fit B0 to all slices 
slm = slmengine(Slices,B0Map,'knots',Nknots,'extrapolation', 'linear');
%evaluate the B0 at the neighbour slices
SlicePos = SliceNo-1:Step:SliceNo+1;
B0Fit = slmeval(SlicePos,slm);
%Check if B0 Fit matches B0 Map values
% figure()
% plot([7:9],B0Map(4:6))
% hold on
% plot(SlicePos,B0Fit)
end

