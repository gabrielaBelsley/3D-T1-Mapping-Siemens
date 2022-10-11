function [ExtrmLogical,PkProm,ValProm] = IdExtrema(Img, Thres)

% IdExtrema Identifies localized extrema in an image Img
% Any point that is more than Thres (above or below) the average of the 8
% surrounding pixels is marked (logical 1) in the output ExtrmLogical 

% input: Img - image
%        Threshold - limit of absDifference with average of sourrounding
%        area that marks an isolated extrema

% output: ExtrmLogical - logical matrix with same dimensions as Img
%                   pixels set to one identify extrema
%         PkProm - matrix listing prominent Peaks
%               each entry (row) lists: Row#, Col#, Peak Height
%         ValProm - matrix listing prominent Valleys
%               each entry (row) list: Row#, Col#, Valley Depth
% external function extrema2

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


%First use extrema2 to identify all extrema: peaks and valleys of B1+ surface
[Pks,PksIndx,Vals,ValsIndx] = extrema2(Img);

% calculate the "prominance" difference from average value of surrounding
% pixels for each of the extrema found

%Eliminate Peak B1+
PkProm= zeros(length(Pks),3);  % 1st column row #, 2nd column col #, 3rd column prominance
cnt = 1;
for j = 1:length(Pks)
    [r,c]  = ind2sub(size(Img),PksIndx(j));
    NoNaNs = sum(isnan(Img(r-1:r+1,c-1:c+1)),'all');
    SurSum = sum(Img(r-1:r+1,c-1:c+1),'all','omitnan')-Img(r,c);%sum of the B1+ of the 8 surrounding pixels - B1+peak
    Prom = Img(r,c)-SurSum/(8-NoNaNs);%subtract the peakB1+ from the average B1+ of the surrounding 8 pixels (=SurSum/(8-NoNaNs))
    
    % ﻿If the absolute difference between the outlier B1+ factor and the
    % average B1+ factor of the surrounding 8 pixels is larger than 0.075, the
    % B1+ pixel is considered a true outlier.
    if (Prom >= Thres) %If the Height Peak is larger than 0.075 (0.075 error in B1+ would give ~0.15 (15%)error in T1)
        PkProm(cnt,1) = r;
        PkProm(cnt,2) = c;
        PkProm(cnt,3) = Prom;
        cnt = cnt+1;
    end
end
PkProm(cnt:length(Pks),:,:) = [];

%Eliminate Valley B1+

ValProm= zeros(length(Vals),3);
cnt = 1;
for j = 1:length(Vals)
    [r,c]  = ind2sub(size(Img),ValsIndx(j));
    NoNaNs = sum(isnan(Img(r-1:r+1,c-1:c+1)),'all');
    SurSum = sum(Img(r-1:r+1,c-1:c+1),'all','omitnan')-Img(r,c);
    Depth = SurSum/(8-NoNaNs)-Img(r,c);
    if (Depth >= Thres)
        ValProm(cnt,1) = r;
        ValProm(cnt,2) = c;
        ValProm(cnt,3) = Depth;
        cnt = cnt+1;
    end
end
ValProm(cnt:length(Vals),:,:) = [];


% Set these isolated peaks and valleys to 1 

ExtremaIdentified = zeros(size(Img));
for j= 1:size(PkProm,1)
    ExtremaIdentified(PkProm(j,1),PkProm(j,2)) = 1;
end
for j= 1:size(ValProm,1)
    ExtremaIdentified(ValProm(j,1),ValProm(j,2)) = 1;
end
ExtrmLogical= logical(ExtremaIdentified);

end % Extrema


