function LiverMask_NoHoles = removeHolesLiverMaskFromB1Map(B1Map_B0GradZCorr,nr_greEPI,nc_greEPI,nslices_greEPI)

%removeHolesLiverMaskFromB1Map remove holes from the mask defined where the B1+ map was calculated, i.e.
%where the B1+ map is different from NaN

%07042021
% Rationale for this function: Because many pixels had B0=0Hz in one of the
% neighbour or target slices, the B1+ (output of calculateB1MapB0GradZCorr)
% was not calculated and thus set to NaN. If we use the
% maskForB1Map_gradientB0Z in removeFatBandB1 the removeFatBandB1 will
% extrapolate the B1+ to all pixels where the B1+ was NaN but that was 1 in
% the mask. It extrapolates the B1+ based on the curvature of the B1+ map.
% Solution: constructed a mask from the B1+ (output of
% calculateB1MapB0GradZCorr) map pixels which are different from NaN (line
% LiverMask = ~isnan(B1Map_B0GradZCorr{islice,1})) and filled in the holes
% first in the vertical direction (for each column) and then in the
% horizontal direction (for each row)
%The output mask of this function is the input to removeFatBandB1. 

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

LiverMask_NoHoles = zeros(nr_greEPI,nc_greEPI,nslices_greEPI);
for  islice = 2:nslices_greEPI-1
    LiverMask = ~isnan(B1Map_B0GradZCorr{islice,1});
    if ~all(LiverMask(:)==0)
        [~,col2Fill] =find(LiverMask==1);
        colStart = col2Fill(1,1);
        colEnd = col2Fill(end,1);
        for icol = colStart:colEnd %fill the rows (vertical lines) for each column
            [rowMask] = find(LiverMask(:,icol)==1);
            if ~isempty(rowMask)
                LiverMask_NoHoles(rowMask(1,1):rowMask(end,1),icol,islice) = 1;
            end
        end
        [row2Fill,~] =find(LiverMask_NoHoles(:,:,islice)==1);
        rowStart = min(row2Fill(:));
        rowEnd = max(row2Fill(:));
        for irow = rowStart:rowEnd %fill the columns (horizontal lines) for each row
            [~,colMask] = find(LiverMask_NoHoles(irow,:,islice)==1);
            if ~isempty(colMask)
                LiverMask_NoHoles(irow,colMask(1,1):colMask(1,end),islice) = 1;
            end
        end
    end
end

end

