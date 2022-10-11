function [matrixData] = cell2matrix(cellData,nr,nc,nslices,nFAs)

%cell2matrix transforms a cell variable into a matrix
%     Inputs:
%         cellData - name of cell
%         nr - number of rows
%         nc - number of columns
%         nslices - number of slices
%         nFAs - number of FAs acquired
%
%     Outputs:
%         matrixData: matrix
%
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

matrixData = zeros(nr,nc,nslices,nFAs);
for ifa = 1:nFAs
    for islice = 1:nslices
        matrixData (:,:,islice,ifa) = cellData{islice,ifa}(:,:);
    end
end
end

