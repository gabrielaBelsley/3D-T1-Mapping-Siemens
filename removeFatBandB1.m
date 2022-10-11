function [B1AllCorrections_NoFat] = removeFatBandB1(B1Map_B0GradZCorr,nslices_greEPI,nr_greEPI,nc_greEPI,maskForB1Map_gradientB0Z)

%REMOVEFATBANDB1 remove fat chemical shift from B1+ Map 

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

%============ Remove chemical shift fat from the B1+ map ===================
B1AllCorrections_NoFat = zeros(nr_greEPI,nc_greEPI,nslices_greEPI);
for islice = 2:nslices_greEPI-1
    if any(~isnan(B1Map_B0GradZCorr{islice,1}(:)))
    [B1NoFat] = B1Map_RemoveFatBand2(B1Map_B0GradZCorr{islice,1},maskForB1Map_gradientB0Z(:,:,islice));
    B1AllCorrections_NoFat(:,:,islice)=B1NoFat;
    end
end

end

