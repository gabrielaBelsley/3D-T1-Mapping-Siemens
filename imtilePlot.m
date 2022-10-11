function imtilePlot(Map,slices2Plot,mapTitle,legendColourbar,ColorScheme,colourmapScale,nRows,nCols,flagColorbar,ftsz)

%IMTILEPLOT plot T1/B1+ map using imtile function
%
%     Inputs:
%       Map: Map to plot
%       slices2Plot: slices to plot
%       mapTitle: Title of map
%       legendColourbar: Colorbar legend
%       ColorScheme: choose a colormap: 'parula', 'jet', etc
%       colourmapScale: scale of colormap
%       nRows,nCols: matrix to show the map, e.g. for 15 slices use nRows,nCols=3,5
%       flagColorbar: boolean, show or omit colorbar
%       ftsz: fontsize of colorbar and legends
%
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


%Date 09/02/2021

%Plot Map using imtile
out = imtile(Map, 'Frames',slices2Plot, 'GridSize', [nRows nCols]);
%figure('Position', [150 100 2000 2000])
figure('Color',[1 1 1]) % change the background to white
ax = axes;
imshow(out,colourmapScale);
 % originalSize = get(gca, 'Position');
cmap =colormap(ax,ColorScheme);
 cmap(1, :) = [0 0 0]; %for black background instead of blue
 colorbar
colormap(cmap);
if flagColorbar
    cbar = colorbar(ax);
    cbar.Location='southoutside';
    cbar.FontSize = ftsz;
    cbar.FontName = 'Arial';
    cbar.Color = [0 0 0];%[1 1 1] white
    cbar.Label.String = legendColourbar;
    cbar.Label.Color = [0 0 0];%[1 1 1] white
    cbar.Label.FontName = 'Arial';
    cbar.Label.FontSize = ftsz;   
    cbar.XTick = colourmapScale(1):200:colourmapScale(2);
    % cbar.XTick = 0:50:250;
%    originalSize(3) = originalSize(3)-0.025;
%      set(ax, 'Position', originalSize);
end


%Add a), b) or c) to figure
% xt = -0.01;
% yt = 1.05;
% str = {'(a)'};
% text(xt,yt,str,'Units','normalized','fontsize',ftsz,'FontName','Arial','FontWeight','bold')
title(ax,mapTitle,'FontName','Arial','FontSize',ftsz)



end

