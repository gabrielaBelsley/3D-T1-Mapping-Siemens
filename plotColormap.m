    function plotColormap(map,mapTitle,scaleLim,ColorScheme,unitsMap,zoomFactor,flagBlackBackground,ax,flagColorbar,textLabel,scaleUnit,flagYXlabel,varargin)

    % plotColormap(map,mapTitle,scaleLim,varargin)
    % colormapPlot function that plots the colormap in the map variable with a
    % colorbar scale with limits defined within scaleLim

   %   Input:
    %       map                 map to plot with lower resolution
    %       mapTitle            title
    %       scaleLim            [min max] scale
    %       ColorScheme         colormap scheme 'parula', 'grey'
    %       unitsMap            Map units, e.g. T1 in ms
    %       zoomFactor          factor to zoom in the map
    %       flagBlackBackground set to 1 for black blackground
    %       ax                  figure axis
    %       flagColorbar        set to 1 to  include colorbar
    %       textLabel           label (a), (b)
    %       scaleUnit           steps in colormap, e.g. 100ms for T1
    %       flagYXlabel         label for x, y axis
    

    %   Author: Gabriela Belsley, 17/01/2019
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

if ~exist('ax','var')
    ax=gca;
end

if ~isempty(scaleLim)
    % Problem with imagesc: it distorts the image to fit in the figure
    % window: e.g. if image dimensions are not equal in x,y: 104*128, then
    % it will stretch the x dimension vertically such that the image has
    % same dimension in the window for the x(vertical=row) and y(horizontal=col) directions.

    imagesc(squeeze(map),scaleLim)
    if exist('ColorScheme','var')
        colormap(ax,ColorScheme)
    end
else
    imagesc(squeeze(map))
    % imshow(squeeze(map),[],'InitialMagnification','fit')
    if exist('ColorScheme','var')
        colormap(ax,ColorScheme)
    end
end

ftsz=14;
    
if ~(isempty(mapTitle))

    title(ax,mapTitle,'FontName','Arial','FontSize',ftsz)
end


%axis off
ax.FontSize = ftsz;
ax.FontName = 'Arial';
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
%ax.XTick=140:5:160;
%ax.YTick=125:15:180;
axis image
%posfig = get(gca, 'Position'); % gives x left, y bottom, width, height
% originalSize = get(gca, 'Position');

if exist('flagColorbar','var') && flagColorbar ==1
    c = colorbar(ax);
    c.FontName = 'Arial';
    c.FontSize = ftsz;
    c.Location = 'southoutside';
  

    %c.Label.Interpreter = 'latex';
    c.Label.FontName = 'Arial';
    c.Label.FontSize = ftsz;
%  set(ax, 'Position', originalSize);
    if exist('scaleLim','var')
        if exist('unitsMap','var')
            c.Label.String = unitsMap;
            if ~isempty(scaleLim)
                %range = scaleLim(1):(round((scaleLim(2)-scaleLim(1))/,2)):scaleLim(2);
                range = scaleLim(1):scaleUnit:scaleLim(2);
                set(c,'ytick',range);
%                 cbarPos =  c.Position; %gets the positon and size of the color bar
%                 set(c, 'Position',  [cbarPos(1), cbarPos(2), cbarPos(3), posfig(4)])

            end
        elseif (scaleLim(1,1)>90 && scaleLim(1,2)>750)
            c.Label.String = 'T_1 (ms)';
            range = scaleLim(1):100:scaleLim(2);
            set(c,'ytick',range);
        elseif (scaleLim(1,2)>1 && scaleLim(1,2)<2.1)
            c.Label.String = 'B_{1}+ Factor';
            range = scaleLim(1):0.3:scaleLim(2);
            set(c,'ytick',range);
            
        else
            c.Label.String = '';
        end
    end
end


% %Add a), b) or c) to figure
if exist('textLabel','var')
xt = -0.01;
yt = 1.04;
str = {textLabel};
text(xt,yt,str,'Units','normalized','fontsize',ftsz,'FontName','Arial','FontWeight','bold')
end
%impixelinfo
%axis off
if exist('flagYXlabel','var')
    ylabel(flagYXlabel,'FontName','Arial','FontSize',ftsz)
    axis on
    %ax.XAxis.Visible = 'off';
    [nrows,~] = size(map);
    ax.YTick = 1:3:nrows;
    ax.XTick = [ ];
    ax.FontSize = ftsz;
    %ax.XTick = 17:6:nrows;
    % xlabel('Phase-Encoding Row Number','FontName','TimesNewRoman','FontSize',12)
end
%  axis image:aspect ratio of the display will match the aspect ratio of
%  the image data matrix.
%  axis image %: does not show x axis ticks, need axis on above
 
 

if exist('zoomFactor','var')
    zoom(zoomFactor)
end


if exist('flagBlackBackground','var') && flagBlackBackground ==1
    cmap = colormap(ax);
    cmap(1, :) = [0 0 0]; %for black
    colormap(ax,cmap);
end

%alternative way of setting background (=all NaN) to black
%         imagesc(squeeze(map),'AlphaData',double(~isnan(map)),scaleLim)%, it changes the aspect ratio
%        set(gca,'color',0*[1 1 1]);



    %     if strcmp(mapTitle,'B0')
    %     fprintf('The maximum off-resonance is %f Hz \n',max(map(:)));
    %     fprintf('The minimum off-resonance is %f Hz \n',min(map(:)));
    %     if strcmp(mapTitle,'B1')
    %         fprintf('The maximum B1 correction is %f\n',max(map(:)));
    %         fprintf('The minimum B1 correction is %f\n',min(map(:)));
    %     end
    %     https://uk.mathworks.com/help/matlab/ref/fprintf.html



    % code used in B1MapScript_InVivoPhantom.m to evaluate the intensity at a
    % pixel input by the user: goal: check that the intensity does not hit the
    % maximum value = 2^12-1
    if ~isempty(varargin)
        hold on
        figureHandle=varargin{1};
        axMap = varargin{2};
        map2 = varargin{3};

        cnt=0;
        pixel = [4,5];
        while ~isempty(pixel)%ishghandle(figureHandle)
            cnt=cnt+1;
            pixel = input(['Pixel([x,y]) ',num2str(cnt),': ']);
            if isempty(pixel)
                break
            else
                fprintf([mapTitle,': %d \n'],map(pixel(2),pixel(1)));
                fprintf(1, '\n');
                if ~isempty(map2)
                    mapTitle2= varargin{4};
                    fprintf([mapTitle2,': %d \n'],map2(pixel(2),pixel(1)));
                    fprintf(1, '\n');
                end
                figure(figureHandle)
                plot(axMap,pixel(1),pixel(2),'o', 'MarkerFaceColor', 'r')
                text(pixel(1),pixel(2),[' ', num2str(cnt),': ',num2str(map(pixel(2),pixel(1)),3)], 'Color','r')
                % str = {['Pixel:',num2str(Row),',',num2str(Col)];'Value: '};
                % annotation('textbox',dim,'String',str,'FitBoxToText','on');
                %         t = text(1,1,{['Pixel:',num2str(Row),',',num2str(Col)],'Value: '})
                %         drawnow();
                %         delete(t)
            end
        end
        hold off
    end



