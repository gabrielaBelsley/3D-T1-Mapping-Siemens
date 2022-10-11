function [B1MapNoExtrema] = EliminateIsolatedExtrema(B1Map, LiverMask, Threshold, verbose)

% EliminateIsolatedExtrema
% The reasoning here is that the B1+ Map should be relatively smooth
% So that isolated points that are more than a certain Threshold (0.075)
% above the average of their surroundings should not exist. 
% They are identified, set to NaN and then filled in (painted) using
% interpolation.
% Calls the function IdExtrema to identify the extrema
% The routine iterates twice

%   Input variables
%          B1Map: B1+ map to clean up
%          LiverMask: The logical mask that defines the liver area
%          Threshold: value above or below the average B1+ value of the
%                     surrounding pixels which identifies this as an extrema
%          verbose: Logical flag if = 1, graphs are output

%    Output: cleaned up B1+ Map without extrema/outliers

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021


[ExtB1Map,Pks,Vals] = IdExtrema(B1Map, Threshold);

% if verbose
%     figure()
%     h = imagesc(B1Map,[0.5 1.5]);
%     set(h, 'AlphaData', ~isnan(B1Map))
%     set(gca, 'Color', [0, 0, 0])
%     colorbar
%     hold on
%     for j = 1:length(Pks)
%         scatter(Pks(j,2),Pks(j,1),[],'m','^','filled')
%     end
%     for j= 1:length(Vals)
%         scatter(Vals(j,2),Vals(j,1),[],'m','v', 'filled')
%     end
%     hold off
%     title('B1Map Extrema 1st Iteration')
% end

B1MapNoPks = B1Map;
B1MapNoPks(ExtB1Map) = NaN;

if verbose
    figure()
    h = imagesc(B1MapNoPks,[0.5 1.5]);
    set(h, 'AlphaData', ~isnan(B1MapNoPks))
    set(gca, 'Color', [0, 0, 0])
    colorbar
    title('B1Map Extrema -> NaN 1st Iteration')
end

% Extrapolate to fill in the NaN values (outlier B1+ factor values) with
% inpaint_nans function: interpolation of the Laplacian evaluated at non-outlier B1+ factor values 
B1MapNoPks = LiverMask.*inpaint_nans(B1MapNoPks,2);
B1MapNoPks(B1MapNoPks == 0) = NaN;


% Try a second iteration of the extrema to see if it can clean up anything else

[ExtB1Map2,~,~] = IdExtrema(B1MapNoPks, Threshold);
if ~isempty(ExtB1Map2)
    B1MapNoPks2 = B1MapNoPks;
    B1MapNoPks2(ExtB1Map2) = NaN;
    B1MapNoPks2 = LiverMask.*inpaint_nans(B1MapNoPks2,2);
    B1MapNoPks2(B1MapNoPks2 == 0) = NaN;
    B1MapNoExtrema=B1MapNoPks2;
else
   B1MapNoExtrema=B1MapNoPks;
end
if verbose
    % compare Original B1+ Map to new B1+ Map
    screenSize = get(0, 'screenSize');
    figure('position', [screenSize(3:4)*.1, 600, 800]);
    subplot(2,1,1)
    h1 = imagesc(B1Map, [0.5 1.5]);
    set(h1, 'AlphaData', ~isnan(B1Map))
    set(gca, 'Color', [0, 0, 0])
    colorbar
    title ('B1MapFat')
    subplot(2,1,2)
    h2 = imagesc(B1MapNoExtrema, [0.5 1.5]);
    set(h2, 'AlphaData', ~isnan(B1MapNoExtrema))
    set(gca, 'Color', [0, 0, 0])
    colorbar
    title ('B1Map with Extrema Removed')
end

end

