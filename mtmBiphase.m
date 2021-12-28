function [rBphase] = mtmBiphase(zBspec, rBic1, rF, rFbnd, nNW, nM, bPlot) 
% ------------------------------mtmBiphase--------------------------------%
% Purpose: This function calculates the biphase over a specified integral
% from bispectral estimates generated using the multitaper technique
%
% Inputs:
%       - zBspec:    complex bispectrum
%       - rBic1:     bicoherence (Haubrich formulation, upper bound of 1) 
%       - rF:        frequency matrix [Hz,Hz]
%       - rFbnd:     frequency bounds [f11 f12 f21 f22]
%       - nNW:       time bandwidth product (for Harvey, NW=3)
%       - nM:        number of segments averaged for zBspec 
%       - bPlot:     boolean for plotting the bphase up to 0.25 Hz
%
% SEE ALSO: mtmBispecd.m, plotMTMBispecdTACC.m 
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       08/20/18        KAnarde             Original code
%
%% MAKE FILTERS -----------------------------------------------------------

% make symmetric filter for bispectrum matrix 
nCol   = size(rBic1,1);
rSym   = NaN(nCol,nCol); 
cntRow = 1;
for i = 1 : nCol
    rSym(1:cntRow,i) = 1;
    cntRow = cntRow + 1;
end

% calculate 95% confidence limit for zero bicoherence 
nCL = sqrt( 1 / (nNW^2*nM) );

% create filter for significant bicoherence 
rCL = zeros(size(rBic1));
rCL(rBic1 >= nCL) = 1; 

% final filter for confidence limits and symmetry
rFilt  = rSym .* rCL;

%% CALCULATE BIPHASE ------------------------------------------------------

% apply filters
zBfilt = zBspec .* rFilt;

% plot the biphase
if bPlot
   tmpBphase = rad2deg(atan2(imag(zBfilt),real(zBfilt)));
   figure
   imAlpha=ones(size(tmpBphase));
   imAlpha(isnan(tmpBphase))=0;
   imagesc(rF, rF, tmpBphase,'AlphaData',imAlpha); grid on
   set(gca,'color','white');
   title('$\beta (^o)$', 'interpreter', 'Latex', 'FontSize', 14)
   xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
   ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
   colormap([cmap('slateblue',50,20,0);flipud(cmap('orange',50,20,0))]); 
   colorbar
   caxis([-180 180])
   xlim([0 0.25]); ylim([0 0.25]);
   axis('square')
   set(gca, 'YDir', 'normal')
end

% isolate frequency band of interest
[~, iF11] = min(abs(rF - rFbnd(1)));
[~, iF12] = min(abs(rF - rFbnd(2)));
[~, iF21] = min(abs(rF - rFbnd(3)));
[~, iF22] = min(abs(rF - rFbnd(4)));
zB = zBfilt(iF21:iF22, iF11:iF12);

% sum over rows and columns
zBsum = sum(sum(zB,'omitnan'));

% calculate the biphase (converts a complex number to polar coordinates)
%rBphase = atand(imag(zBsum)/real(zBsum)); % this only gives [-90 90]
rBphase = atan2(imag(zBsum),real(zBsum)); % [-180 180], angle() is the same
rBphase = rad2deg(rBphase);

end