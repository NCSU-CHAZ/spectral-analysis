function [rVelQC] = qcADV_Harvey(rVel, rCorr, rSNR, nFsamp, nSNR, ...
                                  bPlot, sDirName, sFileName)
% --------------------------------qcADV-----------------------------------%
% Purpose: This function identifies and replace outliers in ADV time series 
% via a correlation threshold and SNR (QC criteria). Methods are detailed 
% in Elgar et al. (2005)
% 
% Inputs:
%       - rVel:       3 column vector (X, Y, Z) of evenly sampled time  
%                     series of ADV velocity 
%       - rCorr:      3 column vector (X, Y, Z) of correlation (%) time  
%                     series corresponding to rVel
%       - rSNR:       3 column vector (X, Y, Z) of SNR (dB) time  
%                     series corresponding to rVel
%       - nFsamp:     used in minimum correlation threshold calculation
%       - nSNR:       3 column vector (X, Y, Z) of SNR threshold from 
%                     instrument manufacturer; for Nortek ADVs, we use 
%                     nSNR = 8-15 (e.g. [8 8 8])
%       - bPlot:      boolean for plotting
%       - sDirName:   string for name of new directory containing figures 
%                     ([] for bPlot = false)
%       - sFileName:  string to be appended to figure file name ([] for
%                     bPlot = false)
% Outputs:
%       - rVelQC:     QC'ed velocity series (X, Y, Z)
%
% SEE ALSO: none
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       05/11/17        KAnarde             Original code
%       04/27/18        " "                 Modified to window both the 
%                                           input velocity and pressure 
%                                           timeseries for Harvey
%
%------------------------------preamble-----------------------------------%

disp('-----------------------------------------------------------')
disp('--------------------------qcADV----------------------------')
disp('-----------------------------------------------------------') 

% calculate minimum correlation threshold for surf zone studies 
% (Elgar et al. (2001 & 2005)) 
nCorrMin = (0.3 + 0.4*sqrt(nFsamp/25)) * 100;    

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~call subroutines~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% 1) identify & replace points below correlation and SNR thresholds
[rVelQC] = findReplace(rVel, rCorr, rSNR, nSNR);

%-------------------------------------------------------------------------%
%% findReplace
    function [rVelQC] = findReplace(rVel, rCorr, rSNR, nSNR) 
        
        % set smoothed time series to original time series
        [VelXorg, VelXsm] = deal(rVel(:,1));        
        [VelYorg, VelYsm] = deal(rVel(:,2));        
        [VelZorg, VelZsm] = deal(rVel(:,3));

        % identify logical indices for data < QC criteria 
        lVelXSpike = rCorr(:,1) < nCorrMin | rSNR(:,1)  < nSNR(1);
        lVelYSpike = rCorr(:,2) < nCorrMin | rSNR(:,2)  < nSNR(2);
        lVelZSpike = rCorr(:,3) < nCorrMin | rSNR(:,3)  < nSNR(3);

        % spline over the outlier indices
        iSmpPts = (1:length(VelXsm))';   % sample points for interp1
        VelXsm(lVelXSpike) = interp1(iSmpPts(~lVelXSpike), ...
                                        VelXorg(~lVelXSpike), ...
                                        iSmpPts(lVelXSpike),'PCHIP');                 
        VelYsm(lVelYSpike) = interp1(iSmpPts(~lVelYSpike), ...
                                        VelYorg(~lVelYSpike), ...
                                        iSmpPts(lVelYSpike),'PCHIP');         
        VelZsm(lVelZSpike) = interp1(iSmpPts(~lVelZSpike), ...
                                        VelZorg(~lVelZSpike), ...
                                        iSmpPts(lVelZSpike),'PCHIP');

        % save final QC'ed arrays                     
        rVelQC(:,1) = VelXsm;         
        rVelQC(:,2) = VelYsm;         
        rVelQC(:,3) = VelZsm;
        
        if bPlot

            % identify spikes associated with low correlations
            iVelXSpikeCr = find(rCorr(:,1) < nCorrMin);
            iVelYSpikeCr = find(rCorr(:,2) < nCorrMin);
            iVelZSpikeCr = find(rCorr(:,3) < nCorrMin);

            % identify spikes associated with low SNR
            iVelXSpikeSNR = find(rSNR(:,1)  < nSNR(1));
            iVelYSpikeSNR = find(rSNR(:,2)  < nSNR(2));
            iVelZSpikeSNR = find(rSNR(:,3)  < nSNR(3));

            % identify percent of data that are < QC criteria
            nPctXcr    = length(iVelXSpikeCr) / length(VelXsm) * 100;
            nPctYcr    = length(iVelYSpikeCr) / length(VelYsm) * 100;
            nPctZcr    = length(iVelZSpikeCr) / length(VelZsm) * 100;
            nPctXsnr   = length(iVelXSpikeSNR) / length(VelXsm) * 100;
            nPctYsnr   = length(iVelYSpikeSNR) / length(VelYsm) * 100;
            nPctZsnr   = length(iVelZSpikeSNR) / length(VelZsm) * 100;
        
            % check if figure directory has been created 
            if exist(sDirName, 'dir') == 0      % 7 specifies a folder
                mkdir(sprintf('%s', sDirName))  % save figs to new folder
            end

            % plot original vs. qc'ed time series
            % velocity X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            f(1) = figure; 
            subplot(3,1,1)
            hAx(1) = gca; 
                        plot(VelXorg);   
                    hold on
                        plot(VelXsm); 
                    hold on
                        plot(iVelXSpikeCr, VelXorg(iVelXSpikeCr), 'o');
                    hold on
                        plot(iVelXSpikeSNR, VelXorg(iVelXSpikeSNR), 'x');
                        xlabel('sample #');
                        ylabel('velocity X (m/s)')
                        plotFancyAxis
                        if ~isempty(iVelXSpikeSNR)
                            legend('raw', 'despiked', '< Corr', '< SNR') 
                        else
                            legend('raw', 'despiked', '< Corr')
                        end
                        title(sprintf(...
                                 '%s, %s, %0.2f%% < Corr, %0.2f%% < SNR', ...
                                 sDirName, sFileName, nPctXcr, nPctXsnr))
                    hold off

            % velocity Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            subplot(3,1,2)
            hAx(2) = gca;   
                        plot(VelYorg);   
                    hold on
                        plot(VelYsm); 
                    hold on
                        plot(iVelYSpikeCr, VelYorg(iVelYSpikeCr), 'o');
                    hold on
                        plot(iVelYSpikeSNR, VelYorg(iVelYSpikeSNR), 'x');
                        xlabel('sample #');
                        ylabel('velocity Y (m/s)')
                        plotFancyAxis
                        if ~isempty(iVelYSpikeSNR)
                            legend('raw', 'despiked', '< Corr', '< SNR') 
                        else
                            legend('raw', 'despiked', '< Corr')
                        end
                        title(sprintf(...
                                 '%s, %s, %0.2f%% < Corr, %0.2f%% < SNR', ...
                                 sDirName, sFileName, nPctYcr, nPctYsnr))
                    hold off

            % velocity Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            subplot(3,1,3)
            hAx(3) = gca; 
                        plot(VelZorg);   
                    hold on
                        plot(VelZsm); 
                    hold on
                        plot(iVelZSpikeCr, VelZorg(iVelZSpikeCr), 'o');
                    hold on
                        plot(iVelZSpikeSNR, VelZorg(iVelZSpikeSNR), 'x');
                        xlabel('sample #');
                        ylabel('velocity Z (m/s)')
                        plotFancyAxis
                        if ~isempty(iVelZSpikeSNR)
                            legend('raw', 'despiked', '< Corr', '< SNR') 
                        else
                            legend('raw', 'despiked', '< Corr')
                        end
                        title(sprintf(...
                                 '%s, %s, %0.2f%% < Corr, %0.2f%% < SNR', ...
                                 sDirName, sFileName, nPctZcr, nPctZsnr))
                    hold off

            % link axis for zooming
            linkaxes(hAx, 'x')

            % save fig
            saveas(f, [pwd sprintf('/%s/%s-ts.fig', sDirName, sFileName)]);
            close(f)
        end
    end
%-------------------------------------------------------------------------%
end