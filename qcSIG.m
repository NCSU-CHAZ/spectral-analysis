function [rVelQC, nPctXcr] = qcSIG(rVel, rCorr, nFsamp, bPlot, sDirName, sFileName)
% --------------------------------qcSIG-----------------------------------%
% Purpose: This function identifies and replace outliers in signature time 
% series via a correlation threshold and SNR (QC criteria). Methods are 
% detailed in Elgar et al. (2005)
% 
% Inputs:
%       - rVel:       3 column vector (X, Y, Z) of evenly sampled time  
%                     series of ADV velocity 
%       - rCorr:      3 column vector (X, Y, Z) of correlation (%) time  
%                     series corresponding to rVel
%       - nFsamp:     used in minimum correlation threshold calculation
%       - bPlot:      boolean for plotting
%       - sDirName:   string for name of new directory containing figures 
%                     ([] for bPlot = false)
%       - sFileName:  string to be appended to figure file name ([] for
%                     bPlot = false)
% Outputs:
%       - rVelQC:     QC'ed velocity series (X, Y, Z)
%
% SEE ALSO: qcADV, qcADV_Harvey
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       11/30/18        KAnarde             Original code
%
%------------------------------preamble-----------------------------------%

disp('-----------------------------------------------------------')
disp('--------------------------qcSIG----------------------------')
disp('-----------------------------------------------------------') 

% calculate minimum correlation threshold for surf zone studies 
% (Elgar et al. (2001 & 2005)) 
nCorrMin = (0.3 + 0.4*sqrt(nFsamp/25)) * 100;    

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~call subroutines~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% 1) identify & replace points below correlation threshold
[rVelQC, nPctXcr] = findReplace(rVel, rCorr);

%-------------------------------------------------------------------------%
%% findReplace
    function [rVelQC, nPctXcr] = findReplace(rVel, rCorr) 
        
        % set smoothed time series to original time series
        [VelXorg, VelXsm] = deal(rVel(:,1));        
        [VelYorg, VelYsm] = deal(rVel(:,2));        
        [VelZorg, VelZsm] = deal(rVel(:,3));
        [VelZ2org, VelZ2sm] = deal(rVel(:,4));
        
        % find the average correlation for all 5 beams
        rCorrMean = mean(rCorr,2);

        % identify logical indices for data < QC criteria 
        lVelXSpike = rCorrMean < nCorrMin;
        lVelYSpike = rCorrMean < nCorrMin;
        lVelZSpike = rCorrMean < nCorrMin;
        lVelZ2Spike = rCorrMean < nCorrMin;

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
        VelZ2sm(lVelZ2Spike) = interp1(iSmpPts(~lVelZ2Spike), ...
                                VelZ2org(~lVelZ2Spike), ...
                                iSmpPts(lVelZ2Spike),'PCHIP');                            

        % save final QC'ed arrays                     
        rVelQC(:,1) = VelXsm;         
        rVelQC(:,2) = VelYsm;         
        rVelQC(:,3) = VelZsm;
        rVelQC(:,4) = VelZ2sm;
        
        % identify spikes associated with low correlations
        iVelSpikeCr = find(rCorrMean < nCorrMin);

        % identify percent of data that are < QC criteria
        nPctXcr    = length(iVelSpikeCr) / length(VelXsm) * 100;
        nPctYcr    = length(iVelSpikeCr) / length(VelYsm) * 100;
        nPctZcr    = length(iVelSpikeCr) / length(VelZsm) * 100;
        %nPctZ2cr   = length(iVelZ2SpikeCr) / length(VelZ2sm) * 100;
        
        if bPlot
 
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
                        plot(iVelSpikeCr, VelXorg(iVelSpikeCr), 'o');
                    xlabel('sample #');
                    ylabel('velocity X (m/s)')
                    plotFancyAxis
                    legend('raw', 'despiked', '< Corr')
                    title(sprintf('%s, %s, %0.2f%% < Corr', sDirName, ...
                          sFileName, nPctXcr))
                    hold off

            % velocity Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            subplot(3,1,2)
            hAx(2) = gca;   
                        plot(VelYorg);   
                    hold on
                        plot(VelYsm); 
                    hold on
                        plot(iVelSpikeCr, VelYorg(iVelSpikeCr), 'o');
                    xlabel('sample #');
                    ylabel('velocity Y (m/s)')
                    plotFancyAxis
                    legend('raw', 'despiked', '< Corr')
                    title(sprintf('%s, %s, %0.2f%% < Corr', sDirName, ...
                          sFileName, nPctYcr))
                    hold off

            % velocity Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            subplot(3,1,3)
            hAx(3) = gca; 
                        plot(VelZorg);   
                    hold on
                        plot(VelZsm); 
                    hold on
                        plot(iVelSpikeCr, VelZorg(iVelSpikeCr), 'o');
                    xlabel('sample #');
                    ylabel('velocity Z (m/s)')
                    plotFancyAxis
                    legend('raw', 'despiked', '< Corr')
                    title(sprintf('%s, %s, %0.2f%% < Corr', sDirName, ...
                          sFileName, nPctZcr))
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