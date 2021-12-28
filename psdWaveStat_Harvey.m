function [waveStat] = szWaveStat_Harvey(rP, dtP, nFminIF, nFminSS, nFmaxSS, n2n, ...
                       bPlot, sDirName, nHpt, nWin, sFileName, nZbeg, nZend)
% --------------------------szWaveStat_Harvey-----------------------------%
% Purpose: This function generates wave statistics from bottom mounted 
% pressure transducers deployed durring Hurricane Harvey using spectral 
% analysis. The method employed herein is based on Dean and 
% Dalyrymple [1991] & Jones and Monismith (2007). Runs psd.f written and 
% provided by Robert L. Parker. Modified from szWaveStat so that it is not
% adaptive for finding Fmax, incorporates prolate spheroidal tapers, and
% also calculates statistics using LPF and BPF in the time domain.
% 
% Inputs:
%       - rP:         evenly sampled time series (vector) of gauge pressure 
%                     (dbar, corrected for local variations in 
%                     atmospheric pressure); time series should be windowed 
%                     evenly (e.g. 12:00 on 1/24 - 12:00 on 1/27) if there 
%                     will be multiple bursts (see nWin below for single
%                     burst < 1 hr)
%       - dtP:        time series of datetime or datenum values ass with rP
%       - nFminIF:    min frequency used for windowing infragravity band
%       - nFminSS:    min frequency used for windowing the sea-swell band 
%       - nFmaxSS:    max frequency that the PTs can detect 
%                     (recommended value = 0.7 Hz)
%       - n2n:        to maintain stationarity and reduce computation time, 
%                     we compute a psd for a single burst of 2^n samples at 
%                     the start of each (hour) window: 
%                           (e.g. at 16 Hz, 3600×16 = 57,600 samp/hr, so 
%                           max n to maintain stationarity, while long 
%                           enough to represent one hr of data, is 14: 
%                           2^14 = 16,384 samples, 16,384/16 = 1024 s 
%                           ~ 17 min; for 2 Hz --> 3600*2=7200 samp/hr,  
%                           n=11, 2^11=2048 samp, 2048/2=1024 s)
%       - bPlot:     boolean for plotting
%       - sDirName:  string for name of new directory containing figures 
%                    ([] for bPlot = false)
%       - nHpt:      height of pressure transducer off (above) the bed (m)
%                    (this is typically zero, do not include negative
%                    numbers...burial is accounted for using poroelastic
%                    theory)
%       - nFinf:     minimum infragravity frequency to separate from tides
%                    (rec. value = 0.005 Hz)
%       - nWin:      number of samples/window (e.g. for a 1-hr window
%                    (default) nWin = 3600 sec * 16 Hz); if the time 
%                    series is < 1 hr and there is only one burst, nWin is  
%                    simply the length of the input time series
%       - sFileName: string to be appended to figure file name ([] for
%                    bPlot = false)
%       - nZbeg:     depth below the bed surface at beginning of time 
%                    series(positive downwards) - will linearly interpolate
%                    to end of time series using
%       - nZend:     depth below the bed surface at end of time series 
%
% Outputs:
%   WaveStat (structure)
%       - Hrms:     root-mean-square wave heigths for each burst [m]
%       - Hm0:      zeroth moment wave heigths for each burst [m]
%       - Fpeak:    peak frequency defined by max(Snn) [Hz]
%       - Tm01:     mean wave period [sec]
%       - Tm02:     mean wave period def by the zero-upcrossing [sec]
%       - SnnWin:   water surface elevation variance spectra [m^2/Hz]
%       - SnnNorm:  normalized Snn [1/Hz]
%       - Pburst:   array of input data organized by burst [Pa]
%       - h:        mean water depth [m]
%       - f:        frequency array [Hz]
%
% SEE ALSO: runPSD.m, psd.f, szWaveStat
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       02/20/17        KAnarde             Original code
%       10/21/17        " "                 Updates for Harvey
%
%------------------------------preamble-----------------------------------%

disp('-----------------------------------------------------------')
disp('------------------------szWaveStat-------------------------')
disp('-----------------------------------------------------------') 

% nested variables for subroutines
sPath = '/Users/KatherineAnardeWheels/Research/Software/BDW_Scripts/psd';

% constants and conversion parameters
nRho     = 1025;   % seawater density [kg/m3]
nG       = 9.81;   % gravitational constant [m/s2]
nDbar2Pa = 10000; 

% boolean for plotting QC steps
bPlotQC  = false;

% define dimensions for preallocating arrays
nFsamp = 1 / mean(seconds(diff(dtP)));  % [Hz] sampling frequency
if ~exist('nWin','var') || isempty(nWin)% # of samples per window 
    nWin = nFsamp * 3600;               % default = 1-hour windows
end
nBurst = 2^n2n;                         % # of samples per burst 
nTotal = length(rP);                    % total # of indices 
nPSD   = floor(nTotal/nWin);            % number of PSDs (bursts)

% throw error for bad choice of n2n
if nBurst > nTotal
    error('2^n is too large for pressure time series')
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~call subroutines~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% 1) window pressure data to burst intervals 
[rPb, dtPb] = winPress(rP, dtP);

% 2) calculate dynamic pressure, rPd [Pa], and the mean water depth, nH [m]
[rPd, nH] = dynPress(rPb, nHpt);

% 3a) calculate Snn
[rSnnPS, rSnnMS, ~, ~, rF] = genSnn(rPb, bPlotQC, dtPb);

% 3b) look for change points (stationarity) using Rahim & Thomson method
%[stCP] = changePt(rPb, dtPb);

% 4) account for attenuation of the pressure signal by burial and water
%    depth
[rSnnPSkp, rSnnMSkp, rKp, rEkz] = depthAtt(rSnnPS, rSnnMS, ...
                                               rF, nH, nHpt, nZbeg, nZend);

% 6) window to min and max frequency
[rSnnPSss, rSnnPSinf, rSnnMSinf, rSnnMSss, ~] = ...
                               winSnn(rSnnPS, rSnnMS, rKp, rEkz, rF, nFminSS, ...
                                      nFminIF, nFmaxSS, nHpt, nH);
                                  
% 7) time domain - low and bandpass filters  
[rFSElpf, rFSEbpf] = timeDomain(rPd, dtPb, nH);
    
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~wave statistics~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%%% prolate spheroidal tapers
%%% infragravity band ----------------------------------------------------
% spectral moments (X(F) must be a vector)
psM0inf   = trapz(rF(:,1), rSnnPSinf);     % [m^2] = variance (sigma^2)
psM1inf   = trapz(rF(:,1), rSnnPSinf .* rF);
psM2inf   = trapz(rF(:,1), rSnnPSinf .* rF.^2);

psNrmsInf = sqrt(psM0inf);                 % std deviation (sigma)
waveStat(1).psHrmsInf = 2 * psNrmsInf;     % rms wave height
waveStat(1).psHm0Inf  = 4.004 * psNrmsInf; % 0-moment wave height (~H1/3)
waveStat(1).psTm01Inf = psM0inf ./ psM1inf;       % mean wave period 
waveStat(1).psTm02Inf = sqrt(psM0inf ./ psM2inf); % mean wave period: 0-up

% peak frequency (in each column)
[~, iMaxPSinf] = max(rSnnPSinf, [], 1);
waveStat(1).psFpeakInf = rF(iMaxPSinf);

% deep water wavelength
waveStat(1).psL0inf = (nG * waveStat(1).psTm01Inf.^2) / (2*pi); 

% normalize spectra for comparison of spectral shapes (area=1)
waveStat(1).psSnnNormInf = rSnnPSinf ./ repmat(psM0inf, nBurst/2+1, 1);

%%% sea-swell band -------------------------------------------------------
% spectral moments 
psM0ss   = trapz(rF(:,1), rSnnPSss);           % [m^2] = variance (sigma^2)
psM1ss   = trapz(rF(:,1), rSnnPSss .* rF);
psM2ss   = trapz(rF(:,1), rSnnPSss .* rF.^2);

psNrmsSS = sqrt(psM0ss);                       % std deviation (sigma)
waveStat(1).psHrmsSS = 2 * psNrmsSS;           % rms wave height
waveStat(1).psHm0SS  = 4.004 * psNrmsSS;       % 0-moment wave height (~H1/3)
waveStat(1).psTm01SS = psM0ss ./ psM1ss;       % mean wave period 
waveStat(1).psTm02SS = sqrt(psM0ss ./ psM2ss); % mean wave period by zero-upcrs

% peak frequency (in each column)
[~, iMaxPSss] = max(rSnnPSss, [], 1);
waveStat(1).psFpeakSS = rF(iMaxPSss);

% deep water wavelength
waveStat(1).psL0ss = (nG * waveStat(1).psTm01SS.^2) / (2*pi); 

% normalize spectra for comparison of spectral shapes (area=1)
waveStat(1).psSnnNormSS = rSnnPSss ./ repmat(psM0ss, nBurst/2+1, 1);

%%% multisine tapers
%%% infragravity band ----------------------------------------------------
% spectral moments (X(F) must be a vector)
msM0inf   = trapz(rF(:,1), rSnnMSinf);     % [m^2] = variance (sigma^2)
msM1inf   = trapz(rF(:,1), rSnnMSinf .* rF);
msM2inf   = trapz(rF(:,1), rSnnMSinf .* rF.^2);

msNrmsInf = sqrt(msM0inf);                 % rms value of the surface elev
waveStat(1).msHrmsInf = 2 * msNrmsInf;     % rms wave height
waveStat(1).msHm0Inf  = 4.004 * msNrmsInf; % 0-moment wave height (~H1/3)
waveStat(1).msTm01Inf = msM0inf ./ msM1inf;       % mean wave period 
waveStat(1).msTm02Inf = sqrt(msM0inf ./ msM2inf); % mean wave period: 0-up

% peak frequency (in each column)
[~, iMaxMSinf] = max(rSnnMSinf, [], 1);
waveStat(1).msFpeakInf = rF(iMaxMSinf);

% deep water wavelength
waveStat(1).msL0inf = (nG * waveStat(1).msTm01Inf.^2) / (2*pi); 

% normalize spectra for comparison of spectral shapes (area=1)
waveStat(1).msSnnNormInf = rSnnMSinf ./ repmat(msM0inf, nBurst/2+1, 1);

%%% sea-swell band -------------------------------------------------------
% spectral moments 
msM0ss   = trapz(rF(:,1), rSnnMSss);          % [m^2] = variance^2
msM1ss   = trapz(rF(:,1), rSnnMSss .* rF);
msM2ss   = trapz(rF(:,1), rSnnMSss .* rF.^2);

msNrmsSS = sqrt(msM0ss);                      % rms value of the surface elev
waveStat(1).msHrmsSS = 2 * msNrmsSS;          % rms wave height
waveStat(1).msHm0SS  = 4.004 * msNrmsSS;      % 0-moment wave height (~H1/3)
waveStat(1).msTm01SS = msM0ss ./ msM1ss;       % mean wave period 
waveStat(1).msTm02SS = sqrt(msM0ss ./ msM2ss); % mean wave period by zero-upcrs

% peak frequency (in each column)
[~, iMaxMSss] = max(rSnnMSss, [], 1);
waveStat(1).msFpeakSS = rF(iMaxMSss);

% deep water wavelength
waveStat(1).msL0ss = (nG * waveStat(1).msTm01SS.^2) / (2*pi);

% normalize spectra for comparison of spectral shapes (area=1)
waveStat(1).msSnnNormSS = rSnnMSss ./ repmat(msM0ss, nBurst/2+1, 1);

%%% save other variables to structure ------------------------------------
waveStat(1).Pb      = rPb;
waveStat(1).psSnnInf= rSnnPSinf;
waveStat(1).psSnnSS = rSnnPSss;
waveStat(1).msSnnInf= rSnnMSinf;
waveStat(1).msSnnSS = rSnnMSss;
waveStat(1).dtPb    = dtPb;
waveStat(1).FSElpf  = rFSElpf;
waveStat(1).h       = nH;
waveStat(1).f       = rF;

%%% time domain statistics 
%%% infragravity band ----------------------------------------------------
waveStat(1).HsLPF = 4 * std(rFSElpf);
waveStat(1).HsBPF = 4 * std(rFSEbpf);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
if bPlot
    
    % preallocate figure array
    f = NaN(size(nPSD));
    
    % find max of Y value for plots
    nYmax = max(cat(1, rSnnPSinf(:), rSnnPSss(:), rSnnMSss(:), rSnnMSinf(:)));
    
    % check if figure directory has been created 
    if exist(sDirName, 'dir') == 0      % 7 specifies a folder
        mkdir(sprintf('%s', sDirName))  % save figs to new folder
    end
    
    nScrsz = get(groot, 'ScreenSize');
    
    for n = 1 : nPSD
        % plot Snn spectra (prolate spheroidal)
        f(n) = figure;
        subplot(1,2,1)
              loglog(rF(:,n), rSnnPSss(:,n), '-.', ...
                     rF(:,n), rSnnPSinf(:,n), '-.', ...
                     rF(:,n), rSnnPSkp(:,n), ...
                     rF(:,n), rSnnPS(:,n), '--');
        legend('PS Snn_{SS}', 'PS Snn_{inf}', 'PS Snn_{Kp}', ' PS Snn')
        xlabel('frequency (Hz)'); 
        ylabel('water surface elevation variance [m^2/Hz]')
        set(gca, 'Xlim', [0 8] , 'Ylim', [0 nYmax]) 
        plotFancyAxis
        title(sprintf('%s, h = %f m, Hm0_{Inf} = %f m, Hm0_{SS} = %f m',...
                      dtPb(1,n), nH(1,n), ...
                      waveStat(1).psHm0Inf(1,n), waveStat(1).psHm0SS(1,n)))
                  
        subplot(1,2,2)
              loglog(rF(:,n), rSnnMSss(:,n), '-.', ...
                     rF(:,n), rSnnMSinf(:,n), '-.', ...
                     rF(:,n), rSnnMSkp(:,n), ...
                     rF(:,n), rSnnMS(:,n), '--');
        legend('MS Snn_{SS}', 'MS Snn_{inf}', 'MS Snn_{Kp}', ' MS Snn')
        xlabel('frequency (Hz)'); 
        ylabel('water surface elevation variance [m^2/Hz]')
        set(gca, 'Xlim', [0 8] , 'Ylim', [0 nYmax]) 
        plotFancyAxis
        title(sprintf('Hm0_{Inf} = %f m, Hm0_{SS} = %f m',...
                      waveStat(1).msHm0Inf(1,n), waveStat(1).msHm0SS(1,n))) 
        
        % make figure large and save          
        set(f(n), 'Position', [1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);  
        saveas(f(n), [pwd sprintf('/%s/%s-burst%d.fig', sDirName, ...
                      sFileName, n)]);
        close(f(n))
    end
end
        
% ----------------------------------------------------------------------- %
%% winPress
    function [rPb, dtPb] = winPress(rP, dtP)
    %%% window pressure time series
    
    [rPwin] = deal(NaN(nWin,nPSD));    % preallocate array
    cntWin  = 0;                       % initialize counter for windowing
    
    % convert rP from dbar to Pa
    rP      = rP * nDbar2Pa;
    
        for i = 1 : nPSD
            % window pressure and datetime arrays
            rPwin(:,i) = rP( (1 + cntWin) : (cntWin + nWin) );
            dtPwin(:,i)= dtP((1 + cntWin) : (cntWin + nWin) );
    
            % increase window counter
            cntWin = cntWin + nWin;
        end
        
    % further window pressure and datetime arrays to the # of burst samples
    rPb  = rPwin(1:nBurst,:);
    dtPb = dtPwin(1:nBurst,:);   
    
    end
% ----------------------------------------------------------------------- %
%% dynPress
    function [rPd, nH] = dynPress(rPb, nHpt)
    %%% dynamic pressure is found by subtracting the mean (hydrostatic)
    %%% pressure from the sea (gauge) pressure (absolute - atmospheric)
    %%% hydrostatic approximation: n(-h)= (P(-h)-rho*g*h)/(rho*g*Kp(-h)) 
    
    % calculate dynamic pressure along each column for the burst samples
    rPd = rPb - repmat(mean(rPb,1), nBurst, 1);

    % calculate mean water level (m), accounting for nHpt (height of PT
    % above the bottom); NOTE: the mean of each burst is also given in the  
    % log file from psd.f using linear least squares fit
    nH = ( mean(rPb,1) / (nRho * nG) ) + nHpt; 
    
    end
% ----------------------------------------------------------------------- %
%% genSnn    
    function [rSnnPS, rSnnMS, rSnnPSwht, rSnnMSwht, rF] = genSnn(rPb, bPlotQC, dtPb)
    %%% calculate Spp and (Snn) using both multisine and prolate spheroidal 
    %%% tapering methods
    
    % preallocate arrays
    [rF, rSnnPS, rSnnMS, rSnnPSwht, rSnnMSwht, rSppPS, rSppMS, ...
     rSppPSwht, rSppMSwht] = deal(NaN(nBurst/2 + 1, nPSD));
    
    for i = 1 : nPSD 
        
        % run psd.f
        % (NOTE: because psd.f removes the mean prior to calc. of the PSD, 
        % this is the dynamic pressure power spectra, Spp)
        % multisine tapers (adaptive), no prewhitening
        bPwht = false;
        sFile = 'MS';
        nPrltTBP = []; % time-bandwidth product for PS tapers
        [rSppMS(:,i), rF(:,i)] = runPsd(rPb(:,i), nFsamp, sFile, sPath,...
                                          nPrltTBP, bPwht, false, []);
        % calculate Snn (no depth attenuation)
        rSnnMS(:,i) = rSppMS(:,i) / (nRho^2 * nG^2);                              
          
        % prolate spheroidal tapers (non-adaptive), no prewhitening
        bPwht = false;
        sFile = 'PS';
        nPrltTBP = 4; % time-bandwidth product for PS tapers
        [rSppPS(:,i)] = runPsd(rPb(:,i), nFsamp, sFile, sPath,...
                                          nPrltTBP, bPwht, false, []);
        % calculate Snn (no depth attenuation)
        rSnnPS(:,i) = rSppPS(:,i) / (nRho^2 * nG^2); 

        % multisine tapers (adaptive), prewhitening
        bPwht = true;
        sFile = 'MSwht';
        nPrltTBP = []; % time-bandwidth product for PS tapers
        [rSppMSwht(:,i)] = runPsd(rPb(:,i), nFsamp, sFile, sPath,...
                                          nPrltTBP, bPwht, false, []);
        % calculate Snn (no depth attenuation)
        rSnnMSwht(:,i) = rSppMSwht(:,i) / (nRho^2 * nG^2);
        
        % prolate spheroidal tapers (non-adaptive), prewhitening
        bPwht = true;
        sFile = 'PSwht';
        nPrltTBP = 4; % time-bandwidth product for PS tapers
        [rSppPSwht(:,i)] = runPsd(rPb(:,i), nFsamp, sFile, sPath,...
                                          nPrltTBP, bPwht, false, []);
        % calculate Snn (no depth attenuation)
        rSnnPSwht(:,i) = rSppPSwht(:,i) / (nRho^2 * nG^2);
        
        % pwelch (default settings using 8 segments (for n=15 at RRU-1, 
        % each window is about 4.2 minutes), Hamming windowed, 50%
        % overlap)
        [rSppPW(:,i), rFpw(:,i)] = pwelch(detrend(rPb(:,i)), [], [], [], nFsamp);
        rSnnPW(:,i) = rSppPW(:,i) / (nRho^2 * nG^2);
        
        if bPlotQC
            % added for ADV sensitivity testing
            if (60 < i) && (i < 80)
                
%             % plot windowed water level
%             rNb = rPb(:,i) / (nRho * nG);
%             figure; plot(rNb);   
%                     plotFancyAxis
%                     xlabel('sample #')
%                     ylabel('water level (m)')
%                     title(sprintf('Burst %d', i))
            
            % plot spectra
            figure; loglog(rF(:,i), rSnnMS(:,i), ...
                           rF(:,i), rSnnMSwht(:,i), ...
                           rF(:,i), rSnnPS(:,i), ...
                           rF(:,i), rSnnPSwht(:,i), ...
                           rFpw(:,i), rSnnPW(:,i)); 
                    legend('multisine tapers', ...
                           'multisine tapers - prewhitened', ...
                           'prolate spheroidal tapers', ...
                           'prolate spheroidal tapers - prewhitened', ...
                           'welch - 50% overlap')
                    xlabel('frequency [Hz]'); ylabel('S_{\eta\eta} [m^2/Hz]')
                    plotFancyAxis
                    title(sprintf('%s, Burst %d', dtPb(1,i), i))
            end
        end                                                                                   
    end
    end
% ----------------------------------------------------------------------- %
%% changePt
    function [stCP] = changePt(rPb, dtPb)
    %%% implement Rahim and Thomson method for identifying change points
    %%% from overlapping multi-taper spectrogram to evaluate stationarity 
    %%% of windows. 
    
    % spectrogram program parameters ------------------------------------ %
    % full window size
    nSampSz = nBurst; % e.g. 32768 samples (2^15)
    % block size
    %nBlckSz = [2^(n2n - 4), 2^(n2n - 3), 2^(n2n - 2), 2^(n2n - 1)]; %[2048, 4096, 8192, 16384]
    nBlckSz = [2^(n2n - 1), 2^(n2n - 1), 2^(n2n - 1), ...
        2^(n2n - 1), 2^(n2n - 1), 2^(n2n - 1)]; % for testing time-bandwidth
    % offset for blocks with 75% overlap
    nOffset = nBlckSz * 0.25; %0.25; 
    nBlocks = (nSampSz - nBlckSz) ./ nOffset + 1;
    nTB = 1 : 6;
    
    % preallocate structure
    stCP = struct('Snn', {}, 'SnnInf', {}, 'f', {}, 'win', {}, 'dt', {});
    %[rSnnPSinf, rSnnPSss, rSnnMSinf, rSnnMSss] = deal(zeros(nBurst/2 + 1, nPSD));
    
    % create spectrogram arrays for each block size --------------------- %
    for i2n = 2 : 6 %: length(nBlckSz)
        
         % preallocate array
        [rCPwin, dtCPwin] = deal(NaN(nBlckSz(i2n), nBlocks(i2n), 1)); 
        %[rCPwin] = NaN(nBlckSz(i2n), nBlocks(i2n), size(rPb,2));
    
        for iBurst = 1 : size(rPb,2)
        
            cntOff = 0; % initialize counter for offset
            
            for iBlck = 1 : nBlocks(i2n)
            
                % create overlapping blocks
                rCPwin(:, iBlck, iBurst) = rPb( (1 + cntOff) : (nBlckSz(i2n) + cntOff), iBurst);
                dtCPwin(:, iBlck, iBurst) = datenum(dtPb((1 + cntOff) : (nBlckSz(i2n) + cntOff), iBurst));

                % increase offset counter
                cntOff = cntOff + nOffset(i2n);
                
                % compute the PSD for each block using prolate spheroidal 
                % tapers with a time-bandwidth parameter = 4 (# of freq.
                % bins over which the spectrum is averaged) 
                bPwht = false;
                sFile = 'CP';
                nPrltTBP = nTB(i2n); % time-bandwidth product for PS tapers
                [rSppCP, rFCP] = runPsd(rCPwin(:, iBlck, iBurst), ...
                                    nFsamp, sFile, sPath, nPrltTBP, ...
                                    bPwht, false, []);
                                
                % find the windowing indices for the infragravity spectra
                [~, iFminIF] = min(abs(rFCP - nFminIF));
                [~, iFminSS] = min(abs(rFCP - nFminSS));
                
                % save to structure
                stCP{i2n}.Snn(:, iBlck, iBurst) = rSppCP / (nRho^2 * nG^2);
                stCP{i2n}.f(:, iBlck, iBurst) = rFCP;
                stCP{i2n}.SnnInf(:, iBlck, iBurst) = ...
                             stCP{i2n}.Snn(iFminIF:iFminSS, iBlck, iBurst);
                stCP{i2n}.fInf(:, iBlck, iBurst) = ...
                             stCP{i2n}.f(iFminIF:iFminSS, iBlck, iBurst);
                
                % find the dynamic pressure for plotting                   
                rPdwin = (rCPwin(:, iBlck, iBurst) - ...
                         repmat(mean(rCPwin(:, iBlck, iBurst), 1), ...
                         nBlckSz(i2n), 1))/(nRho * nG);
                                     
                     if bPlot
                        % plot each of the burst time series and the spectra
                        figure;
                        subplot(1,2,1)
                        plot(dtCPwin(:, iBlck, iBurst), rPdwin)
                        ylabel('free surface elevation (m)')
                        datetick('x','mm/dd HH PM')
                        plotFancyAxis

                        subplot(1,2,2)
                        loglog(stCP{i2n}.f(:, iBlck, iBurst), ...
                               stCP{i2n}.Snn(:, iBlck, iBurst))
                        ylabel('wave power spectra, S_{\eta\eta}(f) [m^2/Hz]')
                        set(gca, 'XLim', [0 0.05])
                        plotFancyAxis
                     end
            end
         
            if bPlot
                % plot the spectrogram for each burst
                figure;
                subplot(1,2,1)
                surf(repmat(dtCPwin(1, :, iBurst), size(stCP{i2n}.f, 1), 1), ...
                     stCP{i2n}.f(:, :, iBurst), stCP{i2n}.Snn(:, :, iBurst), ...
                     'EdgeColor', 'none')
                axis xy; view(2)
                datetick('x','mm/dd HH PM')
                axis tight
                % formatting
                c = colorbar;
                ylabel(c, 'S_{\eta\eta}(f) [m^2/Hz]', ...
                       'FontSize', 13);
                ylabel('frequency (Hz)');
                title(sprintf('change point spectrogram, block size = %d', nBlckSz(i2n)));
                % make it fancy
                plotFancyAxis()
                set(gca, 'YGrid', 'off', 'YLim', [0 0.05]) %, 'CLim', [min(rSn30m(:)) max(rSn30m(:))]);

                % plot all the spectra
                subplot(1,2,2)
                loglog(stCP{i2n}.f(:, :, iBurst), stCP{i2n}.Snn(:, :, iBurst))
                xlabel('frequency')
                ylabel('S_{\eta\eta}(f) [m^2/Hz]')
                plotFancyAxis()
                set(gca, 'XLim', [0 0.05])
            end
        end
        
        % save spectrogram variables to structure for each block size
        stCP{i2n}.win = rCPwin;
        stCP{i2n}.dt  = dtCPwin;
        
        %%% change point statistics 
        %%% NOTE: these are a rough estimate and will be slightly off
        %%% because they are not corrected for the various transfer
        %%% functions
        %%% infragravity band --------------------------------------------
        % spectral moments (X(F) must be a vector)
        cpM0inf   = trapz(stCP{i2n}.fInf(:,1), stCP{i2n}.SnnInf);     % [m^2] = variance (sigma^2)
        cpM1inf   = trapz(stCP{i2n}.fInf(:,1), stCP{i2n}.SnnInf .* stCP{i2n}.fInf);
        cpM2inf   = trapz(stCP{i2n}.fInf(:,1), stCP{i2n}.SnnInf .* stCP{i2n}.fInf.^2);

        cpNrmsInf = sqrt(cpM0inf);                 % rms value of the surface elev
        stCP{i2n}.HrmsInf = 2 * cpNrmsInf;         % rms wave height
        stCP{i2n}.Hm0Inf  = 4.004 * cpNrmsInf;     % 0-moment wave height (~H1/3)
        stCP{i2n}.Tm01Inf = cpM0inf ./ cpM1inf;       % mean wave period 
        stCP{i2n}.Tm02Inf = sqrt(cpM0inf ./ cpM2inf); % mean wave period: 0-up

        % peak frequency (in each column)
        [~, iMaxCPinf] = max(stCP{i2n}.SnnInf, [], 1);
        stCP{i2n}.msFpeakInf = stCP{i2n}.fInf(iMaxCPinf);
        
    end

    end    
%% depthAtt
    function [rSnnPSkp, rSnnMSkp, rKp, rEkz] = depthAtt(rSnnPS, rSnnMS, ...
                                                rF, nH, nHpt, nZbeg, nZend)
    %%% dynamic pressure is attenuated with increasing water depth and 
    %%% therefore Snn needs to be modified by a frequency dependent depth  
    %%% attenuation factor derived from linear wave theory; for a non-bottom   
    %%% mounted PT: Kp(f) = cosh(k*nHpt)/cosh(k*rH), where k is wavenumber
    %%% and Kp is referred to as the pressure response factor; reduces to
    %%% 1/cosh(k*rH) for a bottom mounted pressure transducer.
    
    % preallocate arrays
    [rK0, rKi, rK] = deal(NaN(nBurst/2 + 1, nPSD));
    
    % linearly interpolate burial depth
    nZ = [nZbeg nZend];
    nXq= 1 : 1/(nPSD -1) : 2;
    rZ = interp1(nZ, nXq);
    rZ = repmat(rZ, nBurst/2 + 1, 1);
    
    for i = 1 : nPSD
    % solve for all frequencies except 0 (avoid error with fzero )
        for j = 2 : length(rF(:,i))

            % calculate k from the dispersion function, not taking into 
            % account doppler shift effects (Jones and Monismith (2007))
            fun =  @(k) sqrt(nG * k * tanh(k * nH(1,i))) - ...
                        (rF(j,i) * 2 * pi);

            % initial estimation of k using Beji (2013) empirical formula
            rK0(j,i) = (2 * pi * rF(j,i))^2 / nG ;  % deep water wavenumber
            rKi(j,i) = ((rK0(j,i) * nH(1,i)) * ...
                        (1 + (rK0(j,i) * nH(1,i))^1.09 * ...
                        exp( -1*(1.55 + (1.3 * rK0(j,i) * nH(1,i)) + ...
                        (0.216 * (rK0(j,i) * nH(1,i))^2 )) ))) ...
                        / sqrt(tanh(rK0(j,i) * nH(1,i)) );

            % calculate k given inital guess for k
            rK(j,i) = fzero(fun, rKi(j,i));
            
        end
    end  
    
    % calculate Kp; if nHpt = 0 (bottom mounted PT), numerator = 1
    rKp      = cosh(rK * nHpt) ./ cosh(rK .* repmat(nH, nBurst/2+1, 1));
    rKp(1,:) = 1; % set zeroth term to one (zeroth term ignored above)
    
    % calculate depth attenuation factor by burial (e^kz)
    rEkz      = exp(rK .* rZ); 
    rEkz(1,:) = 1; % set zeroth term to one (zeroth term ignored above)
       
    % apply the transfer function to Snn
    rSnnPSkp    = rSnnPS .* (rEkz.^2 ./ rKp.^2);
    rSnnMSkp    = rSnnMS .* (rEkz.^2 ./ rKp.^2);
    
    end
% ----------------------------------------------------------------------- %
%% winSnn
    function [rSnnPSss, rSnnPSinf, rSnnMSinf, rSnnMSss, nFmaxKp] = ...
                         winSnn(rSnnPS, rSnnMS, rKp, rEkz, rF, nFminSS, ...
                                nFminIF, nFmaxSS, nHpt, nH)
    %%% window Snn with the user specified max and min frequencies

    % preallocate arrays
    [iFminSS, iFminIF] = deal(NaN(1,nPSD));
    
    % preallocate arrays
    [rSnnPSinf, rSnnPSss, rSnnMSinf, rSnnMSss] = deal(zeros(nBurst/2 + 1, nPSD));
    [nMinKp, nFmaxKp, nFmaxKpl] = deal(NaN(1,nPSD));
    
    for i = 1 : nPSD
        % find the windowing indices
        [~, iFminIF(1,i)] = min(abs(rF(:,i) - nFminIF));
        [~, iFminSS(1,i)] = min(abs(rF(:,i) - nFminSS));
        [~, iMinKp]       = min(abs(rF(:,i) - nFmaxSS));
        nMinKp(1,i)  = rKp(iMinKp,i); 
        nFmaxKp(1,i) = rF(iMinKp,i); 
    
        % maximum frequency that Kp should be applied to as predicted by 
        % linear wave theory for a submerged pressure sensor (nKmax is
        % the max linear wavenumber that should be applied to PT, assc w/
        % shortest detectable linear wave)
        nKmax = pi / (nH(1,i) - nHpt); 
        nFmaxKpl(1,i) = (1/(2*pi)) * sqrt(nG *nKmax*tanh(nKmax*nH(1,i)));
        
        % minimum Kmin predicted by linear wave theory 
        nMinKpl = cosh(nKmax * nHpt) / cosh(nKmax * nH(1,i));
    
        % check that Kmin is > than what is pred. by linear wave theory; if
        % not, set to a value 1% greater than that predicted by linear
        % wave theory value and calculate new nFmaxKp
        if ~(nMinKp(1,i) > nMinKpl)
            sprintf('ERROR: Kmin = %f < Kminl = %f for burst %d', ...
                    nMinKp(1,i), nMinKpl, nPSD)
            nMinKp(1,i) = nMinKpl*1.01;
            [~, iMinKp] = min(abs(rKp(:,i)-nMinKp(1,i)));
            nFmaxKp(1,i) = rF(iMinKp,i);
            sprintf('New Kmin set to %f', nMinKp(1,i))
                
            % check that FmaxKp is < what is pred. by linear wave theory 
            if ~(nFmaxKp(1,i) < nFmaxKpl(1,i))
                sprintf('ERROR: FmaxKp = %f > FmaxKpl = %f', ...
                    nFmaxKp(1,i), nFmaxKpl(1,i))
            else
                sprintf('New FmaxKp = %f < FmaxKpl = %f', ...
                    nFmaxKp(1,i), nFmaxKpl(1,i))
            end
        end
    
        % recalculate rSnnKp after setting Kp=kpmin > FmaxKp and window to 
        % the incident wind wave cutoff frequencies (nFmin and nFmaxKp*) &
        % the infragravity wave cutoff frequencies (nFinf and nFmin)
        % (*NOTE: I can't see why I need to extrapolate the tail when I 
        % have no idea what the real spectrum looks like at higher
        % frequencies, therefore I just window to the ~minimum wavelength
        % that the PT can detect!)
        rKp(iMinKp(end):end,i) = nMinKp(1,i);
        
        rSnnPSss(iFminSS(1,i):iMinKp,i) = rSnnPS(iFminSS(1,i):iMinKp,i) .* ...
         (rEkz(iFminSS(1,i):iMinKp,i).^2 ./ rKp(iFminSS(1,i):iMinKp,i).^2);
        
        rSnnMSss(iFminSS(1,i):iMinKp,i) = rSnnMS(iFminSS(1,i):iMinKp,i) .* ...
         (rEkz(iFminSS(1,i):iMinKp,i).^2 ./ rKp(iFminSS(1,i):iMinKp,i).^2);
      
        rSnnPSinf(iFminIF(1,i):iFminSS(1,i),i) = ...
         rSnnPS(iFminIF(1,i):iFminSS(1,i),i) .* ...
         (rEkz(iFminIF(1,i):iFminSS,i).^2 ./ rKp(iFminIF(1,i):iFminSS(1,i),i).^2);
        
        rSnnMSinf(iFminIF(1,i):iFminSS(1,i),i) = ...
         rSnnMS(iFminIF(1,i):iFminSS(1,i),i) .* ...
         (rEkz(iFminIF(1,i):iFminSS,i).^2 ./ rKp(iFminIF(1,i):iFminSS(1,i),i).^2);
        
    end
    end
% ----------------------------------------------------------------------- %
%% timeDomain
    function [rFSElpf, rFSEbpf] = timeDomain(rPd, dtPb, nH)
    %%% calculate wave height as 4x the standard deviation in the time
    %%% domain; note, I don't apply the pressure response factor here since
    %%% it is a function of frequency (k) - also, it approaches unity for
    %%% low frequency waves. We filter water surface elevation using a 
    %%% zero-phase lpf and bpf (there may be transients at the beginning 
    %%% and end if you use an IIR; spot check as needed).
    
    %%% 1/20/18 : I don't use the LPF or BPF wave statistics because they
    %%% introduce too much leakage from other frequencies - retaining this
    %%% code for future needs using the signal processing toolbox 
    
    % free surface elevation (mean removed), minus the pressure correction 
    % and burial factors
    rFSE = rPd / (nRho * nG); 
        
    % design LPF filters -----------------------------------------      
    % passband frequency = start of the pass band (end of infr. band) [Hz]
    % stopband frequency = end of the transition width [Hz]
    % passband ripple = amount of ripple allowed in the pass band [db]
    % stopband attenuation = attenuation over transition width [db]
    % NOTE: use info(dLPF) to get the filter specs and 
    % grpdel(dLPF, 2048, nFsamp) to plot group delay as a function of 
    % frequency
    
    % finite impulse response filters (usally constant filter delay - 
    % linear phase filters)
    dHamLPF = designfilt('lowpassfir', 'FilterOrder', 200, ...
                         'CutoffFrequency', nFminSS, 'SampleRate', nFsamp);
                     
%     dKaiserLPF = designfilt('lowpassfir', 'PassbandFrequency', 0.005, ...
%                'StopbandFrequency', 0.04, 'PassbandRipple', 0.5, ...
%                'StopbandAttenuation', 30, 'SampleRate', nFsamp,  ...
%                'DesignMethod', 'kaiserwin');

    % infinite impulse response filter (frequency dependent filter delay -
    % nonlinear phase filters)                 
%     dButterLPF = designfilt('lowpassiir', 'PassbandFrequency', 0.005 , ...
%                         'StopbandFrequency', 0.04, 'PassbandRipple', 1, ...
%                         'StopbandAttenuation', 30, 'SampleRate', nFsamp,...
%                         'DesignMethod', 'butter');
    
    % NOTE: this introduced transients using the zerophase filter and the
    % frequency dependent group delay (which peaks at 0.05 Hz) --> NOT GOOD
%     dEllipticLPF = designfilt('lowpassiir', 'PassbandFrequency', 0.05, ...
%                       'FilterOrder', 7, 'PassbandRipple', 1, ...
%                       'StopbandAttenuation', 60, 'SampleRate', nFsamp);
                  
    % order of different LPFs (for sensitivity modeling)
%     nNButter = filtord(dButterLPF);
%     nNHam    = filtord(dHamLPF);
%     nNKaiser = filtord(dKaiserLPF);
%     nNEllip  = filtord(dEllipticLPF);
    
    % plot filters (for debugging)
    fvtool(dHamLPF)
    
    % design bandpass filter for sea-swell waves
    % NOTE: this is probably not the best filter for determining sea-swell
    % (high pass FIR is a better approximation) but both methods contain
    % spikes where there are large changes in slope (spikes are larger in 
    % the HPF); as usual, the freq domain estimates of sea swell are best
    dBPF   = designfilt('bandpassiir', 'FilterOrder', 20, ...
               'HalfPowerFrequency1', nFminSS, 'HalfPowerFrequency2', 0.25,...
               'SampleRate', nFsamp);
           
%     dEqHPF = designfilt('highpassfir', 'StopbandFrequency', 0.005, ...
%              'PassbandFrequency', 0.05,'StopbandAttenuation', 65, ...
%              'PassbandRipple', 0.5, 'SampleRate', nFsamp, 'DesignMethod',...
%              'equiripple');

    % plot filters (for debugging)
    %fvtool(dEqHPF)
    fvtool(dBPF)
    
    % preallocate arrays
    [rFSElpf, rFSEbpf] = deal(NaN(size(rFSE)));
    fTS = NaN(size(nPSD));
    
    for i = 1 : nPSD
        % LPF
        rFSElpf(:,i) = filtfilt(dHamLPF, rFSE(:,i)); 

        % BPF
        rFSEbpf(:,i) = filtfilt(dBPF, rFSE(:,i)); 
        
        if bPlot
            % added for ADV sensitivity testing
                
            % check if figure directory has been created 
            if exist(sDirName, 'dir') == 0      
                mkdir(sprintf('%s', sDirName))  % save figs to new folder
            end
    
            % plot LPF and BPF
            fTS(i) = figure; plot(rFSE(:,i))
                hold on
                plot(rFSElpf(:,i))
                hold on
                %plot(rFSEbpf(:,i))
                plotFancyAxis
                xlabel('Sample #')
                ylabel('free surface elevation (m)')
                legend('no filter', 'infragravity')
                %legend('no filter', 'infragravity', 'sea-swell')
                plotFancyAxis
                title(sprintf('%s, Burst %d, h = %f', dtPb(1,i), i, nH(1,i)))
                saveas(fTS(i), [pwd sprintf('/%s/%s-burst%d-TS.fig', ...
                       sDirName, sFileName, i)]);
            close(fTS(i))
        end     
    end  
    end
end