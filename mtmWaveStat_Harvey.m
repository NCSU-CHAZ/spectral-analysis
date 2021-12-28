function [waveStat] = mtmWaveStat_Harvey(rSnn, rF, rCIlow, rCIup, rFtest, ...
                        dtPb, rPb, nFminIF, nFminSS, nFmaxSS, bPlot, sDirName, ...
                        nHpt, sFileName, nZbeg, nZend)
% --------------------------mtmWaveStat_Harvey-----------------------------%
% Purpose: This function generates wave statistics from bottom mounted 
% pressure transducers and an ADV deployed durring Hurricane Harvey using 
% the Thomson multitaper spectral estimator (in R for maxfunctionality). 
% The wave statistics and spectral modifications are adapted from Dean &
% Dalyrymple [1991] and Jones and Monismith (2007). This algorithm differs 
% from szWaveStat_Harvey in that it does not run psd.f but simply
% calculates and plots statistics from the input power spectra, generated 
% using spec.mtm in R (Rahim and Thomson, 2017).
% 
% Inputs:
%       - rSnn:      matrix of elev. var. spectra from spec.mtm (m^2/Hz)
%       - rF:        ... frequencies (Hz) ...
%       - rCIlow:    ... lower confidence intervals from jacknifing...
%       - rCIup:     ... upper confidence interval from jacknifing...
%       - dtPb:      vector of datetime values for windowed pressure ts
%       - rPb:       matrix of windowed pressure time series (Pa)
%       - rFtest:    ... fTest outputs as a function of frequency
%       - nFminIF:   min frequency used for windowing infragravity band
%       - nFminSS:   min frequency used for windowing the sea-swell band 
%       - nFmaxSS:   max frequency that the PTs can detect 
%                    (recommended value = 0.7 Hz)
%       - bPlot:     boolean for plotting
%       - sDirName:  string for name of new directory containing figures 
%                    ([] for bPlot = false)
%       - nHpt:      height of pressure transducer off (above) the bed (m)
%                    (this is typically zero, do not include negative
%                    numbers...burial is accounted for using poroelastic
%                    theory)
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
%       - Snn:      water surface elevation variance spectra [m^2/Hz]
%       - SnnNorm:  normalized Snn [1/Hz]
%       - f:        frequency array [Hz]
%       - L0:       deep water wavelength [m]
%
% CODE FOR MULTITAPER ANALYSIS IN MATLAB:
%
%    % use iteratively to identify the ideal number of tapers per bandwidth
%    [e,v] = dpss(8192,2,10);
%    figure; stem(1:length(v),v,'filled')
%
%    % note eigenvalues are time series specific and this spectra doesn't
%      need to be multiplied by 2 like the R version; CIs are based on dofs
%    [pxx, f2, pxxc] = pmtm(detrend(stCPfin{1,1}.win(:,1)), ...
%                           e(:,1:idx), v(1:idx), ...
%                           length(stCPfin{1,1}.win(:,1)), 16, ...
%                           'DropLastTaper',false, 'onesided', ...
%                           'ConfidenceLevel', 0.95);
%    figure; loglog(f2, pxx, '-k', f2, pxxc, '-r')
%
% SEE ALSO: szWaveStat_Harvey, multitaper (R)
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       02/20/17        KAnarde             Original code
%       10/21/17        " "                 Updates for Harvey
%       02/01/18        " "                 MTM modification
%       09/28/18        " "                 Fixed confidence limit problem
%
%------------------------------preamble-----------------------------------%

disp('-----------------------------------------------------------')
disp('-----------------------mtmWaveStat-------------------------')
disp('-----------------------------------------------------------') 

% constants and conversion parameters
nRho   = 1025;   % seawater density [kg/m3]
nG     = 9.81;   % gravitational constant [m/s2]

% define dimensions for looping
nBurst = size(rSnn, 1);  % length of the PSD
nPSD   = size(rSnn, 2);  % number of PSDs
nN     = size(rPb, 1);   % number of samples

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~call subroutines~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% 1) calculate mean water depth
[rPd, rH] = dynPress(rPb, nHpt);

% 2) account for attenuation of the pressure signal by burial & water depth
[rSnnKp, rKp, rEkz] = depthAtt(rSnn, rF, rH, nHpt, nZbeg, nZend);

% 6) window to min and max frequency
[rSnnSS, rSnnInf, iFminIF] = winSnn(rSnn, rKp, rEkz, rF, nFminSS, nFminIF, ...
                              nFmaxSS, nHpt, rH);
    
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~wave statistics~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%%% infragravity band ----------------------------------------------------
% spectral moments (X(F) must be a vector)
rM0inf   = trapz(rF(:,1), rSnnInf);           % [m^2] = variance (sigma^2)
rM1inf   = trapz(rF(:,1), rSnnInf .* rF);
rM2inf   = trapz(rF(:,1), rSnnInf .* rF.^2);

rNrmsInf = sqrt(rM0inf);                      % std deviation (sigma)
waveStat(1).HrmsInf = 2 * rNrmsInf;           % rms wave height
waveStat(1).Hm0Inf  = 4.004 * rNrmsInf;       % 0-moment wave height (H1/3)
waveStat(1).Tm01Inf = rM0inf ./ rM1inf;       % mean wave period 
waveStat(1).Tm02Inf = sqrt(rM0inf ./ rM2inf); % mean wave period: 0-up

% peak frequency (in each column)
[~, iMaxInf] = max(rSnnInf, [], 1);
waveStat(1).FpeakInf = rF(iMaxInf);

% deep water wavelength
waveStat(1).L0inf = (nG * waveStat(1).Tm01Inf.^2) / (2*pi);

% normalize spectra for comparison of spectral shapes (area=1)
waveStat(1).SnnNormInf = rSnnInf ./ repmat(rM0inf, nBurst, 1);

%%% sea-swell band -------------------------------------------------------
% spectral moments 
rM0ss   = trapz(rF(:,1), rSnnSS);           % [m^2] = variance (sigma^2)
rM1ss   = trapz(rF(:,1), rSnnSS .* rF);
rM2ss   = trapz(rF(:,1), rSnnSS .* rF.^2);

rNrmsSS = sqrt(rM0ss);                      % std deviation (sigma)
waveStat(1).HrmsSS = 2 * rNrmsSS;           % rms wave height
waveStat(1).Hm0SS  = 4.004 * rNrmsSS;       % 0-moment wave height (H1/3)
waveStat(1).Tm01SS = rM0ss ./ rM1ss;        % mean wave period 
waveStat(1).Tm02SS = sqrt(rM0ss ./ rM2ss);  % mean wave period: 0-up

% peak frequency (in each column)
[~, iMaxSS] = max(rSnnSS, [], 1);
waveStat(1).FpeakSS = rF(iMaxSS);

% deep water wavelength
waveStat(1).L0ss = (nG * waveStat(1).Tm01SS.^2) / (2*pi);

% normalize spectra for comparison of spectral shapes (area=1)
waveStat(1).SnnNormSS = rSnnSS ./ repmat(rM0ss, nBurst, 1);

% create a new rSnnKp 
rSnnKp = zeros(size(rSnnInf));
rSnnKp(1:iFminIF(1,1),:) = rSnn(1:iFminIF(1,1),:);
rSnnKp = rSnnKp + rSnnInf + rSnnSS;

%%% save other variables to structure ------------------------------------
waveStat(1).SnnInf= rSnnInf;
waveStat(1).SnnSS = rSnnSS;
waveStat(1).Snn   = rSnn;
waveStat(1).SnnKp = rSnnKp;
waveStat(1).f     = rF;
waveStat(1).h     = rH;
waveStat(1).fTest = rFtest;
waveStat(1).CIlow = rCIlow; 
waveStat(1).CIup  = rCIup;
waveStat(1).dtPb  = dtPb;
waveStat(1).Pb    = rPb;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
if bPlot
    
    % preallocate figure array
    f = NaN(size(nPSD));
    
    % find max of Y value for plots
    nYmax = max(cat(1, rSnnInf(:), rSnnSS(:)));
    
    % check if figure directory has been created 
    if exist(sDirName, 'dir') == 0      % 7 specifies a folder
        mkdir(sprintf('%s', sDirName))  % save figs to new folder
    end

    nScrsz = get(groot, 'ScreenSize');
    
    for n = 1 : nPSD
        % plot Snn spectra
        f(n) = figure;
        
        subplot(1,2,1)
        plot(rPd(:,n)/(nRho*nG));
        xlabel('sample number'); 
        ylabel('free surface elevation [m]')
        plotFancyAxis 
        
        subplot(1,2,2)
              loglog(rF(:,n), rSnnSS(:,n), '-.', ...
                     rF(:,n), rSnnInf(:,n), '-.', ...
                     rF(:,n), rSnnKp(:,n), ...
                     rF(:,n), rSnn(:,n), '--');
        legend('Snn_{SS}', 'Snn_{inf}', 'Snn_{Kp}', 'Snn')
        xlabel('frequency (Hz)'); 
        ylabel('S_{\eta\eta} [m^2/Hz]')
        set(gca, 'Xlim', [0 1] , 'Ylim', [0 nYmax]) 
        plotFancyAxis    
        
        % make figure large and save
        set(f(n), 'Position', [1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);
        saveas(f(n), [pwd sprintf('/%s/%s-burst%d.fig', sDirName, ...
                      sFileName, n)]);
        close(f(n))
    end
end
        
% ----------------------------------------------------------------------- %
%% dynPress
    function [rPd, nH] = dynPress(rPb, nHpt)
    %%% dynamic pressure is found by subtracting the mean (hydrostatic)
    %%% pressure from the sea (gauge) pressure (absolute - atmospheric)
    %%% hydrostatic approximation: n(-h)= (P(-h)-rho*g*h)/(rho*g*Kp(-h)) 
    
    % calculate dynamic pressure along each column for the burst samples
    rPd = rPb - repmat(mean(rPb,1), nN, 1);

    % calculate mean water level (m), accounting for nHpt (height of PT
    % above the bottom)
    nH = ( mean(rPb,1) / (nRho * nG) ) + nHpt; 
    
    end
% ----------------------------------------------------------------------- %   
%% depthAtt
    function [rSnnKp, rKp, rEkz] = depthAtt(rSnn, rF, nH, nHpt, nZbeg, nZend)
    %%% dynamic pressure is attenuated with increasing water depth and 
    %%% therefore Snn needs to be modified by a frequency dependent depth  
    %%% attenuation factor derived from linear wave theory; for a bottom   
    %%% mounted PT: Kp(f) = cosh(k*nHpt)/cosh(k*rH), where k is wavenumber
    %%% and Kp is referred to as the pressure response factor; reduces to
    %%% 1/cosh(k*rH) for a bottom mounted pressure transducer.
    
    % preallocate arrays
    [rK0, rKi, rK] = deal(NaN(nBurst, nPSD));
    
    % linearly interpolate burial depth
    nZ  = [nZbeg nZend];
    nXq = 1 : 1/(nPSD-1) : 2;
    rZ  = interp1(nZ, nXq);
    rZ  = repmat(rZ, nBurst, 1);
    
    for i = 1 : nPSD
    % solve for all frequencies except 0 (avoid error with fzero)
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
    rKp      = cosh(rK * nHpt) ./ cosh(rK .* repmat(nH, nBurst, 1));
    rKp(1,:) = 1; % set zeroth term to one (zeroth term ignored above)
    
    % calculate depth attenuation factor by burial (e^kz)
    rEkz      = exp(rK .* rZ); 
    rEkz(1,:) = 1; % set zeroth term to one (zeroth term ignored above)
       
    % apply the transfer function to Snn
    rSnnKp    = rSnn .* (rEkz.^2 ./ rKp.^2);

    end
% ----------------------------------------------------------------------- %
%% winSnn
    function [rSnnSS, rSnnInf, iFminIF] = winSnn(rSnn, rKp, rEkz, rF, ...
                                       nFminSS, nFminIF, nFmaxSS, nHpt, nH)
    %%% window Snn with the user specified max and min frequencies

    % preallocate arrays
    [iFminSS, iFminIF] = deal(NaN(1,nPSD));
    
    % preallocate arrays
    [rSnnInf, rSnnSS] = deal(zeros(nBurst, nPSD));
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
        
        rSnnSS(iFminSS(1,i)+1:iMinKp,i) = rSnn(iFminSS(1,i)+1:iMinKp,i) .* ...
         (rEkz(iFminSS(1,i)+1:iMinKp,i).^2 ./ rKp(iFminSS(1,i)+1:iMinKp,i).^2);
        
        rSnnInf(iFminIF(1,i):iFminSS(1,i),i) = ...
         rSnn(iFminIF(1,i):iFminSS(1,i),i) .* ...
         (rEkz(iFminIF(1,i):iFminSS,i).^2 ./ rKp(iFminIF(1,i):iFminSS(1,i),i).^2);
            
    end
    end
% ----------------------------------------------------------------------- %
end