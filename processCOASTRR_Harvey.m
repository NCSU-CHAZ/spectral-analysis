%function [stPOD, stVECT, stBP, stVECTwin, stWaveStat, stADVWaveStat] = ...
%                                                   processCOASTRR_Harvey()        
% -----------------------processCOASTRR_Harvey----------------------------%
% Purpose: This script post-processes data from pressure transducers 
% (RBR solo + solo D|wave), TCMs, and an ADV (Vector) collected during 
% Hurricane Harvey which made landfall in Rockport, Texas on August 25,
% 2017
% 
% SEE ALSO: plotBARis_Test8, processBARis_Test8, processBARis_Test9 (flume)
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       09/04/17        KAnarde             Original code
%
%------------------------------preamble-----------------------------------%

% start with a clean slate
close all;                    % close all previous figures
%clear;                        % clear workspace variables
clc;                          % clear the command window
%dbstop if error;

%---------------user input for read scripts and subroutines---------------%

% booleans 
bRemoveOld = false; % re-process data (isolate time series, unit conversion, 
                    % calculate sea-pressure, despiking, filtering data)

% specify input file folders for read scripts
sPSfldr    = 'PT';
sTCMfldr   = 'TCM';
sVECTfldr  = 'ADV';
sATMfldr   = 'ATM';
sUSGSfldr  = 'USGS';
sNDBCfldr  = 'NDBC';
sRTKfldr   = {'RTK/BeachProfiles/FolletsWest-Center' ...
              'RTK/BeachProfiles/FolletsWest-East'   ...
              'RTK/BeachProfiles/FolletsWest-West'   ...
              'RTK/BeachProfiles/Matagorda-Center'   ...
              'RTK/BeachProfiles/Matagorda-East'     ...
              'RTK/BeachProfiles/Matagorda-West'};

% specify output .MAT file names for read scripts
sPODSaveOut  = 'RRU';
sVECTSaveOut = 'VECT';
sATMSaveOut  = 'ATM';
sRTKSaveOut  = 'RTK';
sUSGSSaveOut = 'USGS';
sNDBCSaveOut = 'NDBC';
sPPSaveOut   = 'PostProcessed';

% inputs for badPress (corrects RRU-1 pressure time series for offset)
sSTbuck  = '22-Aug-2017 23:01:00'; % start/end of bucket test
sETbuck  = '22-Aug-2017 23:10:00';

% zero point vector for TCMs (based on Test 9 data for 15 cm tether); note
% this needs to be updated with Test 10 data since I used a 10 cm tether
% (however, the correction is minor)
nZpt       = [0.0269 0.0148 0.9874]; % XYZ for 15 cm tether

% sampling frequencing and RRU# for all instruments (RRU-1 [HI backshore], 
% RRU-2 [HI backbarrier], RRU-6 [FI backbarrier], RRU-7 [FI midbarrier])
nFsamp     = [16, 2, 2, 16];
nRRUid     = [1, 2, 6, 7];

% specify time-series when when we know the instruments were subjected to
% wave action (organized by RRU) during the first landfall
% NORMAL                     %:31
sSTime1  = {'25-Aug-2017 03:07:14', '25-Aug-2017 03:07:14', ...
           '25-Aug-2017 03:07:14', '25-Aug-2017 11:00:54'}; 
sETime1  = {'26-Aug-2017 07:37:49', '26-Aug-2017 07:37:49', ...
           '26-Aug-2017 14:31:38', '26-Aug-2017 14:31:38'}; 
% JUST LANDFALL
% sSTime1  = {'25-Aug-2017 17:15:17', '25-Aug-2017 17:15:17', ...
%            '25-Aug-2017 11:00:54', '25-Aug-2017 11:00:54'};
% sETime1  = {'25-Aug-2017 20:00:00', '25-Aug-2017 20:00:00', ...
%            '26-Aug-2017 14:31:38', '26-Aug-2017 14:31:38'};        
% FOR SPECTRA SENSITIVITY MODELING
% sSTime1  = {'25-Aug-2017 10:07:14', '25-Aug-2017 03:07:14', ...
%            '25-Aug-2017 11:00:54', '25-Aug-2017 11:00:54'};
% sETime1  = {'25-Aug-2017 17:15:17', '26-Aug-2017 07:37:49', ...
%             '26-Aug-2017 14:31:38', '26-Aug-2017 14:31:38'};
% including second landfall (NOTE: RRU-1 not inundated during 2nd landfall; 
% keep as first landfall timing for simplicity)
% sSTime1  = {'25-Aug-2017 03:07:21', '25-Aug-2017 03:07:21', ...
%             '25-Aug-2017 11:00:54', '25-Aug-2017 11:00:54'};  
% sETime1  = {'26-Aug-2017 07:37:49', '30-Aug-2017 02:46:39', ...
%             '30-Aug-2017 02:46:39', '30-Aug-2017 02:46:39'};

% instrument elevations (m NAVD88, post-storm)
nVECTele   = -0.746;                     % bed elevation at ADV
nVECTPTele = nVECTele + 0.275 + 0.425;   % ADV (PT sensor 0.7 m above bed)
nPTele     = [0.9315 0.260 0.190 0.500]; % pressure transducer (top)
nPODele    = [1.1490 0.490 0.410 0.710]; % top of casing

% conversion parameters and constants
nDbar2Pa = 10000; %[Pa/dbar]
nRho     = 1025;  %[kg/m3]
nG       = 9.81;  %[m/s2]

% specify custom calibration curve (X = pitch[deg], Y = speed [cm/s])
rTCMcalX = [0; 1; 2; 3; 5; 12.5; 20; 27; 34; 40; 46; 52; 57; 68; 75; 90];
rTCMcalY = [0; 0.5; 5; 10; 20; 30; 40; 50; 60; 70; 80; 90; 100; 125; 150; 200];
rTCMcalY = rTCMcalY/100; % convert cm/s to m/s

% index for ADV start (due to the barometric time series)
iSeapress = 345568; % index for start of barometric time series

%-------------------------------------------------------------------------%
%---------------------------call read scripts-----------------------------%
%-------------------------------------------------------------------------%

% check to see if readPOD has been executed
if exist(sprintf('%s.mat', sPODSaveOut),'file')
    load(sprintf('%s.mat', sPODSaveOut));
else
    % call readPOD
    [stPOD] = readPOD(sPSfldr,sTCMfldr);
    save(sprintf('%s.mat', sPODSaveOut), 'stPOD', '-v7.3')
end

% check to see if readVECT has been executed
if exist(sprintf('%s.mat', sVECTSaveOut),'file')
    load(sprintf('%s.mat', sVECTSaveOut));
else
    % call readVECT (and read .SEN file)
    [stVECT] = readVECT(sVECTfldr, 'true');
    save(sprintf('%s.mat', sVECTSaveOut), 'stVECT', '-v7.3')
end

% check to see if readATM has been executed
if exist(sprintf('%s.mat', sATMSaveOut),'file')
    load(sprintf('%s.mat', sATMSaveOut));
else
    % call readATM
    [stATM] = readATM(sATMfldr);
    save(sprintf('%s.mat', sATMSaveOut), 'stATM', '-v7.3')
end

% check to see if readRTK has been executed
if exist(sprintf('%s.mat', sRTKSaveOut),'file')
    load(sprintf('%s.mat', sRTKSaveOut));
else
    % call readRTK for each site (don't save decimal degrees as text file)
    stBP.FWcent = readRTK(char(sRTKfldr(1)), 'false');
    stBP.FWeast = readRTK(char(sRTKfldr(2)), 'false');
    stBP.FWwest = readRTK(char(sRTKfldr(3)), 'false');
    stBP.MATAcent = readRTK(char(sRTKfldr(4)), 'false');
    stBP.MATAeast = readRTK(char(sRTKfldr(5)), 'false');
    stBP.MATAwest = readRTK(char(sRTKfldr(6)), 'false');
    % save
    save(sprintf('%s.mat', sRTKSaveOut), 'stBP', '-v7.3')
end

% check to see if USGS data has been read
if exist(sprintf('%s.mat', sUSGSSaveOut),'file')
    load(sprintf('%s.mat', sUSGSSaveOut));
else
    % call readUSGS
    [stUSGS] = readUSGS(sUSGSfldr);
    save(sprintf('%s.mat', sUSGSSaveOut), 'stUSGS', '-v7.3')
end

% check to see if NDBC data has been read
if exist(sprintf('%s.mat', sNDBCSaveOut),'file')
    load(sprintf('%s.mat', sNDBCSaveOut));
else
    % call readNDBC
    [stNDBC] = readNDBC(sNDBCfldr);
    save(sprintf('%s.mat', sNDBCSaveOut), 'stNDBC', '-v7.3')
end

%-------------------------------------------------------------------------%
%----------------------------call subroutines-----------------------------%
%-------------------------------------------------------------------------%

if bRemoveOld
    disp('---------------------------------')
    disp('Removing previously saved output:')
    disp('---------------------------------')
    unix(sprintf('rm %s.mat', sPPSaveOut));
    
    % convert atmospheric pressure units
    [stATM] = convUnits(stATM);

    % calculate mean difference in bucket test between RRU-1 & 2 -->
    % correct RRU-1 pressure
    [stPOD] = badPress(stPOD, sSTbuck, sETbuck); 
          
    % correct for local variations in atmospheric pressure
    [stPOD, stVECT] = atmPress(stPOD, stATM, stVECT);  
    
    % calculate water levels and stillwater elevation
    [stPOD, stVECT] = waterlevelTS(stPOD, stVECT, nDbar2Pa, ...
                                            nRho, nG, nPTele, nPODele, ...
                                            nVECTPTele, nVECTele);
                                        
    % window time series to periods when pressure transducers subjected to
    % wave action during first landfall
    [iMinPS1, iMaxPS1, iMinTCM1, iMaxTCM1, iMinVECT, iMaxVECT] = ...
                                 indexTS(stPOD, stVECT, sSTime1, sETime1); 
    
    % wave spectra and statistics
    [stWaveStat, stADVWaveStat] = waveStats(stPOD, stMTMps, ...
                                                stVECT, iMinPS1, iMaxPS1);
                                        
    % calculate TCM tilt from accelerometer components
    [stPOD] = tiltTCM(stPOD, nZpt)
    
    % window TCM and PT time series for comparison within bursts
    [stTCMwin, stPSwin] = windowTS(stPOD, iMinTCM1, iMaxTCM1, iMinPS1, iMaxPS1) 
 
    % create cross-spectra 
    [stCxTiltXPT, stCxTiltYPT, stCxSpeedPT] = cross(stPSwin, stTCMwin);
    
    % lowpass filter the TCM windows
    [stTCMwin] = lpfTCMwin(stTCMwin, nZpt);
    
    % calculate wave statistics for NDBC wave spectra
    [stNDBC] = ndbcWaveStats(stNDBC);
    
    % despike ADV data
    [stVECT] = despikeVel(stVECT, nFsamp, iSeapress, iMinVECT);
    
    % evaluate frequency dependence of wave dissipation
    waveDissipation(stMTMWaveStat0pt003);
     
    % save new variables
    save(sprintf('%s.mat', sPPSaveOut), '-v7.3')    
else    
    load(sprintf('%s.mat', sPPSaveOut)) 
end

disp('-----------------------------------------------------------')
disp('----------------DONE WITH DATA PROCESSING------------------')
disp('-----------------------------------------------------------')

    %---------------------------------------------------------------------%
    %---------------------------subroutines-------------------------------%
    %---------------------------------------------------------------------%

%% convUnits   
    function [stATM] = convUnits(stATM)

    nKpa2Dbar = 0.1; % conversion parameter
        
    % convert atmospheric pressure data (kPa) to dbar
    for i=1 : length(stATM) 
            stATM(i).baropress = stATM(i).baropress * nKpa2Dbar;
            stATM(i).units     = 'dbar';
    end
        
    end
    % ------------------------------------------------------------------- %
%% badPress   
    function [stPOD] = badPress(stPOD, sSTbuck, sETbuck)
    %%% The pressure transducer at RRU-1 recorded at a lower (and
    %%% incorrect) pressure throughout the deployment. We corrected for the
    %%% erroneously low pressure reading through a linear shift, calculated
    %%% by differencing the pressure recorded during the bucket test on
    %%% August 22, 2017.
    
    [iMin, iMax] = deal(NaN(2,1));
  
    % calculate the mean difference in pressure between RRU-1 and RRU-2
    % during the bucket test: 23:01-23:10
    iMin(1) = find(stPOD(1).PSdt == sSTbuck);
    iMax(1) = find(stPOD(1).PSdt == sETbuck);
    iMin(2) = find(stPOD(2).PSdt == sSTbuck);
    iMax(2) = find(stPOD(2).PSdt == sETbuck);
    
    % downsample 16 Hz time series at RRU-1
    rPressDS= downsample(stPOD(1).pressure(iMin(1):iMax(1)), 8);
    nOffset = mean(stPOD(2).pressure(iMin(2):iMax(2)) - rPressDS);

    % add offset
    stPOD(1).pressure = stPOD(1).pressure + nOffset;
        
    end
    % ------------------------------------------------------------------- %
%% atmPress    
    function [stPOD, stVECT] = atmPress(stPOD, stATM, stVECT)
    %%% The RBRs and ADV record absolute pressure. The following code uses 
    %%% barometric pressure, from a subaerially mounted HOBO (@ RRU-6 & 7), 
    %%% to calculate sea/gauge pressure (recorded every 30 seconds). It 
    %%% also works for NOAA barometric pressure data (RRU-1 & 2, recorded
    %%% every 6 minutes). The code only applies to when the barometric 
    %%% pressure time series is shorter than the PT time series 
    %%% (needs to be made more robust in the future). Also, .csv files are
    %%% terrible with the datetime format; typically have to "customize"
    %%% field in excel.
    
    % RBRs --------------------------------------------------
    % preallocate arrays
    [iStart, iEnd, nTpt, nTatm] = deal(NaN(size(stPOD,2),1));                                     

    for i=1 : length(stPOD)
        % preallocate array
        [stPOD(i).seapress] = deal(NaN(size(stPOD(i).pressure)));
        
        % calculate the sampling period (s) of the PT
        nTpt(i) = mean(seconds(diff(stPOD(i).PSdt)));
        
        % calculate the sampling period (s) of the baropress time series
        nTatm(i) = round(mean(seconds(diff(stATM(i).dt))));
        
        % find index values for baropress start/end time within RBR time
        % series
        iStart(i) = find(stPOD(i).PSdt == stATM(i).dt(1));
        iEnd(i)   = find(stPOD(i).PSdt == stATM(i).dt(end));
        
        % for RBRs, we need to interpolate the baropress time series to 2 
        % or 16 Hz; for some reason the interpolated time series is always
        % one data point short (interpolates up to the last point) for the
        % 16 Hz time series
        rBaropressInt = interp1(stATM(i).baropress, ...
                   (1 : (nTpt(i)/nTatm(i)) : length(stATM(i).baropress)))';
        
        % for 2 Hz RBRs
        if rem(nTpt(i),0.5)==0 % differentiates sampling frequencies
            % calculate sea-pressure (gauge)
            stPOD(i).seapress(iStart(i):iEnd(i),1) = ...
                    stPOD(i).pressure(iStart(i):iEnd(i)) - rBaropressInt;
        
            % save interpolated barometric pressure
            stPOD(i).baropress = NaN(size(stPOD(i).pressure));
            stPOD(i).baropress(iStart(i):iEnd(i)) = rBaropressInt;
            
        % for 16 Hz RBRs
        else
            % calculate sea-pressure (gauge)
            stPOD(i).seapress(iStart(i):iEnd(i)-1,1) = ...
                    stPOD(i).pressure(iStart(i):iEnd(i)-1) - rBaropressInt;
        
            % save interpolated barometric pressure
            stPOD(i).baropress = NaN(size(stPOD(i).pressure));
            stPOD(i).baropress(iStart(i):iEnd(i)-1) = rBaropressInt;
        end
    end
    
    % ADV --------------------------------------------------                                   

    % calculate the sampling period 
    nTptADV = mean(seconds(diff(stVECT.dt)));
 
    % find index values for baropress start time within ADV time series
    % using the barometric time series from Follets
    iStartADV = find(stVECT.dt == stATM(4).dt(1));
    iEndATM   = find(stATM(4).dt == stVECT.dt(end));
    
    % since the ADV shut off early, cutoff baropress time series
    rBaropress = stATM(4).baropress(1:iEndATM);  
    
    % for the ADV, we interpolate the baropress time series to 16 Hz
    rBaropressInt = interp1(rBaropress, ...
                           (1 : (nTptADV/nTatm(4)) : length(rBaropress)))';
        
    % calculate sea-pressure (gauge)
    % NOTE: the offset for the ADV is what I'm guessing is the "calibration
    % in air" prior to deployment
    stVECT.seapress = NaN(size(stVECT.pressure));
    stVECT.seapress(iStartADV:end,1) = ...
                    stVECT.pressure(iStartADV:end) - rBaropressInt + 10.14;
        
    % save interpolated barometric pressure
    stVECT.baropress = NaN(size(stVECT.pressure));
    stVECT.baropress(iStartADV:end) = rBaropressInt;
    
    end
    % ------------------------------------------------------------------- %
%% waterlevelTS  
    function [stPOD, stVECT] = waterlevelTS(stPOD, stVECT, nDbar2Pa, ...
                                            nRho, nG, nPTele, nPODele, ...
                                            nVECTPTele, nVECTele)
    %%% Calculates the water surface elevation (m WGS84) and stillwater 
    %%% elevation (5 minute running mean of WSE) for the RRUs and ADV. 
    %%% NOTE: the USGS calculates SWE by LPF'ing the dynamic pressure using 
    %%% a 4th order butterworth filter with a cutoff frequency of 1 min 
    %%% (forward and backward)--> add the mean back and use the hydrostatic
    %%% approximation
    
    % PTs
    for i=1 : length(stPOD)
        
        % calculate the sampling period (s) of the PT
        nTpt = mean(seconds(diff(stPOD(i).PSdt)));
        
        % convert seapressure to water surface elevation (m NAVD88)
        stPOD(i).wse = (stPOD(i).seapress * nDbar2Pa / (nRho * nG)) + ...
                        nPTele(i);
        
        % calculate water depth above well cap (for TCM analysis)
        stPOD(i).depthTCM = stPOD(i).wse - nPODele(i);
        
        % calculate water depth above PT (for plotting time series, hydrodstatic approximation)
        stPOD(i).depthPT = stPOD(i).wse - nPTele(i);
        
        % calculate still-water as the 5-minute running mean of WSE
        if nTpt == 16  % 16 Hz RBRs
            %stPOD(i).swe = medfilt1(stPOD(i).wse, 4800);
            stPOD(i).swe = movmean(stPOD(i).wse, 4800);
        else           % 2 Hz RBRs
            stPOD(i).swe = movmean(stPOD(i).wse, 600);
        end
        
        % 5-min moving mean of depth above well cap (to track submergence 
        % of TCM by surge)
        stPOD(i).depthTCMsmooth = stPOD(i).swe - nPODele(i);
        
    end
    
    % ADV     
    % convert seapressure to water surface elevation (m NAVD88)
    stVECT.wse = (stVECT.seapress * nDbar2Pa / (nRho * nG)) + nVECTPTele;
        
    % calculate still-water as the 5-minute running mean of WSE (ADV
    % measured at 16 Hz)
    stVECT.swe = movmean(stVECT.wse, 4800);
    
    % calculate total water depth above bottom
    stVECT.depth = stVECT.wse - nVECTele;
    
    end
    % ------------------------------------------------------------------- %
%% indexTS   
    function [iMinPS1, iMaxPS1, iMinTCM1, iMaxTCM1, iMinVECT, iMaxVECT] = ...
                        indexTS(stPOD, stVECT, sSTime1, sETime1) 
    %%% find the indices for points when the pressure sensors are submerged
    %iStart = find(stPOD(1).depthPT > 0.05, 1);
    %iEnd = find(stPOD(1).depthPT > 0.05, 1, 'last');
    
    % preallocate cell arrays (I use cells here because it's helpful 
    % to distinguish between a POD index (i) and a time index {i,j})
    [iMinPS1, iMaxPS1, iMinTCM1, iMaxTCM1] = deal(cell(size(stPOD)));
    
    % find indices for submerged time series
    for i=1 : length(stPOD)
        % pressure transducers
        iMinPS1{i} = find(stPOD(i).PSdt == sSTime1{i});
        iMaxPS1{i} = find(stPOD(i).PSdt == sETime1{i});
        
        % TCMs
        iMinTCM1{i} = find(stPOD(i).TCMdt == sSTime1{i});
        iMaxTCM1{i} = find(stPOD(i).TCMdt == sETime1{i});
    end
    
    % find indices for ADV
    iMinVECT = find(stVECT.dt == sSTime1{1});
    iMaxVECT = find(stVECT.dt == sETime1{1});
    
    end
    % ------------------------------------------------------------------- %
%% waveStats    
    function [stWaveStat, stADVWaveStat] = waveStats(stPOD, stVECT, ...
                                               iMinPS1, iMaxPS1, iSeapress)
    %%% generate wave spectra and statistics

    % RRU-1, RRU-2, RRU-6 -------------------------------------------------
    nZbeg    = [0 0 0];         % burial depth at the begining
    nZend    = [0.25 0 0];      % " " at the end of the time series (m)
    bPlotWv  = false;
    n2nWaves = [15, 16; 12, 13]; % 34, 68 min
    nFmaxSS  = [0.25 0.5 0.5];  % max f that the pressure sensor can detect
    nFminIF  = 0.003;           % min IG frequency
    nFminSS  = 0.04; 	        % min f for sea-swell 
    nHpt     = 0;               % original PT height above bed 
    
    % PSD: sine multitaper method
    % (first landfall only)
    for iPOD = 1 %: 2 
        for n = 3 %: length(n2nWaves)  
            sDirName  = sprintf('WS-RRU%d-LF1-%d', nRRUid(iPOD), n2nWaves(iPOD,n));
            sFileName = sprintf('RRU%d', nRRUid(iPOD));
            nWin      = 2^n2nWaves(iPOD,n); 
            stWaveStat{iPOD,n} = ...
                    psdWaveStat_Harvey(stPOD(iPOD).seapress(iMinPS1{iPOD}:iMaxPS1{iPOD}), ...
                                     stPOD(iPOD).PSdt(iMinPS1{iPOD}:iMaxPS1{iPOD}), ...
                                     nFminIF, nFminSS, nFmaxSS(iPOD), n2nWaves(iPOD,n), ...
                                     bPlotWv(iPOD), sDirName, nHpt, nWin, sFileName, ...
                                     nZbeg(iPOD), nZend(iPOD)); %#ok<*AGROW>
        end
    end

    % MTM: slepian tapers (using inputs from package multitaper in R)
    % (NW=3, K=4)
     for i = 1 %: 2 
                
            % updated on 5/10/19 (single sided spectra) for n, hpf, 68 min 
            % window and again on 9/5/19 for new hpf (subtraction)
            cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs
            load('OUT-mtmRRU1_HPF_v2_34minNW5K9')
            
            % load original pressure data in pascals (used to calculate h)
            cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/INPUTS           
            rPb = load('RRU1_34min_psin.txt');

            % just get the old dt (didn't change)
            cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVE_STATS/
            load('mtmWaveStats_34min_0pt003.mat');
            dtPb = stMTMWaveStat34min{1,1}.dtPb;
                
            % input variables
            sDirName  = '';
            sFileName = 'RRU';
            rSnn      = snn;
            rF        = f;
            rCIlow    = CIlow;
            rCIup     = CIup;
            rFtest    = fTest;
                
            stMTMWaveStatRRU1{1,i} = mtmWaveStat_Harvey(snn, f, CIlow, ...
                           CIup, fTest, dtPb, rPb, nFminIF, nFminSS, ...
                           nFmaxSS(i), bPlotWv, sDirName, nHpt, sFileName, ...
                           nZbeg(i), nZend(i));     
                       
            %stMTMWaveStat34min{2,2} = stMTMWaveStatRRU2{1,2};
     end
     
  % ADV -------------------------------------------------------------------
    for i = 3%4 : length(n2nWaves)
        % input variables
        nZbeg    = 0;     % burial depth at the beginning 
        nZend    = 0;     % " " at the end of the time series (m)
        bPlotWv  = false;
        n2nWaves = [12, 13, 14, 15, 16]; % 4, 8, 17, 34, 68 min
        nFmaxSS  = 0.6;   % maximum frequency (visual) before rise in SnnKp
        nFminIF  = 0.003; % minimum infragravity frequency
        nFminSS  = 0.04;  % minimum frequency for the sea-swell spectra
        sDirName  = sprintf('WS-ADV-N%d', n2nWaves(i));
        sFileName = 'ADV';
        nWin      = 2^n2nWaves(i); 
        nHpt      = 0.7;  % original PT height above bed 
        % NOTE: seapress is NaNs until row 345568 
        %iSeapress = 345568;
        
        % PSD: sine multitaper method
        stADVWaveStat{1,i} = szWaveStat_Harvey(stVECT.seapress(iSeapress:end), ...
                                stVECT.dt(iSeapress:end), ...
                                nFminIF, nFminSS, nFmaxSS, n2nWaves(i), ...
                                bPlotWv, sDirName, nHpt, nWin, sFileName, ...
                                nZbeg, nZend);
    end
    
    % MTM: slepian tapers (using inputs from package multitaper in R)
    % (NW=3, K=4) --- OLD CODE --
     for iPOD = 3 
        for i = 2 %: 2 %length(stMTMps)
                
                % updated on 9/28/18 (double sided spectra) for n+
                cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVE_STATS/ADV/Bispectra/ADVsp-8min
                load('IN-mtmBisp_ADVnPLUS_8minZP.mat', 'snn', 'f', ...
                     'CIlow', 'CIup', 'fTest')
                
                % just get the old pb and dt (didn't change)
                cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/VELOCITY
                load('VECTwin_8min.mat')
                
                % input variables
                nZbeg    = 0;     % burial depth at the beginning 
                nZend    = 0;     % " " at the end of the time series (m)
                bPlotWv  = false;
                nFmaxSS  = 0.25;  % maximum frequency (visual) before rise in SnnKp
                nFminIF  = 0.003; % minimum infragravity frequency
                nFminSS  = 0.04;  % minimum frequency for the sea-swell spectra
                nHpt      = 0.7;  % original PT height above bed 
                sDirName  = '';
                sFileName = 'ADV';
                rSnn      = snn * 2;
                rF        = f;
                rCIlow    = CIlow * 2;
                rCIup     = CIup * 2;
                rFtest    = fTest;
                dtPb      = stVECTwin.dt;
                rPb       = stVECTwin.seapress;
                
                stMTMWaveStat0pt003{i,iPOD} = ...
                    mtmWaveStat_Harvey(rSnn, rF, rCIlow, rCIup, rFtest, ...
                                     dtPb, rPb, nFminIF, nFminSS, ...
                                     nFmaxSS, bPlotWv, ...
                                     sDirName, nHpt, sFileName, ...
                                     nZbeg, nZend);                                  
        end
     end
    
    % updated on 5/10/19 (single sided spectra) for n+, hpf, 68 min window
    cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV
    %load('OUT-mtmADVsp_nPlusHPF_68minNW8K15')
    load('OUT-mtmADVsp_nPlusHPF_v2_68minNW8K15')
    
    % just get the old pb (for calculation of h) and dt (didn't change)
    cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/VELOCITY
    load('VECTwin_68min.mat')

    % input variables
    nZbeg    = 0;     % burial depth at the beginning 
    nZend    = 0;     % " " at the end of the time series (m)
    bPlotWv  = false;
    nFmaxSS  = 0.50;  % maximum frequency (visual) before rise in SnnKp
    nFminIF  = 0.003; % minimum infragravity frequency
    nFminSS  = 0.04;  % minimum frequency for the sea-swell spectra
    nHpt      = 0.7;  % original PT height above bed 
    sDirName  = '';
    sFileName = 'ADV';
    dtPb      = stVECTwin.dt;
    rPb       = stVECTwin.seapress * 10000; % make Pa

    stMTMWaveStat68min{3,3} = ...
        mtmWaveStat_Harvey(snn, f, CIlow, CIup, fTest, dtPb, rPb, ...
                           nFminIF, nFminSS, nFmaxSS, bPlotWv, ...
                           sDirName, nHpt, sFileName, nZbeg, nZend);    
                     
    cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVE_STATS
    save('mtmWaveStats_68min_v2_0pt003.mat', 'stMTMWaveStat68min')                
     
    end 
    % ------------------------------------------------------------------- %
%% tiltTCM    
    function [stPOD] = tiltTCM(stPOD, nZpt)
    %%% calculate tilt from raw accelerometer components and then decompose
    %%% into tilt-easting and tilt-northing using a 3d pointer derived from
    %%% the compass heading (yaw)
    
    for iPOD=1 : length(stPOD)
         rAcc = [stPOD(iPOD).Ax, stPOD(iPOD).Ay, stPOD(iPOD).Az];
         rMag = [stPOD(iPOD).Mx, stPOD(iPOD).My, stPOD(iPOD).Mz];
         [~, ~, stPOD(iPOD).tiltN, stPOD(iPOD).tiltE, ...
                stPOD(iPOD).velocityNcalc, stPOD(iPOD).velocityEcalc] = ...
                                              tcmVectQts(rAcc, rMag, nZpt);
    end
    end
    % ------------------------------------------------------------------- %   
%% despikeVel
    function [stVECT] = despikeVel(stVECT, iSeapress, iMinVECT) %#ok<*INUSD>
    %%% identify and replace outliers in ADV time series via a
    %%% correlation threshold and SNR (QC criteria); reject time series
    %%% based on deviations from linear theory (using coherence plots)
    %%% NOTE: no windows are rejected from Harvey (all look good); rotate
    %%% velocity matrix to cross-shore and alongshore
    
    cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM
    load('VECT.mat') 

    % inputs (user specified)
    nFsamp    = 16;
    sDirName  = 'QCadv';
    bPlotADV  = false;
    n2nADV    = 17; %14; %16; %13   
    nWin      = 2^n2nADV;
    nSNR      = 15;        % SNR threshold from Nortek 
    bWinRRUs  = false;     % true: window ADV ts to match start of RRU-1 ts
    
    if bWinRRUs
        % save windows using this method as 'VECTwinSub'
        iStart  = iMinVECT;    % iMinVECT (corresp. with RRU-1) 
    else
        iStart  = iSeapress;   % full ADV time series
    end
    
    % create windows for QC'ing
    nTotal = length(stVECT.seapress(iStart:end));  % total # of indices 
    nBurst = floor(nTotal/nWin);     % number of bursts
    cntWin = 0;                      % initialize counter 
        
    for iBurst = 1 : nBurst
        
        sFileName = sprintf('Burst%d', iBurst);

        % inputs
        rVel      = [stVECT.velocityE( (iStart + cntWin) : (iStart + cntWin + nWin - 1) ) ...
                     stVECT.velocityN((iStart + cntWin): (iStart + cntWin + nWin - 1)) ...
                     stVECT.velocityUp((iStart + cntWin): (iStart + cntWin + nWin - 1))];
        rCorr     = [stVECT.corrEast((iStart + cntWin): (iStart + cntWin + nWin - 1)) ...
                     stVECT.corrNorth((iStart + cntWin): (iStart + cntWin + nWin - 1)) ... 
                     stVECT.corrUp((iStart + cntWin): (iStart + cntWin + nWin - 1))];
        rSNR      = [stVECT.S2Neast((iStart + cntWin): (iStart + cntWin + nWin - 1)) ... 
                     stVECT.S2Nnorth((iStart + cntWin): (iStart + cntWin + nWin - 1)) ...
                     stVECT.S2Nup((iStart + cntWin): (iStart + cntWin + nWin - 1))];
        rP        = stVECT.seapress((iStart + cntWin): (iStart + cntWin + nWin - 1));

        % run qcADV_Harvey
        [rVelQC] = qcADV_Harvey(rVel, rCorr, rSNR, nFsamp, [nSNR nSNR nSNR], ...
                                bPlotADV, sDirName, sFileName);

        % rotate ENU to cross-shore and alongshore: NOTE (7/1/19): as this
        % code is currently written, I don't think making the alongshore positive
        % downcoast is the typical orientation for u,v -- I think most codes
        % assume u is postive onshore, v positive to the right, as in 
        % Northing and Easting -- the only codes this will affect is something 
        % related to obtaining the correct compass direction (mean wave 
        % direction, spreading, etc.); keep this in mind when using other 
        % peoples codes (may have to change the sign of v)
                               % northing     easting      azimuth from N
        [rRad, rTan] = rt_rotate(rVelQC(:,2), rVelQC(:,1), 50.93); 
        stVECTwin.velU(:,iBurst) = -rTan; % u, cross-shore (positive onshore)
        stVECTwin.velV(:,iBurst) = -rRad; % v, alongshore (positive downcoast)        
        stVECTwin.velW(:,iBurst) = rVelQC(:,3); % positive up
        stVECTwin.dt(:,iBurst) = stVECT.dt((iStart + cntWin): (iStart + cntWin + nWin - 1));
        stVECTwin.seapress(:,iBurst) = rP;
        
        % calculate theta                            v            u
        %stVECTwin.theta(:,iBurst) = cart2compass(rVelQC(:,1),rVelQC(:,2));
       
        % calculate RMS velocity in cross-shore and alongshore directions
        stVECTwin.rmsU(:,iBurst) = rms(stVECTwin.velU(:,iBurst));
        stVECTwin.rmsV(:,iBurst) = rms(stVECTwin.velV(:,iBurst));
        
        % increase window counter
        cntWin = cntWin + nWin;
    end
    
    % set time zone of datetime
    stVECTwin.dt.TimeZone = 'America/Chicago';

    save('VECTwin_137min.mat', 'stVECTwin')
    dlmwrite('IN-VECTwin_U_137min.txt', stVECTwin.velU, 'delimiter', ...
             '\t', 'precision', 12)
    dlmwrite('IN-VECTwin_V_137min.txt', stVECTwin.velV, 'delimiter', ...
             '\t', 'precision', 12)
    
    end
    % ------------------------------------------------------------------- %
%% lowfreqRRU6
    function [rSP, dt] = lowfreqRRU6(stPOD, iMinPS1) %#ok<*DEFNU>
    %%% segment RRU-6 pressure data into windows for spectral analysis in
    %%% the very low frequency band
    
    cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/
    load('RRU.mat')

    % upsample RRU6 to 16 Hz (for some reason the interpolated time 
    % series is always one data point short; interpolates up to the last point) 
    iStart = 345601; % corresponds to 24-Aug-2017 23:00:00
    tmpRRU6= stPOD(3).seapress(iStart:end);   
    nTadv  = 1/16;  % sampling period of ADV
    nTrru  = 1/2;   % sampling period of RRU6
    rRRU6  = interp1(tmpRRU6, (1 : (nTadv/nTrru) : length(tmpRRU6)))';
    tmpDT  = stPOD(3).PSdt(iStart:end);
    dtRRU6 = (tmpDT(1) : seconds(1/16) : tmpDT(end))';
    
    % create ~137 min windows
    % inputs
    n2nRRU = 17;             
    nWin   = 2^n2nRRU;
    nTotal = length(rRRU6);          % total # of indices 
    nBurst = floor(nTotal/nWin);     % number of bursts (should be >75)
    cntWin = 0;                      % initialize counter 
    rSP = NaN(nWin,nBurst);
    
    for iBurst = 1 : nBurst
        % inputs
        %rSP(:,iBurst) = stPOD(3).seapress( (iStart + cntWin) : (iStart + cntWin + nWin - 1) );                     
        %dt(:,iBurst) = stPOD(3).PSdt( (iStart + cntWin) : (iStart + cntWin + nWin - 1) ); 
        rSP(:,iBurst) = rRRU6( (1 + cntWin) : (1 + cntWin + nWin - 1) ); 
        dt(:,iBurst)  = dtRRU6( (1 + cntWin) : (1 + cntWin + nWin - 1) );
        
        % increase window counter
        cntWin = cntWin + nWin;
    end

    % export to text file for analysis in R (only export the same array
    % size as the ADV...which turned off earlier)
    dlmwrite('RRU6_137min_psin_16Hz.txt', rSP(:,1:39)*10000, 'delimiter', ...
             '\t', 'precision', 12)
    
    end
    % ------------------------------------------------------------------- %
%% ndbcWaveStats    
    function [stNDBC] = ndbcWaveStats(stNDBC, stMTMWaveStat0pt0015)
    %%% generate wave statistics for NDBC wave spectra (sea-swell only); 
    %%% compare NDBC and ADCIRC wave stats against infragravity wave 
    %%% statistics from Follets Island (zero-padded, 17 min), averaged to 
    %%% 1 hr intervals, and check the cross-correlation using the full time
    %%% series and one-hr averaged time series
    
    %%% NOTE: (2/2019) both the NDBC and ADCIRC datetimes were off when
    %%% this analysis was originally completed (ADCIRC by one hour, NDBC
    %%% by 6 hours [UTC]) - so much of this code is incorrect but keeping 
    %%% for reference)
    
    % input options
    bNDBC    = true;   % if false, process ADCIRC data
    bMATA    = false;  % if false, compare against Follets Island ADV (IG only)
    bSIG     = false;  % if true, analyze SIG band, not IG
    iNDBC    = 1;      % 1 is Freeport, 2 is Galveston
    if bSIG 
        iNDBClow = 577; % 577 = 8/25 03:00 (ADV analysis), 578 = 8/25 04:00 (RRU)
    else
        iNDBClow = 573; % 573 = 8/24 23:00 (ADV analysis), 578 = 8/25 04:00 (RRU)
    end
    iNDBCup  = 662;   %#ok<*NASGU> % 662 = 8/28 16:00 (ADV analysis), 605 = 8/26 07:00 (RRU)
    iADClow  = 41;
    iADCup   = 130; % 130 = 8/28 16:00, 80 = 8/26 14:00
    iNode    = 3083010; %3083010 = 8.11 m FI, 2542559 = 8.34 m MATA
    nG       = 9.81;
    
    if bNDBC
        cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/NDBC
        load('NDBC.mat') %#ok<*LOAD>
        
        % NDBC SEA-SWELL STATISTICS ----------------------------------------- %
        % spectral moments 
        rM0ss   = trapz(stNDBC{1,iNDBC}.f(:,1), ...
                        stNDBC{1,iNDBC}.psd); % [m^2] = variance (sigma^2)
        rM1ss   = trapz(stNDBC{1,iNDBC}.f(:,1), ...
                        stNDBC{1,iNDBC}.psd .* stNDBC{1,iNDBC}.f);
        rM2ss   = trapz(stNDBC{1,iNDBC}.f(:,1), ...
                        stNDBC{1,iNDBC}.psd .* stNDBC{1,iNDBC}.f.^2);
        rM4ss   = trapz(stNDBC{1,iNDBC}.f(:,1), ...
                        stNDBC{1,iNDBC}.psd .* stNDBC{1,iNDBC}.f.^4);

        rNrmsSS = sqrt(rM0ss);                         % std deviation (sigma)
        stNDBC{1,iNDBC}.HrmsSS = 2 * rNrmsSS;          % rms wave height
        stNDBC{1,iNDBC}.Hm0SS  = 4.004 * rNrmsSS;      % 0-moment wave height (~H1/3)
        stNDBC{1,iNDBC}.Tm01SS = rM0ss ./ rM1ss;       % mean wave period 
        stNDBC{1,iNDBC}.Tm02SS = sqrt(rM0ss ./ rM2ss); % zero-upcrs
        stNDBC{1,iNDBC}.L0 = (nG * stNDBC{1,iNDBC}.Tm01SS.^2) / (2*pi); % deep water wavelength
        
        % peak frequency (in each column)
        [~, iMaxSS] = max(stNDBC{1,iNDBC}.psd, [], 1);
        stNDBC{1,iNDBC}.FpeakSS = stNDBC{1,iNDBC}.f(iMaxSS)';
        
        % spectral width (Cartwright and Longuet-Higgins, 1956)
        stNDBC{1,iNDBC}.sw1 = sqrt( 1 - (rM2ss).^2 ./ (rM0ss./rM4ss) );
        % spectral width (Longuet-Higgins, 1975), wave group statistics
        stNDBC{1,iNDBC}.sw2 = sqrt( (rM0ss.*rM2ss ./ rM1ss.^2) - 1 );
        % peakedness (Goda, 1970) - wave grouping of wind waves
        stNDBC{1,iNDBC}.qp = 2./(rM0ss.^2) .* trapz(stNDBC{1,iNDBC}.f(:,1), ...
                        stNDBC{1,iNDBC}.psd.^2 .* stNDBC{1,iNDBC}.f);
        % spectral peakedness (Medina and Hudspeth, 1987)
        stNDBC{1,iNDBC}.qe = (2*rM1ss)./(rM0ss.^3) .* trapz(stNDBC{1,iNDBC}.f(:,1), ...
                        stNDBC{1,iNDBC}.psd.^2);
        
        % SWELL and SEA ONLY -------------------------------------------- %
        % spectral moments swell
        rM0s   = trapz(stNDBC{1,iNDBC}.f(1:19,1), ...
                       stNDBC{1,iNDBC}.psd(1:19,:)); 
        rNrmsS = sqrt(rM0s);                         
        stNDBC{1,iNDBC}.Hm0Swell = 4.004 * rNrmsS;      
        % spectral moments sea
        rM0sea = trapz(stNDBC{1,iNDBC}.f(20:end,1), ...
                       stNDBC{1,iNDBC}.psd(20:end,:)); 
        rNrmsSea = sqrt(rM0sea);     
        stNDBC{1,iNDBC}.Hm0Sea  = 4.004 * rNrmsSea;
        
        % subset of NDBC used for linear regressions
        rHm0Swell = stNDBC{1,iNDBC}.Hm0Swell(iNDBClow:iNDBCup);
        rHm0SS  = stNDBC{1,iNDBC}.Hm0SS(iNDBClow:iNDBCup);
        rHm0Sea = stNDBC{1,iNDBC}.Hm0Sea(iNDBClow:iNDBCup);
        rTpSS   = 1 ./ stNDBC{1,iNDBC}.FpeakSS(iNDBClow:iNDBCup);
        rTm01SS = stNDBC{1,iNDBC}.Tm01SS(iNDBClow:iNDBCup);
        rTm02SS = stNDBC{1,iNDBC}.Tm02SS(iNDBClow:iNDBCup);
        dtSS    = stNDBC{1,iNDBC}.dt(iNDBClow:iNDBCup);
        rL0     = stNDBC{1,iNDBC}.L0(iNDBClow:iNDBCup);
        rSW1    = stNDBC{1,iNDBC}.sw1(iNDBClow:iNDBCup);
        rSW2    = stNDBC{1,iNDBC}.sw2(iNDBClow:iNDBCup);
        rQp     = stNDBC{1,iNDBC}.qp(iNDBClow:iNDBCup);
        rQe     = stNDBC{1,iNDBC}.qe(iNDBClow:iNDBCup);
    
    else
        cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ADCIRC/Hindcast43-NHC_best_track
        load('ADCIRC.mat')
        
        % isolate ADCIRC data for Hs (sea-swell only), Tm0, Tp, dt
        rHm0SS   = rHsADCIRC(iNode, iADClow:iADCup);
        rTm01SS  = rTm0ADCIRC(iNode, iADClow:iADCup);
        rTpSS    = rTpsADCIRC(iNode, iADClow:iADCup);
        % radiation stress gradient (wave induced force per unit surface
        % area) - x and y are in the output coordinate system (x positive
        % upcoast, y positive cross-shore)
        rRadX    = rRadXADCIRC(iNode, iADClow:iADCup);
        rRadY    = rRadYADCIRC(iNode, iADClow:iADCup);
        dtSS     = rTimeADCIRC(iADClow:iADCup);
        rL0      = (nG * rTm01SS.^2) / (2*pi); % deep water wavelength

    end
    
    % average infragravity statistics from ~17 min to ~1 hr intervals using
    % a moving mean 
    n = 4; % use a 68 minute window
    
    % OLD method (full IG band): now only use for dt vector
    cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS
    load('mtmWaveStats_17min_0pt0015_nw3k4.mat')
    
    % new method: subdivide IG band
    cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/ADV/PowerSpectra/MTM_R
    load('mtmADVsp_17min_ZP.mat', 'snn', 'f')
    iLFlow   = [7 26 57];   % low cutoffs (0.0029, 0.0122, 0.0273) 
    iLFhigh  = [25 56 83];  % high cutoffs (0.0117, 0.0269, 0.04) 

    % individual bands
    for i = 1 : length(iLFlow)
        rM0snnLF = trapz(f(iLFlow(i):iLFhigh(i),1), snn(iLFlow(i):iLFhigh(i),:));
        rHm0LF(i,:) = 4.004 * sqrt(rM0snnLF); % each row corresponds to a band
    end

    % Matagorda RRU data
    if bMATA
        rHm0Infa = movmean(stMTMWaveStat0pt0015{1,1}.Hm0Inf, n);
        
        % linearly interpolate the moving average and subsample at each hr
        rHm0Infa = interp1(datenum(stMTMWaveStat0pt0015{1,1}.dtPb(1,:)), ...
                           rHm0Infa, datenum(dtSS));

    % Follets Island ADV
    else
        if ~bSIG
            rHm0Infa = movmean(rHm0LF(2,:), n);
            %rHm0Infa = movmean(stMTMWaveStat0pt0015{1,3}.Hm0Inf, n);
            rHm0Infa = interp1(datenum(stMTMWaveStat0pt0015{1,3}.dtPb(1,:)), ...
                              rHm0Infa, datenum(dtSS));

%         % lets try the 68 minute windows
%         cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/SensitivityTests/LowFreqSignals/
%         load('mtmADVsp_68min.mat', 'snn', 'f', 'dt')
%         iLow   = 26;   % corresponds to 0.003 Hz
%         iHigh  = 329;   % corresponds to 0.04 Hz
%         rM0snn = trapz(f(iLow:iHigh,1), snn(iLow:iHigh,:));
%         rHm0   = 4.004 * sqrt(rM0snn); % 4*std
%         rHm0Infa = interp1(datenum(dt(1,:)), rHm0, datenum(dtSS));
        
        else
            % calculate SIG Hm0 using the zero-padded spectra 
            cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/SensitivityTests/LowFreqSignals/
            load('mtmADVspZP_68minSUB.mat', 'snn', 'dt', 'f')
            iLFlow   = 4;   % corresponds to low SIG cutoff  (4 = 0.00036) 
            iLFhigh  = 25;  % corresponds to high SIG cutoff (25 = 0.0029) 
            rM0snnLF = trapz(f(iLFlow:iLFhigh,1), snn(iLFlow:iLFhigh,:));
            rHm0LF   = 4.004 * sqrt(rM0snnLF); % 4*std
            rHm0LF   = interp1(datenum(dt(1,:)), rHm0LF, datenum(dtSS));
        end
        
        % tides
        cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM
        load('TIDE.mat') % SLP predicted tides (m)
        rTide    = interp1(datenum(stTide.dt), stTide.pred, datenum(dtSS));

        % surge (SWE - tidal residual)
        nVECTele = -0.746; % ADV elevation
        rSWE   = interp1(datenum(stMTMWaveStat0pt0015{1,3}.dtPb(1,:)), ...
                          movmean(stMTMWaveStat0pt0015{1,3}.h+nVECTele, 5), ...
                          datenum(dtSS)); 
        rSurge = rSWE - rTide;
        
        % surge from ADCIRC
        %rSurge = rWseADCIRC(iNode,iADClow:iADCup)-rTide';
      % ----------------------------------------------------------------- %
%     
%     % OLD METHOD - linear regression ------------------------------------
%     % use polyfit to compute a linear regression that predicts y from x,
%     % then polyval to use the regression to predict y; then compute the
%     % residuals
%     % case 1: Hinf (y) vs. H0 (x)
%     pH0HinfMS = polyfit(rH0NDBC, rHm0Infmsa, 1);
%     pH0HinfPS = polyfit(rH0NDBC, rHm0Infpsa, 1);
%     pFitMS1 = polyval(pH0HinfMS, rH0NDBC);
%     pFitPS1 = polyval(pH0HinfPS, rH0NDBC);
%     % R2 ca
%     yResidMS1 = rHm0Infmsa - pFitMS1;
%     ssResidMS1 = sum(yResidMS1.^2);
%     ssTotalMS1 = (length(rHm0Infmsa)-1) * var(rHm0Infmsa);
%     r2 = 1 - ssResidMS1/ssTotalMS1;
%     r2adj = 1 - ssResidMS1/ssTotalMS1 * (length(rHm0Infmsa)-1)/ ...
%                         (length(rHm0Infmsa)-length(pH0HinfMS));
%     % figure
%     figure; scatter(rHm0SS, rHm0InfaADV); hold on; plot(rH0NDBC, pFitMS1);
%             hold on; 
%             scatter(rH0NDBC, rHm0Infpsa); hold on; plot(rH0NDBC, pFitPS1);
%             legend('multisine tapers', 'MS - linear regression', ...
%                    'prolate spheroidal tapers', 'PS - linear regression')
%
%     % OLD METHOD - TRY USING THE WHOLE TIME SERIES FOR CROSS-CORRELATION          
%     % 1) resample the signal with the lower sampling rate (NDBC)
%      rT1 = interp1(datenum(datenum(dtSS)), rHm0SS, ...
%            datenum(stMTMWaveStat0pt0015{1,3}.dtPb(1,:)));
%      rT2 = stMTMWaveStat0pt0015{1,3}.Hm0Inf; 
% 
%     % cross-correlate (ADV lags offshore) and normalize using the
%     % 'coeff' option; NOTE: switching the order of the operands
%     % only reverses the sequence
%     [rXcorr, rLag] = xcorr(rT1(1:313), rT2(1:313), 'coeff');
%     figure; plot(rLag, rXcorr)
%     xlabel('17 min windows')
%     ylabel('Xcorr (ADV lags offshore)')
%     ylim([-1 1])
    if ~bSIG
        % 2) cross-correlate (ADV lags offshore) using the one-hr averaged ADV
        % time series and normalize using the 'coeff' option; NOTE: switching 
        % the order of the operands only reverses the sequence
        [rXcorr, rLag] = xcorr(rTm01SS, rHm0Infa, 'coeff');
        figure; plot(rLag/(1/60*60), rXcorr)
        xlabel('time (hr)')
        ylabel('Xcorr (ADV lags offshore)')
        ylim([0 1])
    end
    
    end
    
    if ~bSIG
        if bNDBC
            % ADCIRC doesn't have swell only, Tm02SS, or spectral width option 
            lm3 = fitlm(rHm0Swell, rHm0Infa);
                figure; plot(lm3); plotFancyAxis 
                xlabel('H0_{swell} (m)', 'Interpreter', 'tex', 'FontSize', 15) 
                ylabel('Hm0_{inf} (m)', 'Interpreter', 'tex', 'FontSize', 15) 
            fitlm(rTm02SS, rHm0Infa)
            fitlm(rSW1, rHm0Infa)
            fitlm(rSW2, rHm0Infa)
            fitlm(rQp, rHm0Infa)
            fitlm(rQe, rHm0Infa)
            fitlm(rHm0Sea, rHm0Infa)
        else
            fitlm(rRadX, rHm0Infa)
            fitlm(rRadY, rHm0Infa)   
        end
            lm4 = fitlm(rHm0SS.^2 .* rTpSS, rHm0Infa); 
                figure; plot(lm4); plotFancyAxis
                xlabel('H_0^2*T_p (m)', 'Interpreter', 'tex', 'FontSize', 15) 
                ylabel('Hm0_{inf} (m)', 'Interpreter', 'tex', 'FontSize', 15) 
            fitlm(rHm0SS, rHm0Infa)
            lm2 = fitlm(sqrt(rHm0SS .* rL0), rHm0Infa);
                figure; plot(lm2); plotFancyAxis
                xlabel('(H0*L0)^{1/2} (m)', 'Interpreter', 'tex', 'FontSize', 15) 
                ylabel('Hm0_{inf} (m)', 'Interpreter', 'tex', 'FontSize', 15) 
            fitlm(rTpSS, rHm0Infa)
            lm1 = fitlm(rTm01SS, rHm0Infa);
                figure; plot(lm1); plotFancyAxis
                xlabel('T_0 (sec)', 'Interpreter', 'tex', 'FontSize', 15)  
                ylabel('Hm0_{inf} (m)', 'Interpreter', 'tex', 'FontSize', 15) 
            % non-linear model
%             modelfun = @(b,x) b(1)*x(:,1).^b(2);
%             beta0 = [0.1 2];
%             nlm = fitnlm(rTm01SS, rHm0Infa, modelfun, beta0)
%             Xnew = [0 : 0.5 : 10];
%             [ynew,ynewci] = predict(nlm,Xnew');
            %fitlm([rHm0SS' rTm01SS'], rHm0Infa, 'quadratic') 
            fitlm(rSurge, rHm0Infa)
            fitlm(rSWE, rHm0Infa)
    else
        % SIG
        fitlm(rTide, rHm0LF)
        fitlm(rHm0Swell, rHm0LF)
        fitlm(rHm0Sea, rHm0LF)
        fitlm(rHm0SS.^2 .* rTpSS, rHm0LF) 
        fitlm(rHm0SS, rHm0LF)
        fitlm(sqrt(rHm0SS .* rL0), rHm0LF)
        fitlm(rTpSS, rHm0LF)
        fitlm(rTm01SS, rHm0LF)
        fitlm(rSW1, rHm0LF)
        fitlm(rSW2, rHm0LF)
        fitlm(rQp, rHm0LF)
        fitlm(rQe, rHm0LF)
        %fitlm(rRadX, rHm0LF)
        %fitlm(rRadY, rHm0LF)
        fitlm([rHm0SS' rTm01SS'], rHm0LF, 'quadratic')
    end
          
    end
    % ------------------------------------------------------------------- %
%% create OBS calibration curves 
%%% the following is the calibration of OBS1 and OBS2 using a 5 gallon
%%% bucket and the post-storm sediment samples complete in June 2018; in
%%% April 2019, a new bucket test was completed to refine the curve for low
%%% sediment concentrations (PART 2); note that Will read the data in for
%%% me (and calculated the median) and I just saved the variables to stCAL
%%% after reviewing his code

    function [stVECT] = obsCalCurve(stVECT, iSeapress, stVECTwin)
    
    cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM
    load('CAL.mat')
    load('VECT.mat')
    
    %%% FIRST OBS CALIBRATION --------------------------------------------
    
    % cumulative concentration of sediment added to 3 gallon bucket
    rConcPost = [0 28.16 54.54 81.94 111.84 144.26 169.0 183.44 199.66 218.78 247.42]; % g/3-gal
    rConcPost = rConcPost/3/3.785; % make g/gal--> g/L = kg/m^3
    
    %rAmbOBS1 = 597; % red OBS (located 63-70 cm above the bed during Harvey)
    rAmbOBS2 = 996; % orange OBS (located 35-42 cm above the bed during Harvey)

    % start and end time for calculating the mean 
    sSTime = {'21-Jun-2018 17:14:57', '21-Jun-2018 17:15:22', ...
              '21-Jun-2018 17:15:50', '21-Jun-2018 17:16:22', ...
              '21-Jun-2018 17:16:56', '21-Jun-2018 17:17:28', ...
              '21-Jun-2018 17:18:03', '21-Jun-2018 17:18:37', ...
              '21-Jun-2018 17:19:15', '21-Jun-2018 17:19:53'};
    sETime = {'21-Jun-2018 17:15:17', '21-Jun-2018 17:15:42', ...
              '21-Jun-2018 17:16:10', '21-Jun-2018 17:16:42', ...
              '21-Jun-2018 17:17:16', '21-Jun-2018 17:17:48', ...
              '21-Jun-2018 17:18:23', '21-Jun-2018 17:18:57', ...
              '21-Jun-2018 17:19:35', '21-Jun-2018 17:20:13'};
        
    % find indices for start and end time
    for i=1 : length(sSTime)
    
        % pressure transducers
        iMin = find(stCAL.dtPOST == sSTime{i});
        iMax = find(stCAL.dtPOST == sETime{i});
        rMean(i) = mean(stCAL.OBS2.post(iMin:iMax));
    end

    % % fit quadratic model for OBS2 data
    % rQuad = fitlm(rMean, rConcPost, 'quadratic');
    % stVECT.OBS2c = predict(rQuad, stVECT.OBS2);

    % predict using a nonlinear function to the quadratic goes through zero
    polyfctn = @(b,x) (b(1).*x + b(2).*x.^2);
    rQuad2 = fitnlm(rMean, rConcPost, polyfctn, [0,0]);
    stVECT.OBS2c = predict(rQuad2, stVECT.OBS2);

    % plots to check (plot doesn't work here for the nonlinear model)
    %figure; scatter(rMean, rConcPost); hold on; plot(rQuad1)
    figure; plot(stVECT.dt, stVECT.OBS2c)
    
    %save('VECT.mat', 'stVECT')
    
    %%% SECOND OBS CALIBRATION (THIS IS WHAT I USED FOR MY THESIS) --------

    % fit quadratic model for both the shallow (orange) and deep (red) sensors
    rQuadShall = fitlm(stCAL.BT2.shall, stCAL.BT2.load, 'quadratic');
    figure; scatter(stCAL.BT2.shall, stCAL.BT2.load); hold on; plot(rQuadShall)
    stVECT.BT2.OBS2c = predict(rQuadShall, stVECT.OBS2);
    
    rQuadDeep = fitlm(stCAL.BT2.deep, stCAL.BT2.load, 'quadratic');
    figure; scatter(stCAL.BT2.deep, stCAL.BT2.load); hold on; plot(rQuadDeep)
    stVECT.BT2.OBS1c = predict(rQuadDeep, stVECT.OBS1);
    
    % set anything that is negative to zero
    stVECT.BT2.OBS1c(stVECT.BT2.OBS1c < 0) = 0;
    stVECT.BT2.OBS2c(stVECT.BT2.OBS2c < 0) = 0;
    
    % plot
    figure; plot(stVECT.dt, stVECT.BT2.OBS1c, stVECT.dt, stVECT.BT2.OBS2c)
    plotFancyAxis
    ylabel('concentration (g/L)')
    
    % calculate the suspended load
    stVECT.C = stVECT.BT2.OBS2c + stVECT.BT2.OBS1c*0.28; %(kg/m^2)
    
    %save('VECT.mat', 'stVECT')
    %save('CAL.mat', 'stCAL')
    
    %%% WINDOW ------------------------------------------------------------
    
    % window to 17 minute intervals and remove background value check for background value, remove, 
    iStart = iSeapress;   % full ADV time series
    nWin   = 2^13;        % sampling frequency is 16 Hz
    nTotal = length(stVECT.OBS2c(iStart:end));  % total # of indices 
    nBurst = floor(nTotal/nWin);     % number of bursts
    cntWin = 0;                      % initialize counter 
    
    for iBurst = 1 : nBurst

        % suspended sed concentration
        %rSusSed(:,iBurst) = stVECT.OBS2c( (iStart + cntWin) : (iStart + cntWin + nWin - 1) );
        rSusSed(:,iBurst) = stVECT.C( (iStart + cntWin) : (iStart + cntWin + nWin - 1) );
        
        % increase window counter
        cntWin = cntWin + nWin;
    end  
    
    % save to the window'ed structure
    cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/VELOCITY/
    load('VECTwin_8min.mat')
    stVECTwin.C = rSusSed; 
    stVECTwin.units = 'seapress = dbar, vel = m/s, dt = LST, OBS = g/L, C = kg/m^2';
    save('VECTwin_8min.mat', 'stVECTwin')
    
    % export to text file for analysis in R
    dlmwrite('IN-VECTwin_SL_8min.txt', rSusSed, 'delimiter', ...
             '\t', 'precision', 12)
    
end
