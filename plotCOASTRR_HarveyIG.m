% -------------------------plotCOASTRR_HarveyIG---------------------------%
% Purpose: This read script creates plots related to IG wave dynamics (and 
% more generally water level, sensitivity testing, etc.) using post-
% processed data collected during Hurricane Harvey. All publication related
% figures are created individually and located within the manuscript
% folder.
%
% SEE ALSO: processCOASTRR_Harvey.m, plotCOASTRR_HarveySIG.m, 
%           plotCOASTRR_HarveyMorpho.m
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       09/15/17        KAnarde             Original code
%       11/18/18        KAnarde             Removed old code, IG wave only
%
%------------------------------preamble-----------------------------------%

% specify custom calibration curve (X = pitch[deg], Y = speed [cm/s])
rTCMcalX = [0; 1; 2; 3; 5; 12.5; 20; 27; 34; 40; 46; 52; 57; 68; 75; 90];
rTCMcalY = [0; 0.5; 5; 10; 20; 30; 40; 50; 60; 70; 80; 90; 100; 125; 150; 200];
rTCMcalY = rTCMcalY/100; % convert cm/s to m/s

% conversion parameters
nKpa2Dbar = 0.1;
nDbar2Kpa = 10;
nDbar2Pa  = 10000;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plots~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%% checkPitchRoll
% check the heading, pitch, and roll of the Vector throughout the
% deployment to see when/if the Vector became unhinged

figure
% heading
subplot(3,1,1)
plot(stVECT.dtSEN, stVECT.heading)
title('Vector Heading, Roll, Pitch')
ylabel('heading (deg)')

% roll
subplot(3,1,2)
plot(stVECT.dtSEN, stVECT.roll)
ylabel('roll (deg)')

% pitch
subplot(3,1,3)
plot(stVECT.dtSEN, stVECT.pitch)
ylabel('pitch (deg)')

%% checkVECTvel
% check the xyz velocity of the Vector throughout the deployment to see 
% dwhen/if the Vector became unhinged

figure
% east
subplot(4,1,1)
plot(stVECT.dt, stVECT.velocityE)
title('Vector Velocity')
ylabel('East (m/s)')

% north
subplot(4,1,2)
plot(stVECT.dt, stVECT.velocityN)
ylabel('North (m/s)')

% up
subplot(4,1,3)
plot(stVECT.dt, stVECT.velocityUp)
ylabel('Up (m/s)')

% pressure
subplot(4,1,4)
plot(stVECT.dt, stVECT.pressure)
ylabel('pressure (dbar)')

%% checkPressure
% quick pressure plots from transducers

figure
plot(stPOD(1).PSdt, stPOD(1).pressure, ...
     stPOD(2).PSdt, stPOD(2).pressure, ...
     stPOD(3).PSdt, stPOD(3).pressure, ...
     stPOD(4).PSdt, stPOD(4).pressure)
ylabel('pressure (dbar)')
title('Raw pressure (w/RRU-1 correction)')
legend('RRU-1', 'RRU-2', 'RRU-6', 'RRU-7')

figure
plot(stPOD(1).PSdt, stPOD(1).seapress, ...
     stPOD(2).PSdt, stPOD(2).seapress, ...
     stPOD(3).PSdt, stPOD(3).seapress, ...
     stPOD(4).PSdt, stPOD(4).seapress)
ylabel('pressure (dbar)')
title('Raw sea-pressure (w/RRU-1 correction)')
legend('RRU-1', 'RRU-2', 'RRU-6', 'RRU-7')

figure
plot(stPOD(1).PSdt, stPOD(1).depthTCMsmooth, ...
     stPOD(2).PSdt, stPOD(2).depthTCMsmooth, ...
     stPOD(3).PSdt, stPOD(3).depthTCMsmooth, ...
     stPOD(4).PSdt, stPOD(4).depthTCMsmooth)
ylabel('depth (m)')
title('Depth of water column above well cap (5-min moving mean)')
legend('RRU-1', 'RRU-2', 'RRU-6', 'RRU-7')

figure
plot(stPOD(1).PSdt, stPOD(1).depthPT, ...
     stPOD(2).PSdt, stPOD(2).depthPT)
     %stPOD(3).PSdt, stPOD(3).depthPT, ...
     %stPOD(4).PSdt, stPOD(4).depthPT)
ylabel('depth (m)')
title('Depth of water column above pressure transducer')
legend('RRU-1', 'RRU-2')

figure
plot(stPOD(1).PSdt, stPOD(1).swe, ...
     stPOD(2).PSdt, stPOD(2).swe, ...
     stPOD(3).PSdt, stPOD(3).swe, ...
     stPOD(4).PSdt, stPOD(4).swe)
ylabel('stillwater elevation (m NAVD88)')
title('Stillwater elevation (5-min moving mean of WSE)')
legend('RRU-1', 'RRU-2', 'RRU-6', 'RRU-7')

figure
plot(stVECT.dt, stVECT.seapress, stVECT.dt, stVECT.pressure)
ylabel('pressure (dbar)')
title('ADV pressure vs. seapressure')
legend('sea-pressure', 'pressure')

figure
plot(stVECT.dt, stVECT.swe)
ylabel('stillwater elevation (m NAVD88)')
title('ADV stillwater elevation (5-min moving mean of WSE)')

%% landfall time series plots (wse and velocity) at Matagorda Peninsula

figure
subplot(2,1,1)
hAx(1) = gca;
plot(stPOD(1).PSdt, stPOD(1).wse)
ylabel('WSE (m NAVD88)')
title('Matagorda Peninsula')
legend('Backshore')
plotFancyAxis()

subplot(2,1,2)
hAx(2) = gca;
plot(stPOD(2).PSdt, stPOD(2).wse)
ylabel('WSE (m NAVD88)')
legend('Back-barrier')
plotFancyAxis()

linkaxes(hAx, 'x')

%% movie of stillwater elevation accross center beach profile (Hog)
% NOTE: the Backshore definitely accreted post-storm so actual water depth
% would be larger during the storm at this location; however, I think the
% channel throat stayed approximately the same elevation.

nInt = 15;
cntFigs = round(minutes(datetime(sETime1{1})-datetime(sSTime1{1}))/nInt); % # figs
hFigs = zeros(cntFigs, 1); % preallocate array for figs (by handle)
nScrsz = get(groot, 'ScreenSize');   % get screen size

cntTms = {0, 0}; % set counters for time steps (1 minute increments)

% center profile (immediately following storm) at Matagorda Peninsula
rCsorted = sortrows([stBP.HIcent(1).lat stBP.HIcent(1).ele]);
rCsorted_lat = rCsorted(:,1);
rCsorted_ele = rCsorted(:,2);

% east profile (immediately following storm) at Matagorda Peninsula
rEsorted = sortrows([stBP.HIeast(1).lat stBP.HIeast(1).ele]);
rEsorted_lat = rEsorted(:,1);
rEsorted_ele = rEsorted(:,2);

% colormap
cmap = parula(5);

for tms = 1 : cntFigs
    hFigs(tms) = figure;
    %set(hFigs(tms), 'Visible', 'off')
    subplot(2,1,1)
    % plot time series
    plot(stPOD(1).PSdt(iMinPS1{1}:iMaxPS1{1}), stPOD(1).wse(iMinPS1{1}:iMaxPS1{1}), 'Color', cmap(1,:))
    hold on
    plot(stPOD(2).PSdt(iMinPS1{2}:iMaxPS1{2}), stPOD(2).wse(iMinPS1{2}:iMaxPS1{2}), 'Color', cmap(3,:))
    hold on
    % plot points
    plot(stPOD(1).PSdt(iMinPS1{1} + cntTms{1}), stPOD(1).wse(iMinPS1{1} + cntTms{1}), 'k*')
    hold on
    plot(stPOD(2).PSdt(iMinPS1{2} + cntTms{2}), stPOD(2).wse(iMinPS1{2} + cntTms{2}), 'k*')
    % formatting
    ylabel('water surface elevation (m NAVD88)')
    title('Water surface elevation - Matagorda Peninsula')
    legend('Backshore', 'Back-barrier')
    plotFancyAxis()
    hold off

    % plot center profile
    subplot(2,1,2)
    plot(rCsorted_lat, rCsorted_ele, 'k', 'LineWidth', 1)
    hold on
    % location of RRU-2
    text(28.616407, 0.260,'RRU-2', 'HorizontalAlignment', 'right')
    hold on
    % plot east profile
    plot(rEsorted_lat, rEsorted_ele, ':k', 'LineWidth', 1)
    hold on
    % location of RRU-1
    text(28.612403, 0.9315,'RRU-1', 'HorizontalAlignment', 'right')
    legend('center profile (RRU-2)', 'east profile (RRU-1)')
    hold on
    % stillwater level (RRU-1)
    hLineRRU1 = refline(0, stPOD(1).swe(iMinPS1{1} + cntTms{1}));
    hLineRRU1.Color = cmap(1,:);
    hold on
    % stillwater level (RRU-2)
    hLineRRU2 = refline(0, stPOD(2).swe(iMinPS1{2} + cntTms{2}));
    hLineRRU2.Color = cmap(3,:);
    % formatting
    ylabel('elevation (m NAVD88)')
    lgd = legend('Center Profile', 'East Profile', 'Backshore (WSE)', 'Back-barrier (WSE)');
    ldg.Location = 'southeast';
    plotFancyAxis()
    set(gca, 'Ylim', [-2 5])
    hold off

    % make figure large
    set(hFigs(tms), 'Position', ...
        [1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);
    
    % add to counters
    cntTms{1} = cntTms{1} + (16*60*nInt); % 16 Hz * 60 sec * 15 min (RRU-1)
    cntTms{2} = cntTms{2} + (2*60*nInt); % 2 Hz * 60 sec * 15 min (RRU-2)
end

% create video object 
v = VideoWriter(sprintf('%s.mp4', 'Test'),'MPEG-4');
v.FrameRate = 15;        % fast movie (does not resize)

% set visibility on
%set(hFigs, 'Visible', 'on')

% loop through figures and make movie frames
for j=1 : length(hFigs)
    figure(hFigs(j));                % make the handle a figure
    mMov(j) = getframe(hFigs(j));    %#ok<SAGROW> % convert figure to a movie frame
end

% tidy up
open(v)
writeVideo(v, mMov)
close(v)

%close figures
close all

%% landfall time series plots (lpf fse, speed, bearing) at Matagorda Peninsula 
%%% 68 minute intervals

nScrsz = get(groot, 'ScreenSize');   % get screen size

% temp = stTCMwin{1,5}.bearing(:,14);
% for n = 1 : length(stTCMwin{1,5}.bearing(:,14))
%     if stTCMwin{1,5}.bearing(n,14) < 45
%          temp(n) = 360 + stTCMwin{1,5}.bearing(n,14);
%     end
% end

% FSE, bearing, speed
for i = 11 %: size(stWaveStat{1,5}.FSElpf, 2)
    
    hFig(i) = figure; %#ok<SAGROW>
    subplot(3,2,1)
    hAx(1) = gca;
    plot(stWaveStat{1,5}.dtPb(:,i), stWaveStat{1,5}.FSElpf(:,i))
    ylabel('FSE_{lpf} (m)')
    title(sprintf('Backshore - %s', stWaveStat{1,5}.dtPb(1,i)))
    legend(sprintf('h = %f', stWaveStat{1,5}.h(1,i)))
    plotFancyAxis()

    subplot(3,2,3)
    hAx(2) = gca;
    plot(stTCMwin{1,5}.dt(:,i), stTCMwin{1,5}.bearingLPF(:,i))
    ylabel('bearing_{lpf} (deg)')
    %ylabel('bearing (deg)')
    plotFancyAxis()

    subplot(3,2,5)
    hAx(3) = gca;
    plot(stTCMwin{1,5}.dt(:,i), stTCMwin{1,5}.speedLPF(:,i))
    ylabel('speed_{lpf} (m/s)')
    %ylabel('speed (m/s)')
    plotFancyAxis()

    subplot(3,2,2)
    hAx(4) = gca;
    plot(stWaveStat{2,5}.dtPb(:,i), stWaveStat{2,5}.FSElpf(:,i))
    ylabel('FSE_{lpf} (m)')
    title(sprintf('Back-barrier - %s', stWaveStat{2,5}.dtPb(1,i)))
    legend(sprintf('h = %f', stWaveStat{2,5}.h(1,i)))
    plotFancyAxis()

    subplot(3,2,4)
    hAx(5) = gca;
    plot(stTCMwin{2,5}.dt(:,i), stTCMwin{2,5}.bearingLPF(:,i))
    ylabel('bearing_{lpf} (deg)')
    %ylabel('bearing (deg)')
    plotFancyAxis()

    subplot(3,2,6)
    hAx(6) = gca;
    plot(stTCMwin{2,5}.dt(:,i), stTCMwin{2,5}.speedLPF(:,i))
    ylabel('speed_{lpf} (m/s)')
    %ylabel('speed (m/s)')
    plotFancyAxis()
    
    linkaxes(hAx, 'x')

    set(hFig(i), 'Position', ...
        [1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);
    
    print('-depsc2', sprintf('f%d.eps' , hFig(i).Number))
    %saveas(hFig(i), sprintf('FSE_TCM-burst%d.fig', i));
end

% water depth above PT and TCM
for i = 11%1 : size(stPSwin{1,5}.dt, 2)
    
    hFig(i) = figure;
    subplot(2,2,1)
    plot(stPSwin{1,5}.dt(:,i), stPSwin{1,5}.depthPT(:,i))
    ylabel('water depth above PT (m)')
    title(sprintf('Backshore - %s', stPSwin{1,5}.dt(1,i)))
    plotFancyAxis()
    
    subplot(2,2,3)
    plot(stPSwin{1,5}.dt(:,i), stPSwin{1,5}.depthTCM(:,i))
    ylabel('water depth above TCM (m)')
    plotFancyAxis()
    
    subplot(2,2,2)
    plot(stPSwin{2,5}.dt(:,i), stPSwin{2,5}.depthPT(:,i))
    ylabel('water depth above PT (m)')
    title(sprintf('Back-barrier - %s', stPSwin{2,5}.dt(1,i)))
    plotFancyAxis()
    
    subplot(2,2,4)
    plot(stPSwin{2,5}.dt(:,i), stPSwin{2,5}.depthTCM(:,i))
    ylabel('water depth above TCM (m)')
    plotFancyAxis()
    
    set(hFig(i), 'Position',[1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);
    
    %print('-depsc2', sprintf('f%d_depth.eps' , hFig(i).Number))
end

% accelerometer measurements at Backshore location with histogram
for i = 1 : size(stTCMwin{1,5}.dt, 2)
    
    hFig(i) = figure;
    subplot(4,2,1)
    plot(stPSwin{1,5}.dt(:,i), stPSwin{1,5}.depthTCM(:,i))
    ylabel('water depth above TCM (m)')
    title(sprintf('Backshore - %s', stPSwin{1,5}.dt(1,i)))
    plotFancyAxis()
    
    subplot(4,2,3)
    plot(stTCMwin{1,5}.dt(:,i), stTCMwin{1,5}.Ax(:,i))
    ylabel('X acceleration (g)')
    plotFancyAxis()
    
    subplot(4,2,5)
    plot(stTCMwin{1,5}.dt(:,i), stTCMwin{1,5}.Ay(:,i))
    ylabel('Y acceleration (g)')
    plotFancyAxis()
    
    subplot(4,2,7)
    plot(stTCMwin{1,5}.dt(:,i), stTCMwin{1,5}.Az(:,i))
    ylabel('Z acceleration (g)')
    plotFancyAxis()
    
    subplot(4,2,4)
    hist(stTCMwin{1,5}.Ax(:,i), 10)
    set(gca, 'Xlim', [-2 2])
    ylabel('X acceleration')
    plotFancyAxis()
    
    subplot(4,2,6)
    hist(stTCMwin{1,5}.Ay(:,i), 10)
    set(gca, 'Xlim', [-2 2])
    ylabel('Y acceleration')
    plotFancyAxis()
    
    subplot(4,2,8)
    hist(stTCMwin{1,5}.Az(:,i), 10)
    set(gca, 'Xlim', [-2 2])
    ylabel('Z acceleration')
    plotFancyAxis()
    
    set(hFig(i), 'Position',[1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);
    
    print('-depsc2', sprintf('f%d_depth.eps' , hFig(i).Number))
end

% FSE + speed with quivers
for i = 14%1 : size(stWaveStat{1,5}.FSElpf, 2)
    
    hFig(i) = figure;
    subplot(2,2,1)
    hAx(1) = gca;
    plot(datenum(stWaveStat{1,5}.dtPb(:,i)), stWaveStat{1,5}.FSElpf(:,i))
    ylabel('FSE_{lpf} (m)')
    title(sprintf('Backshore - %s', stWaveStat{1,5}.dtPb(1,i)))
    legend(sprintf('h = %f', stWaveStat{1,5}.h(1,i)))
    plotFancyAxis()

    subplot(2,2,3)
    hAx(2) = gca;
    plot(datenum(stTCMwin{1,5}.dt(:,i)), stTCMwin{1,5}.speedLPF(:,i))
    hold on
    %u = mod(90 - stTCMwin{1,5}.bearingLPF(:,i), 360);
    %u = cos(u) * pi/180 * 0.1;
    u = cosd(90 - stTCMwin{1,5}.bearing(:,i));
    %v = mod(90 - stTCMwin{1,5}.bearingLPF(:,i), 360);
    %v = sin(v) * pi/180 * 0.1;
    v = sind(90 - stTCMwin{1,5}.bearing(:,i));
    un = u ./ sqrt(u.^2 + v.^2) * 0.001;
    vn = v ./ sqrt(u.^2 + v.^2) * 0.001;
    quiver(datenum(stTCMwin{1,5}.dt(1:1000:end,i)), stTCMwin{1,5}.speedLPF(1:1000:end,i),...
           un(1:1000:end), vn(1:1000:end), 'color', ...
           [0 0 0], 'Marker', 'o', 'MarkerFaceColor', 'k','MarkerSize',...
           3, 'ShowArrowHead', 'on', 'Autoscale', 'off')
%     arrow([datenum(stTCMwin{1,5}.dt(1:1000:end,i)) stTCMwin{1,5}.speedLPF(1:1000:end,i)],...
%           [(datenum(stTCMwin{1,5}.dt(1:1000:end,i)) + un(1:1000:end)) ...
%            stTCMwin{1,5}.speedLPF(1:1000:end,i) + vn(1:1000:end)], ...
%           'Length', 4)
    ylabel('speed_{lpf} (m/s)')
    %ylabel('speed (m/s)')
    plotFancyAxis()

    subplot(2,2,2)
    hAx(3) = gca;
    plot(stWaveStat{2,5}.dtPb(:,i), stWaveStat{2,5}.FSElpf(:,i))
    ylabel('FSE_{lpf} (m)')
    title(sprintf('Back-barrier - %s', stWaveStat{2,5}.dtPb(1,i)))
    legend(sprintf('h = %f', stWaveStat{2,5}.h(1,i)))
    plotFancyAxis()

    subplot(2,2,4)
    hAx(4) = gca;
    plot(datenum(stTCMwin{2,5}.dt(:,i)), stTCMwin{2,5}.speedLPF(:,i))
    hold on
    % change from compass directions to standard cartesian (y=0, counter
    % clockwise)
    u = mod(90 - stTCMwin{2,5}.bearingLPF(:,i), 360);
    u = cos(deg2rad(u)) .* stTCMwin{2,5}.speedLPF(:,i);
    v = mod(90 - stTCMwin{2,5}.bearingLPF(:,i), 360);
    v = sin(deg2rad(v)) .* stTCMwin{2,5}.speedLPF(:,i);
    quiver(datenum(stTCMwin{2,5}.dt(1:1000:end,i)), stTCMwin{2,5}.speedLPF(1:1000:end,i),...
           u(1:1000:end), v(1:1000:end), 'color',...
           [0 0 0], 'Marker', 'o', 'MarkerFaceColor', 'k','MarkerSize',...
           3, 'ShowArrowHead', 'off', 'AutoScaleFactor', 0.1)
    ylabel('speed_{lpf} (m/s)')
    %ylabel('speed (m/s)')
    plotFancyAxis()
    
    linkaxes(hAx, 'x')

    set(hFig(i), 'Position', ...
        [1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);
    
    %print('-depsc2', sprintf('f%d.eps' , hFig(i).Number))
    %saveas(hFig(i), sprintf('FSE_TCM-burst%d.fig', i));
end

% FSE and prolate spheroidal PSD for 17 & 34 minute intervals
for i = 1 : size(stWaveStat{1,3}.FSElpf, 2)
    
    hFig(i) = figure;
    subplot(1,2,1)
    plot(stWaveStat{1,3}.dtPb(:,i), stWaveStat{1,3}.FSElpf(:,i))
    ylabel('FSE_{lpf} (m)')
    title(sprintf('Backshore - %s', stWaveStat{1,3}.dtPb(1,i)))
    legend(sprintf('h = %f m', stWaveStat{1,3}.h(1,i)))
    plotFancyAxis()
    
    subplot(1,2,2)
    loglog(stWaveStat{1,3}.f(:,i), stWaveStat{1,3}.psSnnInf(:,i), ...
         stWaveStat{1,3}.f(:,i), stWaveStat{1,3}.psSnnSS(:,i), 'k')
    ylabel('wave power spectra, S_{\eta\eta} (m^2/Hz)')
    xlabel('frequency (Hz)')
    legend(sprintf('df = %f Hz', stWaveStat{1,3}.f(3,i)-stWaveStat{1,3}.f(2,i)))
    plotFancyAxis()

    set(hFig(i), 'Position', ...
        [1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);
    
    print('-depsc2', sprintf('FSE_PSD-burst%d.eps' , hFig(i).Number))
    saveas(hFig(i), sprintf('FSE_PSD-burst%d.fig', i));
end

%% variance of TCM accelerometer measurements 
% 8 minute intervals

nScrsz = get(groot, 'ScreenSize');   % get screen size

hFigTCM = figure;
subplot(3,1,1)
plot(stTCMwin{1,2}.dt(1,:), var(stTCMwin{1,2}.Ax))
     %stTCMwin{2,2}.dt(1,:), var(stTCMwin{2,2}.Ax))
ylabel('\sigma^2 X_{acc} (g)')
title('Backshore accelerometer variance')
plotFancyAxis()
    
subplot(3,1,2)
plot(stTCMwin{1,2}.dt(1,:), var(stTCMwin{1,2}.Ay))
     %stTCMwin{2,2}.dt(1,:), var(stTCMwin{2,2}.Ay))
ylabel('\sigma^2 Y_{acc} (g)')
plotFancyAxis()
    
subplot(3,1,3)
plot(stTCMwin{1,2}.dt(1,:), var(stTCMwin{1,2}.Az))
     %stTCMwin{2,2}.dt(1,:), var(stTCMwin{2,2}.Az))
ylabel('\sigma^2 Z_{acc} (g)')
plotFancyAxis()
    
set(hFigTCM, 'Position',[1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);

hFigTCM = figure;
subplot(3,1,1)
plot(stTCMwin{2,2}.dt(1,:), var(stTCMwin{2,2}.Ax))
ylabel('\sigma^2 X_{acc} (g)')
title('Back-barrier accelerometer variance')
plotFancyAxis()
    
subplot(3,1,2)
plot(stTCMwin{2,2}.dt(1,:), var(stTCMwin{2,2}.Ay))
ylabel('\sigma^2 Y_{acc} (g)')
plotFancyAxis()
    
subplot(3,1,3)
plot(stTCMwin{2,2}.dt(1,:), var(stTCMwin{2,2}.Az))
ylabel('\sigma^2 Z_{acc} (g)')
plotFancyAxis()
    
set(hFigTCM, 'Position',[1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);

%% check Nick Lowell's gui codes vs. my code for tilt and bearing

% RRU-1 (minor differences, likely do to interpolation algorithm and zero
% point vector; remember that my codes don't allow for tilt > 90 degrees 
% whereas Nick's do)
iPOD = 1;
figure;
plot(stPOD(iPOD).TCMdt, stPOD(iPOD).velocityNcalc, ...
     stPOD(iPOD).TCMdt, stPOD(iPOD).velocityN/100)
ylabel('velocity North (m/s)')
legend('calculated', 'Matlogger')

temp = mod(stPOD(iPOD).yaw, pi);
figure;
plot(stPOD(iPOD).TCMdt, rad2deg(stPOD(iPOD).point), ...
     stPOD(iPOD).TCMdt, stPOD(iPOD).bearing)
ylabel('bearing (deg)')
legend('calculated', 'Matlogger')

%% WSE + mean water level + precipitation + ADCIRC comparison (AGU FIGURE)

cmap = parula(5);

% WSEs
figure
subplot(2,1,1)
plot(stPOD(1).PSdt, stPOD(1).wse, 'Color', cmap(1,:))
hold on
plot(stPOD(2).PSdt, stPOD(2).wse, 'Color', cmap(3,:))
hold on
% ADCIRC WSEs (SWE + tidal residual)
plot(rTimeADCIRC, rWseADCIRC(2471503,:), 'k')
hold on
plot(rTimeADCIRC, rWseADCIRC(2419870,:), 'k')
hold on
plot(rTimeADCIRC, rWseADCIRC(2471504,:), 'k')
% hold on
% 8 minute mean water levels 
% plot(stWaveStat{2,2}.dtPb(1,:), stWaveStat{2,2}.h(1,:) + 0.260, 'k')
% hold on
% plot(stWaveStat{1,2}.dtPb(1,:), stWaveStat{1,2}.h(1,:) + 0.9315, 'k')
legend('Backshore', 'Back-barrier')
ylabel('water surface elevation (m NAVD88)')
plotFancyAxis()

% rainfall
subplot(2,1,2)
plot(dtUSGSprecip, rUSGSprecip)
ylabel('precipitation (inches)')
plotFancyAxis()

%% Significant wave height for RRU-1, RRU-2 (AGU FIGURE)

% decided to use multi-sine tapers here because the Hs at landfall seemed
% more accurate (based on visual evaluation of wave height ~0.19 m); only
% plot back-barrier between 16:00 and 22:00

% landfall datetime series for differencing 
dtLandRRU1 = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(stWaveStat{1,2}.dtPb(1,:),2), 1);
dtLandRRU2 = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(stWaveStat{2,2}.dtPb(1,:),2), 1);
dtLandADs  = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(rTimeADCIRC,1), 1);                
dtLandRRU1 = hours(stWaveStat{1,2}.dtPb(1,:)' - dtLandRRU1);
dtLandRRU2 = hours(stWaveStat{2,2}.dtPb(1,:)' - dtLandRRU2);
dtLandADs  = hours(rTimeADCIRC - dtLandADs);

figure
yyaxis left
plot(dtLandRRU1, stWaveStat{1,2}.msHm0Inf, ...
     dtLandRRU2(92:134), stWaveStat{2,2}.msHm0Inf(92:134))
     %rTimeADCIRC, rHsADCIRC(2903359,:))
ylabel('significant wave height, H_{m0} (m)')
hold on
yyaxis right
plot(dtLandRRU1, stWaveStat{1,2}.msTm01Inf, ...
     dtLandRRU2(92:134), stWaveStat{2,2}.msTm01Inf(92:134), ...
     dtLandADs, rTm0ADCIRC(2471504,:)) 
ylabel('mean wave period, T_{m0} (m)')
legend('backshore', 'back-barrier')
plotFancyAxis()

%% Water depth above PT + LPF, speed, and arrows (8 min records) [AGU FIG]

for i = 12%1 : size(stPSwin{1,5}.dt, 2)
    
%     hFig(i) = figure;
%     subplot(2,1,1)
%     yyaxis right
%     plot(stPSwin{1,5}.dt(:,i), stPSwin{1,5}.depthPT(:,i))
%     %hold on
%     %plot(stPSwin{1,5}.dt(:,i), stPSwin{1,5}.depthPTlpf(:,i))
%     ylabel('water depth (m)')
%     %legend('raw', 'lpf')
%     title(sprintf('Backshore - %s', stPSwin{1,5}.dt(1,i)))
%     plotFancyAxis()
%     
%     %subplot(2,1,2)
%     yyaxis left
%     plot(stTCMwin{1,5}.dt(:,i), stTCMwin{1,5}.speed(:,i)/100)
%     hold on
%     plot(stTCMwin{1,5}.dt(:,i), stTCMwin{1,5}.speedLPF(:,i), 'k')
% %    hold on
% %     u = cosd(90 - stTCMwin{1,5}.bearing(:,i));
% %     v = sind(90 - stTCMwin{1,5}.bearing(:,i));
% %     un = u ./ sqrt(u.^2 + v.^2) * 0.001;
% %     vn = v ./ sqrt(u.^2 + v.^2) * 0.001;
% %     quiver(datenum(stTCMwin{1,5}.dt(1:1000:end,i)), stTCMwin{1,5}.speedLPF(1:1000:end,i),...
% %            un(1:1000:end), vn(1:1000:end), 'color', ...
% %            [0 0 0], 'Marker', 'o', 'MarkerFaceColor', 'k','MarkerSize',...
% %            3, 'ShowArrowHead', 'on', 'LineStyle', 'none')
%     ylabel('speed (m/s)')
%     plotFancyAxis()
    
    % RRU-2
    hFig(i) = figure;
    subplot(2,1,1)
    yyaxis right
    plot(stPSwin{2,5}.dt(:,i), stPSwin{2,5}.depthPT(:,i))
    hold on
    plot(stPSwin{2,5}.dt(:,i), stPSwin{2,5}.depthPTlpf(:,i), 'k')
    ylabel('water depth (m)')
    legend('raw', 'lpf')
    title(sprintf('Back-barrier - %s', stPSwin{2,5}.dt(1,i)))
    plotFancyAxis()
    
    %subplot(2,1,2)
    yyaxis left
    plot(stTCMwin{2,5}.dt(:,i), stTCMwin{2,5}.velocityN(:,i)/100)
    hold on
    plot(stTCMwin{2,5}.dt(:,i), stTCMwin{2,5}.velNLPF(:,i), 'k')
%     hold on
%     u = cosd(90 - stTCMwin{2,5}.bearing(:,i));
%     v = sind(90 - stTCMwin{2,5}.bearing(:,i));
%     un = stTCMwin{2,5}.velELPF(:,i) ./ sqrt(stTCMwin{2,5}.velELPF(:,i).^2 + stTCMwin{2,5}.velNLPF(:,i).^2) * 0.001;
%     vn = stTCMwin{2,5}.velNLPF(:,i) ./ sqrt(stTCMwin{2,5}.velELPF(:,i).^2 + stTCMwin{2,5}.velNLPF(:,i).^2) * 0.001;
%     quiver(datenum(stTCMwin{2,5}.dt(1:1000:end,i)), stTCMwin{2,5}.speedLPF(1:1000:end,i),...
%            un(1:1000:end), vn(1:1000:end), 'color', ...
%            [0 0 0], 'Marker', 'o', 'MarkerFaceColor', 'k','MarkerSize',...
%            3, 'ShowArrowHead', 'off', 'AutoScale', 'off', 'LineStyle', 'none')
%     ylabel('speed_{lpf} (m/s)')
    ylabel('speed (m/s)')
    plotFancyAxis()
    
    subplot(2,1,2)
    yyaxis right
    plot(stPSwin{2,5}.dt(:,i), stPSwin{2,5}.depthPT(:,i))
    hold on
    plot(stPSwin{2,5}.dt(:,i), stPSwin{2,5}.depthPTlpf(:,i), 'k')
    ylabel('water depth (m)')
    legend('raw', 'lpf')
    title(sprintf('Back-barrier - %s', stPSwin{2,5}.dt(1,i)))
    plotFancyAxis()
    
    %subplot(2,1,2)
    yyaxis left
    plot(stTCMwin{2,5}.dt(:,i), stTCMwin{2,5}.velocityE(:,i)/100)
    hold on
    plot(stTCMwin{2,5}.dt(:,i), stTCMwin{2,5}.velELPF(:,i), 'k')
%     hold on
%     u = cosd(90 - stTCMwin{2,5}.bearing(:,i));
%     v = sind(90 - stTCMwin{2,5}.bearing(:,i));
%     un = stTCMwin{2,5}.velELPF(:,i) ./ sqrt(stTCMwin{2,5}.velELPF(:,i).^2 + stTCMwin{2,5}.velNLPF(:,i).^2) * 0.001;
%     vn = stTCMwin{2,5}.velNLPF(:,i) ./ sqrt(stTCMwin{2,5}.velELPF(:,i).^2 + stTCMwin{2,5}.velNLPF(:,i).^2) * 0.001;
%     quiver(datenum(stTCMwin{2,5}.dt(1:1000:end,i)), stTCMwin{2,5}.speedLPF(1:1000:end,i),...
%            un(1:1000:end), vn(1:1000:end), 'color', ...
%            [0 0 0], 'Marker', 'o', 'MarkerFaceColor', 'k','MarkerSize',...
%            3, 'ShowArrowHead', 'off', 'AutoScale', 'off', 'LineStyle', 'none')
%     ylabel('speed_{lpf} (m/s)')
    ylabel('speed (m/s)')
    plotFancyAxis()
    
    set(hFig(i), 'Position',[1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);
    
    %print('-depsc2', sprintf('f%d_depth.eps' , hFig(i).Number))
end

%% ADCIRC comparisons

% Follets West
figure
% offshore
subplot(2,1,1)
plot(stPOD(4).PSdt, stPOD(4).wse, 'Color', cmap(2,:))
hold on
plot(stPOD(3).PSdt, stPOD(3).wse, 'Color', cmap(1,:))
hold on
% ADCIRC WSEs (SWE + tidal residual)
plot(rTimeADCIRC, rWseADCIRC(2986162,:), '-.k')
hold on
plot(rTimeADCIRC, rWseADCIRC(2904365,:), '-.k')
% formatting
legend('mid-barrier: 0.50 m', 'back-barrier: 0.190 m', ...
       'ADCIRC node 2986162: 0.773 m', 'ADCIRC node 2904365: 0.277 m')
ylabel('water surface elevation (m NAVD88)')
plotFancyAxis()

% onshore
subplot(2,1,2)
plot(stVECT.dt, stVECT.wse, 'Color', cmap(3,:))
hold on
plot(stVECT.dt, stVECT.swe, 'k')
hold on
% ADCIRC WSE (SWE + tidal residual)
plot(rTimeADCIRC, rWseADCIRC(3014219,:), '-.k')
legend('ADV: -0.75 m', 'ADV SWE', 'ADCIRC node 3014219: -1.45 m')
ylabel('water surface elevation (m NAVD88)')
plotFancyAxis() 

% Matagorda Peninsula
figure
subplot(1,2,1)
yyaxis left
plot(rTimeADCIRC, rWseADCIRC(2471504,:))
ylabel('water surface elevation (m NAVD88)')
yyaxis right
plot(rTimeADCIRC, rTm0ADCIRC(2471504,:))
ylabel('mean wave period (s)') 
title('ADCIRC node 2471504 - surf zone')

subplot(1,2,2)
yyaxis left
plot(rTimeADCIRC, rHsADCIRC(2471504,:))
ylabel('significant wave height, H_{m0} (m)')
yyaxis right
plot(rTimeADCIRC, rSwanDirADCIRC(2471504,:))
ylabel('wave direction (from deg)')

%% NDBC wave spectra and wave statistics

    %%% NOTE: (2/2019) both the NDBC and ADCIRC datetimes were off when
    %%% this analysis was originally completed (ADCIRC by one hour, NDBC
    %%% by 6 hours [UTC]) - so much of this code is incorrect but keeping 
    %%% for reference)
    
% find y limits
nYmin = min(stNDBC.psd(:));
nYmax = max(stNDBC.psd(:));

%figure;
% only make plots for 8/25/17 03:00 AM to 08/26/17 06:00
for i = 577 : 604
    
    hFig(i) = figure; 
    plot(stNDBC.f(:,i), stNDBC.psd(:,i)) 
    xlabel('frequency [Hz]')
    ylabel('S_{\eta\eta} [m^2/Hz]')
    ylim([nYmin nYmax])
    title(sprintf('NDBC wave spectra: %s', stNDBC.dt(:,i)))
    plotFancyAxis
    %hold on
    saveas(hFig(i), sprintf('NDBCspectra-%d.fig', i));
    
end

% offshore significant wave height and mean period
figure
yyaxis left
plot(stNDBC.dt(577:604), stNDBC.Hm0SS(577:604))
ylabel('significant wave height, H_{m0}- ss (m)')
hold on
yyaxis right
plot(stNDBC.dt(577:604), stNDBC.Tm01SS(577:604)) 
ylabel('mean wave period, T_{m0}- ss (s)')

% offshore significant wave height and peak period
figure
yyaxis left
plot(stNDBC.dt(577:604), stNDBC.Hm0SS(577:604))
ylabel('significant wave height, H_{m0}- ss (m)')
hold on
yyaxis right
plot(stNDBC.dt(577:605), 1 ./ stNDBC.FpeakSS(577:605)) 
ylabel('peak wave period, T_{p}- ss (s)')

% spectogram
figure
surf(datenum(repmat(stNDBC.dt(577:604), size(stNDBC.f,1), 1)), ...
     stNDBC.f(:,577:604), stNDBC.psd(:,577:604), 'EdgeColor', 'none')
axis xy; view(2)
datetick('x','mm/dd HH PM')
axis tight
% formatting
c = colorbar;
ylabel(c, 'Surface elevation variance, S_{\eta\eta}(f) [m^2/Hz]','FontSize', 13);
ylabel('Frequency (Hz)');
title('NDBC Wave ');
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0.05 0.25]) 

%% change point wave statistics vs full time series (MTM only)

for iPOD = 1 : 2
    figure;
    subplot(4,1,1)
    % compare against regular, not zero padded
    plot(stMTMWaveStatCP{2,iPOD}.dtPb(1,:), stMTMWaveStatCP{2,iPOD}.Hm0Inf(1,:), ...
         stMTMWaveStat{2,iPOD}.dtPb(1,:), stMTMWaveStat{2,iPOD}.Hm0Inf(1,:));
    ylabel('Hm0_{inf} [m]')
    legend('CP', 'full')
    plotFancyAxis
    title('~17 min wave statistics (MTM only)')

    subplot(4,1,2)
    plot(stMTMWaveStatCP{2,iPOD}.dtPb(1,:), stMTMWaveStatCP{2,iPOD}.Tm01Inf(1,:), ...
         stMTMWaveStat{2,iPOD}.dtPb(1,:), stMTMWaveStat{2,iPOD}.Tm01Inf(1,:));
    plotFancyAxis
    ylabel('Tm01_{inf} (s)')

    subplot(4,1,3)
    plot(stMTMWaveStatCP{2,iPOD}.dtPb(1,:), stMTMWaveStatCP{2,iPOD}.Tm02Inf(1,:), ...
         stMTMWaveStat{2,iPOD}.dtPb(1,:), stMTMWaveStat{2,iPOD}.Tm02Inf(1,:));
    plotFancyAxis
    ylabel('Tm02_{inf} (s)')

    subplot(4,1,4)
    plot(stMTMWaveStatCP{2,iPOD}.dtPb(1,:), 1 ./ stMTMWaveStatCP{2,iPOD}.FpeakInf(1,:), ...
         stMTMWaveStat{2,iPOD}.dtPb(1,:), 1 ./ stMTMWaveStat{2,iPOD}.FpeakInf(1,:));
    plotFancyAxis
    ylabel('Tp_{inf} (s)')
end

% dlmwrite('CPpsin.txt',rCPwin,'delimiter','\t','precision',12)

%% comparison of mtmWaveStatCP 17 min wave spectra between RRU-1 and RRU-2

%%% plots wave statistics from the outputs from Thomson's prolate 
%%% spheroidal taper spectra exported from R, which were derived using
%%% the ~17 min change point identified spectra for each ~34 min win
%%% (stMTMps); note that the first row of stMTMps is the zero padded spectra 
%%% and the second row is the standard nfft; K = 4 and NW = 3

for iBurst = 25 : 34%length(stCPfin{1,2}.Hm0Inf)
    
    hFig = figure;
    
%     % plot each of the burst time series
%     subplot(1,3,1)
%     plot(stMTMWaveStatCP{1,1}.Pb(:,iBurst)/(nRho * nG))
%     ylabel('water surface elevation (m)')
%     title('RRU-1')
%     plotFancyAxis     
%     
%     subplot(1,3,2)
%     plot(stMTMWaveStatCP{1,2}.Pb(:,iBurst)/(nRho * nG))
%     ylabel('water surface elevation (m)')
%     title('RRU-2')
%     plotFancyAxis
    
    % plot normalized power spectra (zero padded)
    %subplot(1,3,3)
    loglog(stMTMWaveStatCP{1,1}.f(:,iBurst), stMTMWaveStatCP{1,1}.SnnNormInf(:,iBurst), ...
           stMTMWaveStatCP{1,2}.f(:,iBurst), stMTMWaveStatCP{1,2}.SnnNormInf(:,iBurst))
    ylabel('S_{\eta\eta}(f) [arb]')
    legend('RRU-1', 'RRU-2')
    xlabel('frequency (Hz)')
    set(gca, 'XLim', [0.005 0.04])
    plotFancyAxis
    title(sprintf('%s, Burst %d', stMTMWaveStatCP{1,1}.dtPb(1,iBurst), iBurst))
    
    %set(hFig, 'Position',[1 nScrsz(4)*3/4 nScrsz(3)*3/4 nScrsz(4)*3/4]);
end

%% Comparison of psd vs mtm 17 min wave stats

for iPOD = 1 : 2
    figure;
    subplot(4,1,1)
    % compare against regular, not zero padded
    plot(stWaveStat{iPOD,3}.dtPb(1,:), stWaveStat{iPOD,3}.psHm0Inf(1,:), ...
           stMTMWaveStat{2,iPOD}.dtPb(1,:), stMTMWaveStat{2,iPOD}.Hm0Inf(1,:));
    ylabel('Hm0_{inf} [m]')
    legend('psd', 'mtm')
    plotFancyAxis
    title('~17 min wave statistics')

    subplot(4,1,2)
    plot(stWaveStat{iPOD,3}.dtPb(1,:), stWaveStat{iPOD,3}.psTm01Inf(1,:), ...
         stMTMWaveStat{2,iPOD}.dtPb(1,:), stMTMWaveStat{2,iPOD}.Tm01Inf(1,:));
    plotFancyAxis
    ylabel('Tm01_{inf} (s)')

    subplot(4,1,3)
    plot(stWaveStat{iPOD,3}.dtPb(1,:), stWaveStat{iPOD,3}.psTm02Inf(1,:), ...
         stMTMWaveStat{2,iPOD}.dtPb(1,:), stMTMWaveStat{2,iPOD}.Tm02Inf(1,:));
    plotFancyAxis
    ylabel('Tm02_{inf} (s)')

    subplot(4,1,4)
    plot(stWaveStat{iPOD,3}.dtPb(1,:), 1 ./ stWaveStat{iPOD,3}.psFpeakInf(1,:), ...
         stMTMWaveStat{2,iPOD}.dtPb(1,:), 1 ./ stMTMWaveStat{2,iPOD}.FpeakInf(1,:));
    plotFancyAxis
    ylabel('Tp_{inf} (s)')
end

%% Significant wave height for RRUs and ADV (OSM & THESIS FIGURES)

    %%% NOTE: (2/2019) both the NDBC and ADCIRC datetimes were off when
    %%% this analysis was originally completed (ADCIRC by one hour, NDBC
    %%% by 6 hours [UTC]) - so much of this code is incorrect but keeping 
    %%% for reference)

% decided to use prolate-spheroidal tapers (MTM) here because the spectral 
% estimate is superior; only plot back-barrier between 16:00 and 22:00

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS
load('mtmWaveStats_17min_0pt003_nw3k4.mat') 
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/NDBC
load('NDBC.mat') 
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ADCIRC/Hindcast43-NHC_best_track
load('ADCIRC.mat') 
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM
load('TIDE.mat') % SLP predicted tides (m)

% landfall datetime series for differencing (zero-padded)
dtLandRRU1 = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(stMTMWaveStat0pt003{1,1}.dtPb(1,:),2), 1);
% dtLandRRU1cp = repmat(datetime('25-Aug-2017 22:00:00'), ...
%                     size(stMTMWaveStatCP0pt003{1,1}.dtPb(1,:),2), 1);                
dtLandRRU2 = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(stMTMWaveStat0pt003{1,2}.dtPb(1,:),2), 1);
dtLandADV  = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(stMTMWaveStat0pt003{1,3}.dtPb(1,:),2), 1);    
dtLandADs  = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(rTimeADCIRC,1), 1); 
dtLandNDBC = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(stNDBC{1,1}.dt,2), 1); 
dtLandTide = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(stTide.dt,2), 1);                
dtLandRRU1 = hours(stMTMWaveStat0pt003{1,1}.dtPb(1,:)' - dtLandRRU1);
% dtLandRRU1cp = hours(stMTMWaveStatCP0pt003{1,1}.dtPb(1,:)' - dtLandRRU1cp);
dtLandRRU2 = hours(stMTMWaveStat0pt003{1,2}.dtPb(1,:)' - dtLandRRU2);
dtLandADV  = hours(stMTMWaveStat0pt003{1,3}.dtPb(1,:)' - dtLandADV);
dtLandADs  = hours(rTimeADCIRC - dtLandADs);
dtLandNDBC = hours(stNDBC{1,1}.dt' - dtLandNDBC);
dtLandTide = hours(stTide.dt' - dtLandTide);

% RRU1 & RRU2
iBlow = 46;
iBup  = 68;

figure
yyaxis left
plot(dtLandRRU1cp(1:23), stMTMWaveStatCP0pt003{1,1}.Hm0Inf(1:23), '.-')
hold on
plot(dtLandRRU1(iBlow:iBup), stMTMWaveStat0pt003{1,1}.Hm0Inf(iBlow:iBup)) 
hold on
plot(dtLandRRU1cp(35:end), stMTMWaveStatCP0pt003{1,1}.Hm0Inf(35:end), '.-')
hold on
plot(dtLandRRU2(iBlow:iBup), stMTMWaveStat0pt003{1,2}.Hm0Inf(iBlow:iBup))
ylabel('significant wave height, Hm0_{inf} (m)')
hold on
yyaxis right
plot(dtLandRRU1cp(1:23), stMTMWaveStatCP0pt003{1,1}.Tm01Inf(1:23), '.-')
hold on
plot(dtLandRRU1(iBlow:iBup), stMTMWaveStat0pt003{1,1}.Tm01Inf(iBlow:iBup)) 
hold on
plot(dtLandRRU1cp(35:end), stMTMWaveStatCP0pt003{1,1}.Tm01Inf(35:end), '.-')
hold on
plot(dtLandRRU2(iBlow:iBup), stMTMWaveStat0pt003{1,2}.Tm01Inf(iBlow:iBup))
ylabel('mean wave period, Tm01_{inf} ()')
xlabel('hours relative to hurricane landfall')
plotFancyAxis

% ADV (infragravity only)
figure
subplot(2,1,1)
yyaxis left
plot(dtLandADV, stMTMWaveStat0pt003{1,3}.Hm0Inf)
ylabel('Hm0_{inf} (m)')
hold on
yyaxis right
plot(dtLandADV, stMTMWaveStat0pt003{1,3}.Tm01Inf)
ylabel('Tm01_{inf} (s)')
set(gca, 'xlim', [-22 66]) 
plotFancyAxis

subplot(2,1,2)
% plot the predicted tide and SWE
nVECTele   = -0.746; % ADV elevation 
nMSL       = 0.1;  % San Luis Pass MSL
plot(dtLandTide, stTide.pred, 'k', 'LineWidth', 2)
ylabel('\eta (m NAVD88)')
hold on
plot(dtLandADV, movmean(stMTMWaveStat0pt003{1,3}.h+nVECTele, 5), 'LineWidth', 2)
hold on 
plot(dtLandTide, repmat(nMSL, size(stTide.dt)), 'LineWidth', 2)
set(gca, 'xlim', [-22 66]) 
%plot(stADVWaveStat{1,5}.dtPb(1,:), stADVWaveStat{1,5}.h)
legend('predicted tide', 'measured SWE', 'MSL')
xlabel('hours relative to hurricane landfall')
plotFancyAxis

% ADV infragravity and offshore (NDBC) sea-swell for landfall only
iBlow = 1;
iBup = 118;

figure
subplot(2,1,1)
yyaxis left
plot(dtLandADV(iBlow:iBup), stMTMWaveStat0pt003{1,3}.Hm0Inf(iBlow:iBup))
ylabel('significant wave height, Hm0_{inf} (m)')
hold on
yyaxis right
plot(dtLandADV(iBlow:iBup), stMTMWaveStat0pt003{1,3}.Tm01Inf(iBlow:iBup))
ylabel('mean wave period, Tm01_{inf} (s)')
plotFancyAxis

% NDBC
iBlow = 573;
iBup = 607;

subplot(2,1,2)
yyaxis left
plot(dtLandNDBC(iBlow:iBup), stNDBC{1,1}.Hm0SS(iBlow:iBup))
ylabel('offshore wave height, H_{0} (m)')
hold on
yyaxis right
plot(dtLandNDBC(iBlow:iBup), stNDBC{1,1}.Tm01SS(iBlow:iBup)) 
ylabel('offshore wave period, T_{0} (s)')
plotFancyAxis

% ADV (infragravity and sea-swell)
figure
subplot(3,1,1)
yyaxis left
plot(dtLandADV, stMTMWaveStat0pt003{1,3}.Hm0Inf)
ylabel('Hm0_{inf} (m)')
hold on
yyaxis right
plot(dtLandADV, stMTMWaveStat0pt003{1,3}.Tm01Inf)
%hold on
%plot(dtLandADV, 1./stMTMWaveStat0pt003{1,3}.FpeakInf)
ylabel('Tm0_{inf} (sec)')
legend('Hm0', 'Tm01')
plotFancyAxis
set(gca, 'xlim', [-22 66]) 

subplot(3,1,2)
yyaxis left
plot(dtLandADV, stMTMWaveStat0pt003{1,3}.Hm0SS)
ylabel('Hm0_{ss} (sec)')
hold on
yyaxis right
plot(dtLandADV, stMTMWaveStat0pt003{1,3}.Tm01SS)
%hold on
%semilogy(dtLandADV, 1./stMTMWaveStat0pt003{1,3}.FpeakSS)
ylabel('Tm0_{ss} (sec)')
legend('Hm0', 'Tm01')
plotFancyAxis
set(gca, 'xlim', [-22 66]) 

subplot(3,1,3)
% plot the predicted tide and SWE
nVECTele   = -0.746; % ADV elevation 
nMSL       = 0.1;  % San Luis Pass MSL

plot(dtLandTide, stTide.pred, 'k', 'LineWidth', 2)
ylabel('\eta (m NAVD88)')
hold on
plot(dtLandADV, movmean(stMTMWaveStat0pt003{1,3}.h+nVECTele, 5), 'LineWidth', 2)
hold on 
plot(dtLandTide, repmat(nMSL, size(stTide.dt)), 'LineWidth', 2)
set(gca, 'xlim', [-22 66]) 
%plot(stADVWaveStat{1,5}.dtPb(1,:), stADVWaveStat{1,5}.h)
legend('predicted tide', 'measured SWE', 'MSL')
xlabel('hours relative to hurricane landfall')
plotFancyAxis

% ADV vs. ADCIRC for sea-swell band
iADoff = 3083010;
iADosz = 3001063;
iADisz = 3014219;
% ADCIRC time bounds
iADClow  = 41;
iADCup   = 130;

% Hs
figure
subplot(4,1,1)
semilogy(dtLandADV, stMTMWaveStat0pt003{1,3}.Hm0SS, ...
     dtLandADs(iADClow:iADCup), rHsADCIRC(iADoff,iADClow:iADCup), ...
     dtLandADs(iADClow:iADCup), rHsADCIRC(iADosz,iADClow:iADCup), ...
     dtLandADs(iADClow:iADCup), rHsADCIRC(iADisz,iADClow:iADCup))
ylabel('Hm0_{ss} (m)')
legend('ADV: -0.75 m', 'ADCIRC offshore', ...
       'ADCIRC outer sz', 'ADCIRC inner sz')
plotFancyAxis() 
set(gca, 'xlim', [-22 66]) 

% Tm0
subplot(4,1,2)
yyaxis left
plot(dtLandADV, stMTMWaveStat0pt003{1,3}.Tm01SS)
ylabel('Tm01_{ss} (m)')
yyaxis right
plot(dtLandADs(iADClow:iADCup), rTm0ADCIRC(iADoff,iADClow:iADCup), ...
     dtLandADs(iADClow:iADCup), rTm0ADCIRC(iADosz,iADClow:iADCup), ...
     dtLandADs(iADClow:iADCup), rTm0ADCIRC(iADisz,iADClow:iADCup))
ylabel('T0_{ss} (m)')
legend('ADV: -0.75 m', 'ADCIRC offshore', ...
       'ADCIRC outer sz', 'ADCIRC inner sz')
plotFancyAxis() 
set(gca, 'xlim', [-22 66]) 

subplot(4,1,3)
%plot(dtLandADV, movmean(1./stMTMWaveStat0pt003{1,3}.FpeakSS,20), ...
plot(dtLandADV, 1./stMTMWaveStat0pt003{1,3}.FpeakSS, ...
     dtLandADs(iADClow:iADCup), rTpsADCIRC(iADoff,iADClow:iADCup), ...
     dtLandADs(iADClow:iADCup), rTpsADCIRC(iADosz,iADClow:iADCup), ...
     dtLandADs(iADClow:iADCup), rTpsADCIRC(iADisz,iADClow:iADCup))
ylabel('Tp_{ss} (m)')
plotFancyAxis() 
set(gca, 'xlim', [-22 66]) 

% surge (zeta)
subplot(4,1,4)
% plot the predicted tide and SWE
nVECTele   = -0.746; % ADV elevation 
%nMSL      = 0.1;  % San Luis Pass MSL
dtSS       = rTimeADCIRC(iADClow:iADCup);

rTide = interp1(datenum(stTide.dt), stTide.pred, datenum(dtSS));
rSWE  = interp1(datenum(stMTMWaveStat0pt003{1,3}.dtPb(1,:)), ...
                movmean(stMTMWaveStat0pt003{1,3}.h+nVECTele, 5), ...
                datenum(dtSS)); 
rSurgeADV = rSWE - rTide;

% surge from ADCIRC
rSurgeOFF = rWseADCIRC(iADoff,iADClow:iADCup)-rTide';
rSurgeOSZ = rWseADCIRC(iADosz,iADClow:iADCup)-rTide';
rSurgeISZ = rWseADCIRC(iADisz,iADClow:iADCup)-rTide';

plot(dtLandADs(iADClow:iADCup), rSurgeADV, ...
     dtLandADs(iADClow:iADCup), rSurgeOFF, dtLandADs(iADClow:iADCup), ...
     rSurgeOSZ, dtLandADs(iADClow:iADCup), rSurgeISZ, ...
     'LineWidth', 2)
ylabel('\zeta (m')
set(gca, 'xlim', [-22 66]) 
legend('ADV: -0.75 m', 'ADCIRC offshore', ...
       'ADCIRC outer sz', 'ADCIRC inner sz')
xlabel('hours relative to hurricane landfall')
plotFancyAxis

% NDBC and ADCIRC for Hs and T
iADoff = 3083010;
iADosz = 3001063;
iADisz = 3014219;
% ADCIRC time bounds
iADClow= 41;
iADCup = 130;
% NDBC time bounds
iBlow  = 573;
iBup   = 662;

% Hs
figure;
subplot(2,1,1)
plot(dtLandNDBC(iBlow:iBup), stNDBC{1,1}.Hm0SS(iBlow:iBup), ...
     dtLandADs(iADClow:iADCup), rHsADCIRC(iADoff,iADClow:iADCup), ...
     dtLandADs(iADClow:iADCup), rHsADCIRC(iADosz,iADClow:iADCup), 'LineWidth', 2)
ylabel('H_0 (m)')
legend('DW', 'OF', 'OSZ')
plotFancyAxis() 
set(gca, 'xlim', [-23 66]) 

subplot(2,1,2)
plot(dtLandNDBC(iBlow:iBup), stNDBC{1,1}.Tm01SS(iBlow:iBup), ...
     dtLandADs(iADClow:iADCup), rTm0ADCIRC(iADoff,iADClow:iADCup), ...
     dtLandADs(iADClow:iADCup), rTm0ADCIRC(iADosz,iADClow:iADCup), 'LineWidth', 2)
hold on
plot(dtLandADs(iADClow:iADCup), repmat((1/0.14),iADCup-iADClow+1,1))
ylabel('T_{0} (s)')
xlabel('hours relative to hurricane landfall')
plotFancyAxis
set(gca, 'xlim', [-23 66]) 

%% s for RRUs and ADV (infragravity band only) 

% again, use zero-padded data
% ADV ---------------------------------------------------------------------
figure
% subplot(1,2,1)
surf(repmat(dtLandADV', size(stMTMWaveStat{1,3}.f,1), 1), ...
     stMTMWaveStat{1,3}.f, stMTMWaveStat{1,3}.SnnInf, 'EdgeColor', 'none')
axis xy; view(2)
%datetick('x','mm/dd HH PM')
axis tight
% formatting
c = colorbar;
ylabel(c, 'water surface elevation variance, S_{\eta\eta}(f) [m^2/Hz]', ...
                'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0.005 0.04]) 

% RRU1 --------------------------------------------------------------------
figure
% subplot(1,2,1)
surf(repmat(dtLandRRU1', size(stMTMWaveStat0pt003{1,1}.f,1), 1), ...
     stMTMWaveStat0pt003{1,1}.f, stMTMWaveStat0pt003{1,1}.SnnInf, 'EdgeColor', 'none')
axis xy; view(2)
%datetick('x','mm/dd HH PM')
axis tight
% formatting
c = colorbar;
ylabel(c, 'water surface elevation variance, S_{\eta\eta}(f) [m^2/Hz]', ...
                'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0.005 0.04]) 

% RRU2 --------------------------------------------------------------------
figure
% subplot(1,2,1)
surf(repmat(dtLandRRU2', size(stMTMWaveStat0pt003{1,2}.f,1), 1), ...
     stMTMWaveStat0pt003{1,2}.f, stMTMWaveStat0pt003{1,2}.SnnInf, 'EdgeColor', 'none')
axis xy; view(2)
%datetick('x','mm/dd HH PM')
axis tight
% formatting
c = colorbar;
ylabel(c, 'water surface elevation variance, S_{\eta\eta}(f) [m^2/Hz]', ...
                'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0.005 0.04])

%% Wave frequency dissipation (OSM Figure)

% bursts corresponding to 16:00 and 22:00 (I looked these up)
% RRUs -------------------------------------------------------------------
iBlow = 46;
iBup  = 68;
    
figure;
subplot(4,1,1)
    plot(dtLandRRU1(iBlow:iBup), ...
     1 ./ stMTMWaveStat{1,1}.Tm01Inf(1,iBlow:iBup), '-k', ...
     dtLandRRU2(iBlow:iBup), ...
     1 ./ stMTMWaveStat{1,2}.Tm01Inf(1,iBlow:iBup), '--k')
ylabel('mean frequency, 1/Tm01 (Hz)')
legend('backshore', 'back-barrier')
plotFancyAxis

iRRU = 1;
subplot(4,1,2)
plot(dtLandRRU1(iBlow:iBup), ...
     stMTMWaveStat{1,iRRU}.infLow(1,iBlow:iBup), '-o', ...
     dtLandRRU1(iBlow:iBup), ...
     stMTMWaveStat{1,iRRU}.infMid(1,iBlow:iBup), '-o', ...
     dtLandRRU1(iBlow:iBup), ...
     stMTMWaveStat{1,iRRU}.infUp(1,iBlow:iBup), '-o')
ylabel('% S\eta\eta_{inf} (arb)')
legend('0.005 < f < 0.017', '0.017 < f < 0.028', '0.028 < f < 0.04')
title('backshore')
plotFancyAxis

iRRU = 2;
subplot(4,1,3)
plot(dtLandRRU2(iBlow:iBup), ...
     stMTMWaveStat{1,iRRU}.infLow(1,iBlow:iBup), '-o', ...
     dtLandRRU2(iBlow:iBup), ...
     stMTMWaveStat{1,iRRU}.infMid(1,iBlow:iBup), '-o', ...
     dtLandRRU2(iBlow:iBup), ...
     stMTMWaveStat{1,iRRU}.infUp(1,iBlow:iBup), '-o')
ylabel('% S\eta\eta_{inf} (arb)')
title('back-barrier')
plotFancyAxis

subplot(4,1,4)
    plot(dtLandRRU1(iBlow:iBup), ...
     stMTMWaveStat{1,1}.h(1,iBlow:iBup), '-k', ...
     dtLandRRU2(iBlow:iBup), ...
     stMTMWaveStat{1,2}.h(1,iBlow:iBup), '--k')
ylabel('mean water depth (m)')
%legend('backshore', 'back-barrier')
xlabel('hours relative to hurricane landfall')
plotFancyAxis

% ADV -------------------------------------------------------------------
% for the 0.003 lower bound
% iBlow = 1;%11;
% iBup = 316;%153; 
iBlow = 1;
iBup = 118;

iRRU = 3;
figure;
% subplot(3,1,1)
%     plot(dtLandADV(iBlow:iBup), stMTMWaveStat0pt003{1,iRRU}.Tm01Inf(iBlow:iBup), '-k')
% ylabel('mean period, Tm01_{inf} (s)')
% plotFancyAxis
% 
% subplot(3,1,2)
plot(dtLandADV(iBlow:iBup), stMTMWaveStat0pt003{1,iRRU}.infLow(iBlow:iBup), ...
     dtLandADV(iBlow:iBup), stMTMWaveStat0pt003{1,iRRU}.infMid(iBlow:iBup), ...
     dtLandADV(iBlow:iBup), stMTMWaveStat0pt003{1,iRRU}.infUp(iBlow:iBup))
ylabel('% S\eta\eta_{inf} (arb)')
legend('0.003 < f < 0.015', '0.015 < f < 0.028', '0.028 < f < 0.04')
plotFancyAxis

% subplot(3,1,3)
%     plot(dtLandADV(iBlow:iBup), movmean(stMTMWaveStat0pt0015{1,iRRU}.h(iBlow:iBup), 5), '-k')
% ylabel('mean water depth (m)')
% xlabel('hours relative to hurricane landfall')
% plotFancyAxis

%% normalized spectra and  comparison for RRUs & ADV (OSM & THESIS FIGS)
% updated on 6/2/18

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS
load('mtmWaveStats_17min_0pt0015_nw3k4.mat') 

% for RRUs only (ADV below)
iBlow = 49;%100; 
iBup  = 68;%120;

%nClow = 0.01416;%0.004;   
nCup  = 0.008;%0.008;
    
% plot normalized spectra RRUs
for iBurst = 65%iBlow : iBup
    figure; 
    loglog(stMTMWaveStat0pt0015{1,1}.f(:,iBurst), ...
           stMTMWaveStat0pt0015{1,1}.SnnInf(:,iBurst), ...%.*sqrt(nG*stMTMWaveStat0pt0015{1,1}.h(:,iBurst)), ...
           stMTMWaveStat0pt0015{1,2}.f(:,iBurst), ...
           stMTMWaveStat0pt0015{1,2}.SnnInf(:,iBurst))%.*sqrt(nG*stMTMWaveStat0pt0015{1,2}.h(:,iBurst)))
    %hold on
    %line([0.005 0.005], get(gca, 'ylim'), 'LineStyle', '--', 'Color', 'black')
    %hold on
    %line([nClow nClow], get(gca, 'ylim'), 'LineStyle', '--', 'Color', 'black')
    hold on
    xlabel('frequency (Hz)')
    ylabel('S\eta\eta_{inf} (arb)')
    set(gca, 'ylim', [0.00001 1], 'xlim', [0.001 0.1])
    line([nCup nCup], get(gca, 'ylim'), 'LineStyle', '--', 'Color', 'black')
    %legend('backshore', 'back-barrier')
    title(sprintf('%s, Burst %d', stMTMWaveStat0pt0015{1,2}.dtPb(1,iBurst), iBurst))
end

% plot normalized spectra ADV
iFmin = 4;  % corresponds to 0.0015 Hz
iFmax = 83; % corresponds to 0.04 Hz
for iBurst = 260 : 265
    figure;
    % I fucked up and didn't calculate the normalized confidence limits
    rM0inf  = trapz(stMTMWaveStat0pt0015{1,3}.f(iFmin:iFmax,iBurst), ...
                    stMTMWaveStat0pt0015{1,3}.SnnInf(iFmin:iFmax,iBurst));
    CIbelow = abs(stMTMWaveStat0pt0015{1,3}.CIlow(iFmin:iFmax,iBurst) ...
                  - stMTMWaveStat0pt0015{1,3}.Snn(iFmin:iFmax,iBurst)) / rM0inf;
    CIabove = ( stMTMWaveStat0pt0015{1,3}.CIup(iFmin:iFmax,iBurst) ...
             - stMTMWaveStat0pt0015{1,3}.Snn(iFmin:iFmax,iBurst) ) / rM0inf ;
    boundedline(stMTMWaveStat0pt0015{1,3}.f(iFmin:iFmax,iBurst)', ...
                stMTMWaveStat0pt0015{1,3}.SnnNormInf(iFmin:iFmax,iBurst)', ...
                [CIbelow CIabove]);
    xlabel('frequency (Hz)')
    ylabel('S\eta\eta_{inf} (arb)')
    set(gca, 'YScale', 'log',  'XScale', 'log', 'ylim', [1 1000], 'xlim', [0.003 0.04])
    plotFancyAxis
    title(sprintf('%s, Burst %d', stMTMWaveStat0pt0015{1,3}.dtPb(1,iBurst), iBurst))
    %       loglog(stMTMWaveStat0pt0015{1,3}.f(:,iBurst), ...
%            stMTMWaveStat0pt0015{1,3}.SnnNormInf(:,iBurst))
end

% normalized  using 0pt0015 spectra
% RRU1 --------------------------------------------------------------------
figure
subplot(1,2,1)
dtTemp = repmat(dtLandRRU1', size(stMTMWaveStat0pt0015{1,1}.f,1), 1);
surf(dtTemp(:,iBlow:iBup), ...
     stMTMWaveStat0pt0015{1,1}.f(:,iBlow : iBup), ...
     stMTMWaveStat0pt0015{1,1}.SnnNormInf(:,iBlow : iBup), ...
     'EdgeColor', 'none')
axis xy; view(2)
%datetick('x','mm/dd HH PM')
axis tight
% formatting
c = colorbar;
ylabel(c, 'S\eta\eta_{norm}(f) [arb]', ...
                'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
title('backshore')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0.0015 0.04]) 

% RRU2 --------------------------------------------------------------------
%figure
subplot(1,2,2)
dtTemp = repmat(dtLandRRU2', size(stMTMWaveStat0pt0015{1,2}.f,1), 1);
surf(dtTemp(:,iBlow:iBup), ...
     stMTMWaveStat0pt0015{1,2}.f(:,iBlow : iBup), stMTMWaveStat0pt0015{1,2}.SnnNormInf(:,iBlow : iBup), 'EdgeColor', 'none')
axis xy; view(2)
%datetick('x','mm/dd HH PM')
axis tight
% formatting
c = colorbar;
ylabel(c, 'S\eta\eta_{norm}(f) [arb]', ...
                'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
title('back-barrier')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0.0015 0.04]) 

% ADV ---------------------------------------------------------------------
iBlow = 1;
iBup = 316;
iFmin = 4;  % corresponds to 0.0015 Hz
iFmax = 83; % corresponds to 0.04 Hz

dtLandADV  = repmat(datetime('25-Aug-2017 22:00:00'), ...
size(stMTMWaveStat0pt0015{1,3}.dtPb(1,:),2), 1);
dtLandADV  = hours(stMTMWaveStat0pt0015{1,3}.dtPb(1,:)' - dtLandADV);

figure
%dtTemp = repmat(dtLandADV', size(stMTMWaveStat0pt0015{1,3}.f,1), 1);
% surf(dtTemp(:,iBlow:iBup), ...
%      stMTMWaveStat0pt0015{1,3}.f(:,iBlow:iBup), ...
%      stMTMWaveStat0pt0015{1,3}.SnnNormInf(:,iBlow:iBup), ...
%      'EdgeColor', 'none')
imagesc(dtLandADV(iBlow:iBup), stMTMWaveStat0pt0015{1,3}.f(iFmin:iFmax,1), ...
        stMTMWaveStat0pt0015{1,3}.SnnNormInf(iFmin:iFmax,iBlow:iBup)) 
set(gca,'YDir','normal')    
% axis xy; view(2)
%datetick('x','mm/dd HH PM')
% axis tight
% formatting
c = colorbar;
ylabel(c, 'S\eta\eta_{norm}(f) [arb]', ...
                'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
title('ADV')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0.003 0.04]) 

%% WSE + mean water level + precipitation + tides (OSM & THESIS FIGURE)
% modified to not plot ADCIRC outputs; replaced with tides

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM
load('RRU.mat') 
load('TIDE.mat')     % SLP predicted tides (m)
load('USGS.mat')
cmap = parula(5);
nMSL = 0.1;          % San Luis Pass MSL

% landfall datetime series for differencing (UTC)
dtLandRRU1 = repmat(datetime('26-Aug-2017 03:00:00'), ...
                    size(stPOD(1).PSdt,1), 1);
dtLandRRU2 = repmat(datetime('26-Aug-2017 03:00:00'), ...
                    size(stPOD(2).PSdt,1), 1);
dtLandRRU1 = hours(stPOD(1).PSdt+hours(5) - dtLandRRU1);
dtLandRRU2 = hours(stPOD(2).PSdt+hours(5) - dtLandRRU2);

% WSEs
figure
plot(stPOD(1).PSdt+hours(5), stPOD(1).wse, 'Color', cmap(1,:))
hold on
plot(stPOD(2).PSdt+hours(5), stPOD(2).wse, 'Color', cmap(3,:))
hold on
plot(stTide.dt+hours(5), stTide.pred, 'k', 'LineWidth', 2)
hold on
plot(stTide.dt+hours(5), repmat(nMSL, size(stTide.dt)))
ylabel('\eta (m NAVD88)')
legend('backshore RRU', 'back-barrier RRU', 'predicted tide', 'MSL')
plotFancyAxis

% WSEs
figure
plot(dtLandRRU1, stPOD(1).wse, 'Color', cmap(1,:))
hold on
plot(dtLandRRU2, stPOD(2).wse, 'Color', cmap(3,:))
ylabel('\eta (m NAVD88)')
legend('backshore RRU', 'back-barrier RRU', 'predicted tide', 'MSL')
plotFancyAxis

% precipitation
figure; 
plot(stUSGS.dtPrecip+hours(5), stUSGS.precip)
ylabel('precipitation (inches)')
plotFancyAxis()

%% pwelch comparison with MTM

% pwelch on burst 52
[rSppPW, rFpw] = pwelch(detrend(stMTMWaveStat0pt003{1,1}.Pb(:,53)), ...
                        [], [], [], 16);                  

figure
loglog(stMTMps{2,1}.f(2:42,53), stMTMps{2,1}.spp(2:42,53), ...
       stMTMps{1,1}.f(3:83,53), stMTMps{1,1}.spp(3:83,53), ...
       rFpw(2:12), rSppPW(2:12))
xlabel('frequency (Hz)')
ylabel('pressure power spectra, Spp (m^2/Hz)')
legend('Slepian multitapers: NW=3, K=4', 'Slepian multitaper: zero-padded', 'Hamming window: 50% overlap')
plotFancyAxis

%% Follets West WSE (OSM FIGURE & THESIS FIGURE)

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM
load('RRU.mat') 
load('VECT.mat')
cmap = parula(5);

% landfall datetime series for differencing (UTC)
dtLandRRU7 = repmat(datetime('26-Aug-2017 03:00:00'), ...
                    size(stPOD(4).PSdt,1), 1);
dtLandRRU6 = repmat(datetime('26-Aug-2017 03:00:00'), ...
                    size(stPOD(3).PSdt,1), 1);
dtLandADV  = repmat(datetime('26-Aug-2017 03:00:00'), ...
                    size(stVECT.dt,1), 1);
dtLandRRU7 = hours(stPOD(4).PSdt+hours(5) - dtLandRRU7);
dtLandRRU6 = hours(stPOD(3).PSdt+hours(5) - dtLandRRU6);
dtLandADV  = hours(stVECT.dt+hours(5) - dtLandADV);

% Follets West
figure
% convert from CST to UTC for plot
plot(stVECT.dt+hours(5), stVECT.wse, 'Color', cmap(1,:))
hold on
plot(stPOD(4).PSdt+hours(5), stPOD(4).wse, 'Color', cmap(4,:))
hold on
plot(stPOD(3).PSdt+hours(5), stPOD(3).wse, 'Color', cmap(3,:))
hold on
%line([datetime('25-Aug-2017 22:00:00') datetime('25-Aug-2017 22:00:00')], ...
%    get(gca, 'ylim'), 'LineStyle', '--', 'Color', 'black')
line([datetime('26-Aug-2017 03:00:00') datetime('26-Aug-2017 03:00:00')], ...
    get(gca, 'ylim'), 'LineStyle', '--', 'Color', 'black')
% formatting
legend('surf zone ADV', 'mid-barrier RRU', 'back-barrier RRU')
ylabel('\eta (m NAVD88)')
plotFancyAxis()

% now in hours relative to landfall
figure
% convert from CST to UTC for plot
plot(dtLandADV, stVECT.wse, 'Color', cmap(1,:))
hold on
plot(dtLandRRU7, stPOD(4).wse, 'Color', cmap(4,:))
hold on
plot(dtLandRRU6, stPOD(3).wse, 'Color', cmap(3,:))
hold on
% formatting
legend('surf zone ADV', 'mid-barrier RRU', 'back-barrier RRU')
ylabel('\eta (m NAVD88)')
plotFancyAxis()

%% sensitivity plots for ADV bispectra

load('IN-mtmBspd_ADVsp_17minZP.mat')
figure;
subplot(2,3,1)
    i=1:18; 
    semilogy(f(:,i), snn(:,i)*2);
    hold on
    semilogy(f(:,1), mean(snn(:,i)*2,2), 'k', 'LineWidth', 3); 
    xlim([0.001, 0.04])
    ylabel('$\sim$17 min Snn(f) (m$^2$/Hz)', 'Interpreter', 'LaTex', 'FontSize', 14)
subplot(2,3,2)
    i=19:82; 
    semilogy(f(:,i), snn(:,i)*2);
    hold on
    semilogy(f(:,1), mean(snn(:,i)*2,2), 'k', 'LineWidth', 3); 
    xlim([0.001, 0.04])
subplot(2,3,3)
    i=83:115; 
    semilogy(f(:,i), snn(:,i)*2);
    hold on
    semilogy(f(:,1), mean(snn(:,i)*2,2), 'k', 'LineWidth', 3); 
    xlim([0.001, 0.04])
    
load('IN-mtmBisp_ADVsp_8minZP.mat')
subplot(2,3,4)
    i=1:36; 
    semilogy(f(:,i), snn(:,i)*2);
    hold on
    semilogy(f(:,1), mean(snn(:,i)*2,2), 'k', 'LineWidth', 3);  
    xlim([0.001, 0.04])
    ylabel('$\sim$8.5 min Snn(f) (m$^2$/Hz)', 'Interpreter', 'LaTex', 'FontSize', 14)
subplot(2,3,5)
    i=37:163; 
    semilogy(f(:,i), snn(:,i)*2);
    hold on
    semilogy(f(:,1), mean(snn(:,i)*2,2), 'k', 'LineWidth', 3);  
    xlim([0.001, 0.04])
    xlabel('frequency (Hz)', 'Interpreter', 'LaTex', 'FontSize', 14)
subplot(2,3,6)
    i=164:234; 
    semilogy(f(:,i), snn(:,i)*2);
    hold on
    semilogy(f(:,1), mean(snn(:,i)*2,2), 'k', 'LineWidth', 3); 
    xlim([0.001, 0.04])  

%% Origin of IG waves (thesis figure); NDBC spectral width

    %%% NOTE: (2/2019) both the NDBC and ADCIRC datetimes were off when
    %%% this analysis was originally completed (ADCIRC by one hour, NDBC
    %%% by 6 hours [UTC]) - so much of this code is incorrect but keeping 
    %%% for reference)

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS
load('mtmWaveStats_17min_0pt003_nw3k4.mat') 
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/NDBC
load('NDBC.mat') 

dtLandADV  = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(stMTMWaveStat0pt003{1,3}.dtPb(1,:),2), 1);    
dtLandNDBC = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(stNDBC{1,1}.dt,2), 1);                
dtLandADV  = hours(stMTMWaveStat0pt003{1,3}.dtPb(1,:)' - dtLandADV);
dtLandNDBC = hours(stNDBC{1,1}.dt' - dtLandNDBC);

% ADV infragravity and offshore (NDBC) sea-swell 
iBlow = 1;
iBup = 316;

figure
subplot(3,1,1)
yyaxis left
plot(dtLandADV(iBlow:iBup), stMTMWaveStat0pt003{1,3}.Hm0Inf(iBlow:iBup))
%plot(dtLandADV(iBlow:iBup), movmean(stMTMWaveStat0pt003{1,3}.Hm0Inf(iBlow:iBup),4))
ylabel('Hm0_{inf} (m)')
hold on
yyaxis right
plot(dtLandADV(iBlow:iBup), stMTMWaveStat0pt003{1,3}.Tm01Inf(iBlow:iBup))
%plot(dtLandADV(iBlow:iBup), movmean(stMTMWaveStat0pt003{1,3}.Tm01Inf(iBlow:iBup),4))
ylabel('Tm01_{inf} (s)')
set(gca, 'xlim', [-23 66])
plotFancyAxis

% NDBC
iBlow = 573; % -23 hours
iBup  = 662; % + 66 hours

subplot(3,1,2)
yyaxis left
plot(dtLandNDBC(iBlow:iBup), stNDBC{1,1}.Hm0SS(iBlow:iBup))
ylabel('H_{0} (m)')
hold on
yyaxis right
plot(dtLandNDBC(iBlow:iBup), stNDBC{1,1}.Tm01SS(iBlow:iBup)) 
ylabel('T_{0} (s)')
set(gca, 'xlim', [-23 66])
plotFancyAxis

subplot(3,1,3)
yyaxis left
plot(dtLandNDBC(iBlow:iBup), stNDBC{1,1}.sw2(iBlow:iBup))
ylabel('v')
hold on
yyaxis right
plot(dtLandNDBC(iBlow:iBup), stNDBC{1,1}.qp(iBlow:iBup)) 
ylabel('Q_{p}')
xlabel('hours relative to hurricane landfall')
set(gca, 'xlim', [-23 66])
plotFancyAxis

%% spectrograms of u,v,n - THESIS FIGURE (updated 4/2/19)

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/VELOCITY
load('IN-mtmCrossCU_ADV_68min.mat') % NW=3, K=4
rSuu = suu * 2; % make single sided
rSvv = svv * 2; % make single sided
iFmin = 1;
iFmax = 165;

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV
load('OUT-mtmADVsp_nPLUS_68min.mat', 'dt', 'snn') % NW=3, K=4
dtLandADV = repmat(datetime('25-Aug-2017 22:00:00', 'TimeZone', 'America/Chicago'), size(dt,2), 1);                
dtLandADV = hours(dt(1,:)' - dtLandADV)';

figure
subplot(1,3,1)
imagesc(dtLandADV, f(iFmin:iFmax,1), rSuu(iFmin:iFmax,:)) 
set(gca,'YDir','normal')    
% formatting
c = colorbar;
ylabel(c, 'S_{uu}(f) [(m/s)^2/Hz]', 'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
title('Suu')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0 0.04]) 

subplot(1,3,2)
imagesc(dtLandADV, f(iFmin:iFmax,1), rSvv(iFmin:iFmax,:)) 
set(gca,'YDir','normal')    
% formatting
c = colorbar;
ylabel(c, 'S_{vv}(f) [(m/s)^2/Hz]', 'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
title('Svv')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0 0.04]) 

subplot(1,3,3)
imagesc(dtLandADV, f(iFmin:iFmax,1), snn(iFmin:iFmax,:))
set(gca,'YDir','normal')    
% formatting
c = colorbar;
ylabel(c, 'S_{\eta\eta}+ (f) [m^2/Hz]', 'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
title('Snn')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0 0.04]) 

% 17 minute ---------------------------------------------------------------
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/VELOCITY
load('IN-mtmCrossCU_ADV_17min.mat') % NW=3, K=4
rSuu = suu * 2; % make single sided
rSvv = svv * 2; % make single sided
iFmin = 4;
iFmax = 42;

cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVE_STATS
load('mtmWaveStats_17min_0pt003_FINAL.mat') % NW=3, K=4
dt = stMTMWaveStat0pt003{2,3}.dtPb;

dtLandADV = repmat(datetime('25-Aug-2017 22:00:00'), size(dt,2), 1);                
dtLandADV = hours(dt(1,:)' - dtLandADV)';

% try different colormaps
%cm=flipud(magma(100));
%cm=inferno(100);
%cm=flipud(plasma(100));
cm=flipud(viridis(100));

figure
h1 = subplot(1,3,1);
imagesc(dtLandADV, f(iFmin:iFmax,1), rSuu(iFmin:iFmax,:)) 
set(gca,'YDir','normal')    
% formatting
colormap(h1,cm) 
c = colorbar;
ylabel(c, 'S_{uu}(f) [(m/s)^2/Hz]', 'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
title('Suu')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0.003 0.04]) 

subplot(1,3,2)
imagesc(dtLandADV, f(iFmin:iFmax,1), rSvv(iFmin:iFmax,:)) 
set(gca,'YDir','normal')    
% formatting
c = colorbar;
ylabel(c, 'S_{vv}(f) [(m/s)^2/Hz]', 'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
title('Svv')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0.003 0.04]) 

subplot(1,3,3)
imagesc(dtLandADV, stMTMWaveStat0pt003{2,3}.f(:,1), stMTMWaveStat0pt003{2,3}.SnnInf) 
set(gca,'YDir','normal')    
% formatting
c = colorbar;
ylabel(c, 'S_{\eta\eta}(f) [m^2/Hz]', 'FontSize', 13);
ylabel('frequency (Hz)');
xlabel('hours relative to hurricane landfall')
title('Snn')
% make it fancy
plotFancyAxis()
set(gca, 'YGrid', 'off', 'YLim', [0.003 0.04]) 

% check when or if V is > 75% of U
rM0u = trapz(f(iFmin:iFmax,1), rSuu(iFmin:iFmax,:));
rM0v = trapz(f(iFmin:iFmax,1), rSvv(iFmin:iFmax,:));
figure; 
plot(dtLandADV, rM0v ./ rM0u)   

%% low-pass filter u, v, and n (infragravity band and below, 17min)

nScrsz = get(groot, 'ScreenSize');  

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/VELOCITY
load('OUT-mtmADVvelU_17min_ZP.mat', 'data', 'zmts', 'dt') % NW=3, K=4
cU = struct2cell(data);
rUzm = zmts;

load('OUT-mtmADVvelV_17min_ZP.mat', 'data', 'zmts') % NW=3, K=4
cV = struct2cell(data);
rVzm = zmts;

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/ADV/PowerSpectra/MTM_R
load('mtmADVsp_17min_ZP.mat', 'zmts')
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/
load('mtmWaveStats_17min_0pt0015_nw3k4.mat') % NW=3, K=4
% I can't find a windowed 17min WSE so here is use the sea-pressure + the
% PT elevation above the bed and then convert to NAVD88 based on the bed
% level elevation
nVECTele   = -0.746;                     % bed elevation at ADV
nVECTPTele = nVECTele + 0.275 + 0.425;   % ADV (PT sensor 0.7 m above bed)
nRho       = 1025;   % seawater density [kg/m3]
nG         = 9.81;   % gravitational constant [m/s2]
       % in Pascals
rWSE = ( stMTMWaveStat0pt0015{1,3}.Pb / (nRho * nG) ) + nVECTPTele;
rWSEzm = zmts;

% this is the same filter used in szWaveStat_Harvey
    % finite impulse response filters (usally constant filter delay - 
    % linear phase filters)
nFminSS = 0.04;
nFsamp  = 16;
dHamLPF = designfilt('lowpassfir', 'FilterOrder', 200, ...
                     'CutoffFrequency', nFminSS, 'SampleRate', nFsamp);
% plot filter (for debugging)
fvtool(dHamLPF)                     

% plot LPF
for iBurst = 66 %:75 %: size(f,2)
    rUlpf   = filtfilt(dHamLPF, cU{iBurst});
    rVlpf   = filtfilt(dHamLPF, cV{iBurst});
    rWSElpf = filtfilt(dHamLPF, rWSE(:,iBurst));
    
    rUlpfzm   = filtfilt(dHamLPF, rUzm(:,iBurst));
    rVlpfzm   = filtfilt(dHamLPF, rVzm(:,iBurst));
    rWSElpfzm = filtfilt(dHamLPF, rWSEzm(:,iBurst));

%     figure;
%     yyaxis left
%     plot(dt(:,iBurst), rUlpf, dt(:,iBurst), rVlpf) 
%     ylabel('velocity_{LPF} (m/s)')
%     yyaxis right
%     plot(dt(:,iBurst), rWSElpf)
%     ylabel('\eta_{LPF} (m NAVD88)')
%     legend('u', 'v', '\eta') 
%     plotFancyAxis
%     
    hFig = figure;
    yyaxis left
    plot(dt(:,iBurst), rUlpfzm, dt(:,iBurst), rVlpfzm) 
    ylabel('velocity_{LPF} (m/s)')
    yyaxis right
    plot(dt(:,iBurst), rWSElpfzm)
    ylabel('\eta_{LPF} (m)')
    legend('u', 'v', '\eta') 
    plotFancyAxis
    
    set(hFig, 'Position',[1 nScrsz(4)*2/4 nScrsz(3)*3/4 nScrsz(4)*2/4]);
end

%% RRU1 and RRU2 wave spectra examples w/ CL (thesis figure)

% NOTE: the confidence limits in stMTMWaveStat0pt003 and 0pt0015 look wrong 
% at low frequencies; re-calculated using ThomsonMTM_NWsensitivity.R which
% doesn't modify the spectrum for burial etc.; shouldn't matter for
% infragravity frequencies

%cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS
%load('mtmWaveStats_17min_0pt0015_nw3k4.mat')

% NOTE: boundedline doesn't work for the full confidence interval time
% series...not sure why
iBurst = 55;
iFmax  = 83; % corresponds to 0.04 Hz

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/RRUs/PowerSpectra/MTM_R
load('mtmRRU1_17min_ZP.mat')
% RRU-1
% CIbelow = abs(stMTMWaveStat0pt0015{1,1}.CIlow(1:iFmax,iBurst) ...
%               - stMTMWaveStat0pt0015{1,1}.SnnInf(1:iFmax,iBurst));
% CIabove = stMTMWaveStat0pt0015{1,1}.CIup(1:iFmax,iBurst) - ...
%           stMTMWaveStat0pt0015{1,1}.SnnInf(1:iFmax,iBurst);
% boundedline(stMTMWaveStat0pt0015{1,1}.f(1:iFmax,iBurst)', ...
%             stMTMWaveStat0pt0015{1,1}.SnnInf(1:iFmax,iBurst)', ...
%             [CIbelow CIabove]);
figure
subplot(1,2,1)
CIbelow = abs(CIlow(1:iFmax,iBurst) - snn(1:iFmax,iBurst));
CIabove = CIup(1:iFmax,iBurst) - snn(1:iFmax,iBurst);
[l,p]=boundedline(f(1:iFmax,iBurst)', snn(1:iFmax,iBurst)', ...
            [CIbelow CIabove]); 
outlinebounds(l,p);  
hold on
load('mtmRRU2_17min_ZP.mat')
% RRU-2
% CIbelow = abs(stMTMWaveStat0pt0015{1,2}.CIlow(1:iFmax,iBurst) ...
%               - stMTMWaveStat0pt0015{1,2}.SnnInf(1:iFmax,iBurst));
% CIabove = stMTMWaveStat0pt0015{1,2}.CIup(1:iFmax,iBurst) - ...
%           stMTMWaveStat0pt0015{1,2}.SnnInf(1:iFmax,iBurst);
% boundedline(stMTMWaveStat0pt0015{1,2}.f(1:iFmax,iBurst)', ...
%             stMTMWaveStat0pt0015{1,2}.SnnInf(1:iFmax,iBurst)', ...
%             [CIbelow CIabove]);
CIbelow = abs(CIlow(1:iFmax,iBurst) - snn(1:iFmax,iBurst));
CIabove = CIup(1:iFmax,iBurst) - snn(1:iFmax,iBurst);
[l,p]=boundedline(f(1:iFmax,iBurst)', snn(1:iFmax,iBurst)', ...
            [CIbelow CIabove]); 
outlinebounds(l,p);        
xlim([0.001 0.04])
ylim([0.00001 2])
set(gca, 'YScale', 'log', 'Xscale', 'log')
plotFancyAxis
axis square
ylabel('S\eta\eta_{inf} [m^2/Hz]', 'FontSize', 18)
   
iBurst = 65;
iFmax  = 83; % corresponds to 0.04 Hz

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/RRUs/PowerSpectra/MTM_R
load('mtmRRU1_17min_ZP.mat')
% RRU-1
% CIbelow = abs(stMTMWaveStat0pt0015{1,1}.CIlow(1:iFmax,iBurst) ...
%               - stMTMWaveStat0pt0015{1,1}.SnnInf(1:iFmax,iBurst));
% CIabove = stMTMWaveStat0pt0015{1,1}.CIup(1:iFmax,iBurst) - ...
%           stMTMWaveStat0pt0015{1,1}.SnnInf(1:iFmax,iBurst);
% boundedline(stMTMWaveStat0pt0015{1,1}.f(1:iFmax,iBurst)', ...
%             stMTMWaveStat0pt0015{1,1}.SnnInf(1:iFmax,iBurst)', ...
%             [CIbelow CIabove]);
subplot(1,2,2)
CIbelow = abs(CIlow(1:iFmax,iBurst) - snn(1:iFmax,iBurst));
CIabove = CIup(1:iFmax,iBurst) - snn(1:iFmax,iBurst);
[l,p]=boundedline(f(1:iFmax,iBurst)', snn(1:iFmax,iBurst)', ...
            [CIbelow CIabove]); 
outlinebounds(l,p); 
hold on
load('mtmRRU2_17min_ZP.mat')
% RRU-2
% CIbelow = abs(stMTMWaveStat0pt0015{1,2}.CIlow(1:iFmax,iBurst) ...
%               - stMTMWaveStat0pt0015{1,2}.SnnInf(1:iFmax,iBurst));
% CIabove = stMTMWaveStat0pt0015{1,2}.CIup(1:iFmax,iBurst) - ...
%           stMTMWaveStat0pt0015{1,2}.SnnInf(1:iFmax,iBurst);
% boundedline(stMTMWaveStat0pt0015{1,2}.f(1:iFmax,iBurst)', ...
%             stMTMWaveStat0pt0015{1,2}.SnnInf(1:iFmax,iBurst)', ...
%             [CIbelow CIabove]);
CIbelow = abs(CIlow(1:iFmax,iBurst) - snn(1:iFmax,iBurst));
CIabove = CIup(1:iFmax,iBurst) - snn(1:iFmax,iBurst);
[l,p]=boundedline(f(1:iFmax,iBurst)', snn(1:iFmax,iBurst)', ...
            [CIbelow CIabove]); 
outlinebounds(l,p);        
xlim([0.001 0.04])
ylim([0.00001 2])
set(gca, 'YScale', 'log', 'Xscale', 'log')
plotFancyAxis
xlabel('frequency (Hz)')
axis square

%% wave rose for TCM bearing throughout the storm

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/
load('RRU.mat')

% break up storm into segments
% NOTE: at '25-Aug-2017 23:15:17' to '26-Aug-2017 00:15:17' the TCM appears
% buried; I'm only really confident through 19:00
sSTime1  = {'25-Aug-2017 10:15:17', '25-Aug-2017 13:15:17', ...
           '25-Aug-2017 18:00:00'};
%           '25-Aug-2017 16:15:17'};
sETime1  = {'25-Aug-2017 13:15:17', '25-Aug-2017 16:15:17', ... 
           '25-Aug-2017 19:00:00'};
%           '25-Aug-2017 19:15:17'};
       
for i=3 %: length(sSTime1) 
    
    % find index
    iMinTCM = find(stPOD(1).TCMdt == sSTime1{i});
    iMaxTCM = find(stPOD(1).TCMdt == sETime1{i});
    
    % create dummy variable
    rI = ones(size(stPOD(1).bearing(iMinTCM:iMaxTCM)));

    % plot
    figure;
    %hist(stPOD(1).bearing(iMinTCM:iMaxTCM))
    plot(stPOD(1).TCMdt(iMinTCM:iMaxTCM),stPOD(1).bearing(iMinTCM:iMaxTCM))
    figure
    %wind_rose(stPOD(1).bearing(iMinTCM:iMaxTCM), rI, ...
    %              'dtype', 'meteo')
    wind_rose(stPOD(1).bearing(iMinTCM:iMaxTCM), rI)
end

%% thesis figure: Hm0 in 3 infragravity bands at FI (17 min) [saturation]

% calculate Hm0 using the zero-padded spectra to better isolate peaks
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/ADV/PowerSpectra/MTM_R
load('mtmADVsp_17min_ZP.mat')
iLFlow   = [7 26 57];   % low cutoffs (0.0029, 0.0122, 0.0273) 
iLFhigh  = [25 83 83];  % high cutoffs (0.0117, 0.0269, 0.04) 

% individual bands
for i = 1 : length(iLFlow)
    rM0snnLF = trapz(f(iLFlow(i):iLFhigh(i),1), snn(iLFlow(i):iLFhigh(i),:));
    rHm0LF(i,:) = 4.004 * sqrt(rM0snnLF); %#ok<SAGROW> % each row corresponds to a band
end

% all of IG band
rM0snnLFt = trapz(f(7:83,1), snn(7:83,:));
rHm0LFt = 4.004 * sqrt(rM0snnLFt);

% try cumtrapz
rM0snnLFctz = cumtrapz(f(7:83,1), snn(7:83,:));
rHm0LFctz(1,:) = 4.004 * sqrt(rM0snnLFctz(19,:));
rHm0LFctz(2,:) = rHm0LFctz(1,:)- (4.004 * sqrt(rM0snnLFctz(50,:)));
rHm0LFctz(3,:) = rHm0LFctz(1,:) - rHm0LFctz(2,:) - (4.004 * sqrt(rM0snnLFctz(end,:)));
%rHm0LFctz(3,:) = 4.004 * sqrt(rM0snnLFctz(end,:));

% plot IG wave heights over each band
figure;
plot(dtLandADV', rHm0LF(1,:), dtLandADV', rHm0LF(2,:), ...
     dtLandADV', rHm0LF(3,:), 'LineWidth', 2)
ylabel('Hm0_{IG} (m/s)')
plotFancyAxis

figure;
plot(dtLandADV', rHm0LF(1,:)+rHm0LF(2,:)+rHm0LF(3,:), dtLandADV', rHm0LFt, 'LineWidth', 2)
ylabel('Hm0_{IG} (m/s)')
plotFancyAxis

figure;
plot(dtLandADV', rHm0LFctz(1,:)+rHm0LFctz(2,:)+rHm0LFctz(3,:), dtLandADV', rHm0LFt, 'LineWidth', 2)
ylabel('Hm0_{IG} (m/s)')
plotFancyAxis

%% thesis figure: percent variance in SIG, IG, SEA, SWELL (68 min)

%%% integrates the normalized full snn spectra and SIG-IG spectra only -- 
%%% both have variance of unity -- in the SIG, IG, and sea-swell and for 
%%% ~68 min wave spectra zero-padded spectra

nClow = [0.0004 0.003 0.04];  % SIG, IG, SS
nCup  = [0.003  0.04  0.6];  

% surf zone
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/ADV/Bispectra/ADVsp-68min
load('IN-mtmBisp_ADVsp_68minZP.mat', 'snn', 'f')
load('IN-mtmBisp_ADVsp_68min.mat', 'dtLandADV')

% normalize full snn spectra (note that this snn has not been modified for
% attenuation of any signals as it should be at SS frequencies - will need
% to modify in the future)
rM0 = trapz(f(:,1), snn);
rSnnNorm = snn ./ repmat(rM0, length(snn), 1);
% trapz(f(:,1), rSnnNorm) % test that variance is indeed 1
    
for iBand = 1 : 3
    
    % find indices for each band
    [~, iFlow(iBand)] = min(abs(f(:,1) - nClow(iBand))); %#ok<SAGROW>
    [~, iFup(iBand)]  = min(abs(f(:,1) - nCup(iBand))); %#ok<SAGROW>
    
end

% SIG
rM0SIG = trapz(f(iFlow(1):iFup(1),1), rSnnNorm(iFlow(1):iFup(1),:));
mean(rM0SIG)

% IG
rM0IG  = trapz(f(iFlow(2)+1:iFup(2),1), rSnnNorm(iFlow(2)+1:iFup(2),:));
mean(rM0IG)

% SS
rM0SS  = trapz(f(iFlow(3)+1:iFup(3),1), rSnnNorm(iFlow(3)+1:iFup(3),:));
mean(rM0SS)

% plot
figure; 
plot(dtLandADV, rM0SIG, dtLandADV, rM0IG, dtLandADV, rM0SS, 'LineWidth', 2)
plotFancyAxis

% back-barrier
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/SensitivityTests/LowFreqSignals
load('mtmRRU6_68min_ZP.mat', 'snn', 'f')

% normalize full snn spectra (note that this snn has not been modified for
% attenuation of any signals as it should be at SS frequencies - will need
% to modify in the future)
rM0 = trapz(f(:,1), snn);
rSnnNorm = snn ./ repmat(rM0, length(snn), 1);
% trapz(f(:,1), rSnnNorm) % test that variance is indeed 1

% SIG
rM0SIG = trapz(f(iFlow(1):iFup(1),1), rSnnNorm(iFlow(1):iFup(1),:));
mean(rM0SIG)

% IG
rM0IG  = trapz(f(iFlow(2)+1:iFup(2),1), rSnnNorm(iFlow(2)+1:iFup(2),:));
mean(rM0IG)

% SS
rM0SS  = trapz(f(iFlow(3)+1:iFup(3),1), rSnnNorm(iFlow(3)+1:iFup(3),:));
mean(rM0SS)

% plot
figure; 
plot(rM0SIG)
plotFancyAxis

% backshore
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/SensitivityTests/LowFreqSignals
load('mtmRRU1_68min.mat', 'snn', 'f')

% normalize full snn spectra (note that this snn has not been modified for
% attenuation of any signals as it should be at SS frequencies - will need
% to modify in the future)
rM0 = trapz(f(:,1), snn);
rSnnNorm = snn ./ repmat(rM0, length(snn), 1);
% trapz(f(:,1), rSnnNorm) % test that variance is indeed 1

% SIG
rM0SIG = trapz(f(iFlow(1):iFup(1),1), rSnnNorm(iFlow(1):iFup(1),:));
mean(rM0SIG)

% IG
rM0IG  = trapz(f(iFlow(2)+1:iFup(2),1), rSnnNorm(iFlow(2)+1:iFup(2),:));
mean(rM0IG)

% SS
rM0SS  = trapz(f(iFlow(3)+1:iFup(3),1), rSnnNorm(iFlow(3)+1:iFup(3),:));
mean(rM0SS)

% plot
figure; 
plot(rM0SIG); hold on; plot(rM0IG) 
plotFancyAxis

%% biphase in the surf zone
% evaluate the biphase for 4 case scenarios

% preamble 
%cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/ADV/Bispectra/ADVsp-8min/4hr/OUT
cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/ADV/Bispectra/ADVsp-8min/4hr/OUT/
stFiles = dir;
stFiles = stFiles(~[stFiles.isdir]); % keep files, not directories
nNW  = 3;  % multitaper input parameters
nM   = 28; % ~8 minute segments
rBst = 1 : 28 : 589; % for plotting
[rBphc1, rBphc2l, rBphc2m, rBphc2u, rBphc3, idt] = deal(zeros(numel(stFiles),1));

% CASE 1: full IG-SWELL division ------------------------------------------
rFbndC1 = [0.04 0.14 0 0.04];

% CASE 1: subdivide IG-SWELL ----------------------------------------------
rFbndC2 = [0.04 0.14 0.000 0.011; ...   % low band
           0.04 0.14 0.012 0.026; ... % middle band
           0.04 0.14 0.027 0.040];     % high band    

% CASE 3: swell-swell division --------------------------------------------
rFbndC3 = [0.04 0.14 0.04 0.14];
                    
for i = 6%1 : numel(stFiles)            
    load(stFiles(i).name,'zBspec', 'rBic1', 'rF');
    
    % calculate biphase for each case
    rBphc1(i)  = mtmBiphase(zBspec, rBic1, rF, rFbndC1, nNW, nM, 'true');
    %rBphc2l(i) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(1,:), nNW, nM, 'false');
    %rBphc2m(i) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(2,:), nNW, nM, 'false');
    %rBphc2u(i) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(3,:), nNW, nM, 'false');
    %rBphc3(i)  = mtmBiphase(zBspec, rBic1, rF, rFbndC3, nNW, nM, 'false');
    
    % find the corresponding index for getting the datetime (average of the
    % interval)
    idt(i) = rBst(i)+13;
    
end    

% get datetime vector
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/ADV/Bispectra/ADVsp-8min/
load('IN-mtmBisp_ADVsp_8minZP', 'dtLandADV')

% plot
figure; 
subplot(1,3,1)
scatter(dtLandADV(idt),rBphc1) 
xlabel('hours relative to hurricane landfall')
ylabel('\beta_{IG,SW}-full (^o)')
plotFancyAxis
axis square

subplot(1,3,2)
scatter(dtLandADV(idt),rBphc2l); hold on;
scatter(dtLandADV(idt),rBphc2m); hold on;
scatter(dtLandADV(idt),rBphc2u); 
xlabel('hours relative to hurricane landfall')
ylabel('\beta_{IG,SW}-sub (^o)')
legend('lower', 'middle', 'upper')
plotFancyAxis
axis square

subplot(1,3,3)
scatter(dtLandADV(idt),rBphc3) 
ylabel('\beta_{SW,SW}-full (^o)')
xlabel('hours relative to hurricane landfall')
plotFancyAxis
axis square
    
% % period D (-11 to -7 hrs)
% load('OUT-mtmBisp_ADVsp_8minB85-112.mat', 'zBspec', 'rBic1', 'rF')
% rBphc1(1)  = mtmBiphase(zBspec, rBic1, rF, rFbndC1, nNW, nM);
% rBphc2l(1) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(1,:), nNW, nM);
% rBphc2m(1) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(2,:), nNW, nM);
% rBphc2u(1) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(3,:), nNW, nM);
% rBphc3(1)  = mtmBiphase(zBspec, rBic1, rF, rFbndC3, nNW, nM);
% 
% % period E (-3 to 1 hrs)
% load('OUT-mtmBisp_ADVsp_8minB141-168.mat', 'zBspec', 'rBic1', 'rF')
% rBphc1(2)  = mtmBiphase(zBspec, rBic1, rF, rFbndC1, nNW, nM);
% rBphc2l(2) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(1,:), nNW, nM);
% rBphc2m(2) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(2,:), nNW, nM);
% rBphc2u(2) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(3,:), nNW, nM);
% rBphc3(2)  = mtmBiphase(zBspec, rBic1, rF, rFbndC3, nNW, nM);
% 
% % period F (50 to 54 hrs)
% load('OUT-mtmBisp_ADVsp_8minB515-543.mat', 'zBspec', 'rBic1', 'rF')
% rBphc1(3)  = mtmBiphase(zBspec, rBic1, rF, rFbndC1, nNW, nM);
% rBphc2l(3) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(1,:), nNW, nM);
% rBphc2m(3) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(2,:), nNW, nM);
% rBphc2u(3) = mtmBiphase(zBspec, rBic1, rF, rFbndC2(3,:), nNW, nM);
% rBphc3(3)  = mtmBiphase(zBspec, rBic1, rF, rFbndC3, nNW, nM);

% figure; 
% subplot(1,3,1)
% scatter([-11 -3 50],rBphc1) 
% xlabel('hours relative to hurricane landfall')
% ylabel('\beta_{IG,SW}-full (^o)')
% plotFancyAxis
% axis square
% 
% subplot(1,3,2)
% scatter([-11 -3 50],rBphc2l); hold on;
% scatter([-11 -3 50],rBphc2m); hold on;
% scatter([-11 -3 50],rBphc2u); 
% xlabel('hours relative to hurricane landfall')
% ylabel('\beta_{IG,SW}-sub (^o)')
% legend('lower', 'middle', 'upper')
% plotFancyAxis
% axis square
% 
% subplot(1,3,3)
% scatter([-11 -3 50],rBphc3) 
% ylabel('\beta_{SW,SW}-full (^o)')
% xlabel('hours relative to hurricane landfall')
% plotFancyAxis
% axis square

%% variance reduction between the backshore and back barrier

cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVE_STATS
%load('mtmWaveStats_17min_0pt003_FINAL.mat') 
load('mtmWaveStats_17min_0pt0015_nw3k4.mat')

% for RRUs only 
iBlow = 49;
iBup  = 68;

iClow = 7;  % 0.0029 Hz
%iCup  = 18; % 0.0083 Hz

%iClow = 19; % 0.0088 Hz
iCup  = 83; % 0.04 Hz

% window backshore data
rF       = stMTMWaveStat0pt0015{1,1}.f(iClow:iCup,1);
rSnnRRU1 = stMTMWaveStat0pt0015{1,1}.SnnInf(iClow:iCup,iBlow:iBup);
rM0RRU1  = trapz(rF, rSnnRRU1);

% window back barrier data
rF       = stMTMWaveStat0pt0015{1,2}.f(iClow:iCup,1);
rSnnRRU2 = stMTMWaveStat0pt0015{1,2}.SnnInf(iClow:iCup,iBlow:iBup);
rM0RRU2  = trapz(rF, rSnnRRU2);

% variance reduction
rLoss = (rM0RRU1 - rM0RRU2) ./ rM0RRU1;

mean(rLoss)
std(rLoss)

%% values for table of wave stats in IG energy paper

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS
load('mtmWaveStats_17min_0pt003_nw3k4.mat') 

% RRU-1
maxHinfRRU1  = max(stMTMWaveStat0pt0015{1,1}.Hm0Inf); 
meanHinfRRU1 = mean(stMTMWaveStat0pt0015{1,1}.Hm0Inf);
stdHinfRRU1  = std(stMTMWaveStat0pt0015{1,1}.Hm0Inf);

maxHssRRU1  = max(stMTMWaveStat0pt0015{1,1}.Hm0SS); 
meanHssRRU1 = mean(stMTMWaveStat0pt0015{1,1}.Hm0SS);
stdHssRRU1  = std(stMTMWaveStat0pt0015{1,1}.Hm0SS);

% RRU-2
maxHinfRRU2  = max(stMTMWaveStat0pt0015{1,2}.Hm0Inf); 
meanHinfRRU2 = mean(stMTMWaveStat0pt0015{1,2}.Hm0Inf);
stdHinfRRU2  = std(stMTMWaveStat0pt0015{1,2}.Hm0Inf);

maxHssRRU2  = max(stMTMWaveStat0pt0015{1,2}.Hm0SS);
meanHssRRU2 = mean(stMTMWaveStat0pt0015{1,2}.Hm0SS);
stdHssRRU2  = std(stMTMWaveStat0pt0015{1,2}.Hm0SS);

% ADV
maxHinfADV  = max(stMTMWaveStat0pt003{1,3}.Hm0Inf); 
meanHinfADV = mean(stMTMWaveStat0pt003{1,3}.Hm0Inf);
stdHinfADV  = std(stMTMWaveStat0pt003{1,3}.Hm0Inf);

maxHssADV  = max(stMTMWaveStat0pt003{1,3}.Hm0SS); 
meanHssADV = mean(stMTMWaveStat0pt003{1,3}.Hm0SS);
stdHssADV  = std(stMTMWaveStat0pt003{1,3}.Hm0SS);

% ADV incoming IG wave heights (bootstrap)
% MTM: slepian tapers (using inputs from package multitaper in R)
% (NW=3, K=4)
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/ADV/Bispectra/ADVsp-17min
load('IN-mtmBisp_ADVnPLUS_17minZP.mat', 'snn', 'f')
iIGlow   = 7;
iIGhigh  = 83;
% variance over the IG band (double sided spectra)
rM0inf   = trapz(f(iIGlow:iIGhigh,1), snn(iIGlow:iIGhigh,:)*2);  
rNrmsInf = sqrt(rM0inf);           % std deviation (sigma)
rHm0Inf  = 4.004 * rNrmsInf;       % 0-moment wave height (H1/3)
max(rHm0Inf)
mean(rHm0Inf)
std(rHm0Inf)

%% infragravity reflectance (updated version in paper plotting code)

% change to the data directory on external hard-drive 
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV/REFLECTION
load('IN-mtmCrossReflect_ADVn_17minZP')

% input parameters
nBeta   = 0.0174;     % local bed slope
rTpIG   = [1/0.008 1/0.020 1/0.034]; % peak IG wave heights
iFmin   = [7 25 57];  % 7 = 0.0029 Hz
iFmax   = [24 56 83]; % corresponds to 0.04 Hz

% load('IN-mtmCrossReflect_ADVn_8minZP')
% iFmin   = [4 14 30];  % 7 = 0.0029 Hz
% iFmax   = [13 29 42]; % corresponds to 0.04 Hz

% power-spectra (double-sided)
% NOTE: I checked using the non-adaptive snn that the following code is correct
% for creating the mean power spectra directly from the eigencoefficients
rSuu = mean(abs(Yku).^2, 2);
rSnn = mean(abs(Ykn).^2, 2);

% cross-spectra
zSnu = mean(Ykn .* conj(Yku), 2);
rSnu = real(zSnu);   % co-spectrum

% msc
% NOTE: I checked against this output against the multitaper msc -->
% correct!
zMscNU = ( abs(zSnu).^2 ) ./ ( rSnn .* rSuu );

% reshape arrays so they are 2d and multiply by 2 to make them one-sided
zMscNU = reshape(zMscNU, [length(zMscNU) size(zMscNU,3)]);
rSuu = reshape(rSuu, [length(rSuu) size(rSuu,3)]) * 2;
rSnn = reshape(rSnn, [length(rSnn) size(rSnn,3)]) * 2;
rSnu = reshape(rSnu, [length(rSnu) size(rSnu,3)]) * 2;

% energy and energy flux density (from Sheremet et al. [2002])
nG    = 9.81; % [m/s2]
rH    = repmat(h, size(rSnn,1), 1);
rEin  = 1/4 * ( rSnn + ( (rH/nG) .* rSuu) + ( (2 * sqrt(rH/nG)) .* rSnu) ); % m^2 s + m^2/s^2 (m^4/s^4) * s
rEout = 1/4 * ( rSnn + ( (rH/nG) .* rSuu) - ( (2 * sqrt(rH/nG)) .* rSnu) );
rFin  = rEin .* sqrt(nG*rH); 
rFout = rEout .* sqrt(nG*rH);

% load the n+ generated spectra for comparison
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/ADV/Bispectra/ADVsp-17min
load('IN-mtmBisp_ADVnPLUS_17minZP.mat', 'snn', 'f')

% now integrate over the infragravity band to get the energy flux and
% calculate the reflection coefficient
[rFinIG, rFoutIG, rR2, rR, rHigIN, rHm0Inf, rBetaH, rHigOUT] = deal(zeros(3,size(rFin,2)));
for i = 1 : length(iFmin)
    rFinIG(i,:)  = trapz(f(iFmin(i):iFmax(i),1), rFin(iFmin(i):iFmax(i),:));
    rFoutIG(i,:) = trapz(f(iFmin(i):iFmax(i),1), rFout(iFmin(i):iFmax(i),:));
    rR2(i,:)     = rFoutIG(i,:) ./ rFinIG(i,:);
    rR(i,:)      = sqrt(rR2(i,:));
    
    % calculate incoming Hig
    rHigIN(i,:)  = 4 * sqrt(trapz(f(iFmin(i):iFmax(i),1), rEin(iFmin(i):iFmax(i),:)));
    
    % calculate significant infragravity wave heights with this method
    rHigOUT(i,:) = 4 * sqrt(trapz(f(iFmin(i):iFmax(i),1), rEout(iFmin(i):iFmax(i),:)));
    
    % calculate HIG+ using n+
    rM0inf   = trapz(f(iFmin(i):iFmax(i),1), snn(iFmin(i):iFmax(i),:)*2);  
    rNrmsInf = sqrt(rM0inf);                % std deviation (sigma)
    rHm0Inf(i,:)  = 4.004 * rNrmsInf;       % 0-moment wave height (H1/3)

    % normalized bed slope parameter
    rBetaH(i,:) = ( (nBeta*rTpIG(i)) / (2*pi) ) * sqrt(nG ./ rHigIN(i,:));
    
end

% plot three different IG bands for reflectance
figure
% subplot(3,1,1)
% scatter(dtLandADV, rR2(1,:)); hold on; 
% scatter(dtLandADV, rR2(2,:)); hold on;
% scatter(dtLandADV, rR2(3,:));
% ylabel('R^2')

subplot(3,1,1)
scatter(dtLandADV, rR(1,:)); hold on; 
scatter(dtLandADV, rR(2,:)); hold on;
scatter(dtLandADV, rR(3,:));
ylabel('R')

subplot(3,1,2)
scatter(dtLandADV, rFinIG(1,:)); hold on; 
scatter(dtLandADV, rFinIG(2,:)); hold on;
scatter(dtLandADV, rFinIG(3,:));
ylabel('F^+ (m^3/s)')

subplot(3,1,3)
scatter(dtLandADV, rFoutIG(1,:)); hold on; 
scatter(dtLandADV, rFoutIG(2,:)); hold on;
scatter(dtLandADV, rFoutIG(3,:));
ylabel('F^- (m^3/s)')

% plot IG wave heights for each method
figure
subplot(1,3,1)
plot(dtLandADV, rHm0Inf(1,:), dtLandADV, rHigIN(1,:))

subplot(1,3,2)
plot(dtLandADV, rHm0Inf(2,:), dtLandADV, rHigIN(2,:))

subplot(1,3,3)
plot(dtLandADV, rHm0Inf(3,:), dtLandADV, rHigIN(3,:))

% plot the normalized bed slope for each band
figure
scatter(rBetaH(1,:), rR(1,:)); hold on; 
scatter(rBetaH(2,:), rR(2,:)); hold on; 
scatter(rBetaH(3,:), rR(3,:)); hold on; 
% van Dongeren theoretical R
rBh = 0:0.1:7;
plot(rBh, 0.2*pi*rBh.^2); hold on
plot(rBh, ones(size(rBh)))
ylim([0 1.8])

%% Hss+/h and Hinf+/h

% change to the data directory on external hard-drive 
cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVE_STATS/
load('mtmWaveStats_17min_0pt003_FINAL') 

dtLandADV  = repmat(datetime('25-Aug-2017 22:00:00'), ...
                    size(stMTMWaveStat0pt003{1,3}.dtPb(1,:),2), 1);             
dtLandADV  = hours(stMTMWaveStat0pt003{1,3}.dtPb(1,:)' - dtLandADV);

% separate Hss+ and Hig+
rHss = stMTMWaveStat0pt003{1,3}.Hm0SS;   % total
rHigIN = stMTMWaveStat0pt003{2,3}.Hm0Inf;  % incoming only
h = stMTMWaveStat0pt003{1,3}.h;

% plot
subplot(4,1,3)
yyaxis left
scatter(dtLandADV, rHss./h)
hold on
scatter(dtLandADV, rHigIN./h)
hold on
yyaxis right
scatter(dtLandADV, rHigIN ./ rHssIN)

% calculate the mean rHss/h ratio for each period of the bispectra for the
% IG wave paper
rHssH = rHss./h;
mean(rHssH(17:29))    % -19 to -15 hrs
mean(rHssH(71:86))    % -3 to 1 hrs
mean(rHssH(253:267))

%% directional spreading - ADV puv method - 17 min windows
%%% NOTE: 6/30/19 - I'm not confident in this code; it needs to be tested 
%%% against directional spectra outputs using the PCA method 
%%% (see Ruessink's code) as well as by the maximum entropy method (see Ad
%%% Reniers code for generating directional spectra)

%%% use the multitaper eigenfunctions to generate cospectra (can't get 
%%% directly from mtm code), so naturally these estimates are non-adaptive

% change to the data directory on external hard-drive 
cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV/DSPREAD
load('IN-mtmDSpect_ADVnuv_17min')

% sea-swell indices
iFmin  = 42;   
iFmax  = 257;  

% power-spectra (double-sided)
% NOTE: I previously checked using the non-adaptive power spectra that the 
% following code is correct for creating the mean power spectra directly 
% from the eigencoefficients
rSuu = mean(abs(Yku).^2, 2);
rSvv = mean(abs(Ykv).^2, 2);
rSnn = mean(abs(Ykn).^2, 2);

% Cpu
zSnu = mean(Ykn .* conj(Yku), 2); % cross-spectra
rSnu = real(zSnu);                % co-spectrum

% Cpv
zSnv = mean(Ykn .* conj(Ykv), 2); % cross-spectra
rSnv = real(zSnv);                % co-spectrum

% reshape arrays so they are 2d and multiply by 2 to make them one-sided
% (doesn't actually matter if they are double-sided)
rSuu = reshape(rSuu, [length(rSuu) size(rSuu,3)]) * 2;
rSvv = reshape(rSvv, [length(rSvv) size(rSvv,3)]) * 2;
rSnn = reshape(rSnn, [length(rSnn) size(rSnn,3)]) * 2;
rSnu = reshape(rSnu, [length(rSnu) size(rSnu,3)]) * 2;
rSnv = reshape(rSnv, [length(rSnv) size(rSnv,3)]) * 2;

% METHOD 1: from Nortek manual 
% calculate a1 and b1 for directional spreading 
rA1 = rSnu ./ (( rSnn .* (rSuu + rSvv) ).^(1/2));
rB1 = rSnv ./ (( rSnn .* (rSuu + rSvv) ).^(1/2));

% integrate over the sea-swell band
%rA1 = trapz(f(iFmin:iFmax,1), rA1(iFmin:iFmax,:));
%rB1 = trapz(f(iFmin:iFmax,1), rB1(iFmin:iFmax,:));

% directional spreading
rR1 = sqrt(rA1.^2 + rB1.^2);
rSigma = sqrt(2 * (1 - rR1));
rDeg = rad2deg(trapz(f(iFmin:iFmax,1), rSigma(iFmin:iFmax,:)));
%rDeg = rad2deg(rSigma);

% METHOD 2: Herbers et al. [1999] that incorporates Cuv (in Appendix)
%%% NOTE: the order of integration here is super confusing, but the numbers
%%% are most reasonable (at least in comparison to Herbers et al. [1999])
%%% when you integrate at the end
zSuv = mean(Yku .* conj(Ykv), 2); % cross-spectra
rSuv = real(zSuv);                % co-spectrum
rSuv = reshape(rSuv, [length(rSuv) size(rSuv,3)]) * 2; % make one-sided

% find the mean direction
rThetaMean = atan( (2 * rSuv) ./ (rSuu-rSvv) ) / 2;
%rThetaMean = trapz(f(iFmin:iFmax,1), rThetaMean(iFmin:iFmax,:));

% fourier coefficients 
rA1 = rSnu ./ (( rSnn .* (rSuu + rSvv) ).^(1/2));
rB1 = rSnv ./ (( rSnn .* (rSuu + rSvv) ).^(1/2));
%rA1 = trapz(f(iFmin:iFmax,1), rA1(iFmin:iFmax,:));
%rB1 = trapz(f(iFmin:iFmax,1), rB1(iFmin:iFmax,:));

% standard deviation of the directional distribution (spread)
rSigma = sqrt( 2 * (1 - ( (rA1 .* cos(rThetaMean)) + (rB1 .* sin(rThetaMean)) ) ) );

% integrate at the end
rDeg2 = rad2deg(trapz(f(iFmin:iFmax,1), rSigma(iFmin:iFmax,:)));

% find mean directional spread for each period used in bispectral analysis
mean(rDeg2(17:29))    % -19 to -15 hrs
mean(rDeg2(71:86))    % -3 to 1 hrs
mean(rDeg2(253:267))

%% check to see if I fucked up with NDBC time zone data

cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/NDBC/Freeport_80m

% read offshore data for all of 2017 (downloaded on 01/23/19 with data in
% UTC: https://www.ndbc.noaa.gov/measdes.shtml)
sFileName = '4201982017_SpectralWave_2017.txt';  
nNumberCols = 52;
fid = fopen(sFileName);
% header = frequency bins
format = ['%*s %*s %*s %*s %*s' repmat('%f', [1 nNumberCols-5])];
stNDBCtest.f = cell2mat(textscan(fid, format, 1))';
% first five columns are YY MM DD HH MM
format = repmat('%f', [1 nNumberCols]);
cDATA = textscan(fid, format);
fclose(fid);

% save date and time to datetime variable 
stNDBCtest.dt = datetime([cDATA{1}, cDATA{2}, cDATA{3}, cDATA{4}, ...
                      cDATA{5}, zeros(size(cDATA{1}))])';       
stNDBCtest.psd = cell2mat(cDATA(6:end))';

cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/NDBC
load('NDBC.mat')

% now plot a psd for the "same" time periods
iNDBClow = 5618; % 573 = 8/24 23:00 
%iNDBCup  = 667; % 662 = 8/28 16:00 

figure;
plot(stNDBCtest.f(:,1), stNDBCtest.psd(:,5620))
hold on
plot(stNDBC{1,1}.f(:,1), stNDBC{1,1}.psd(:,575))

%%% they are the same, therefore I fucked up and have been treating the
%%% NDBC data as LST even though it is UTC (+5 hours for August 2017); aye
%%% carumba!

%% check for shear motions in the IG and VLF band

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV/DSPREAD
%load('IN-mtmDSpect_ADVnuv_8minZP.mat')
%load('IN-mtmDSpect_ADVnuv_17minZP.mat')
%load('IN-mtmDSpect_ADVnuv_68min.mat') 
load('IN-mtmDSpect_ADVnuv_HPF_68min.mat') 
%load('IN-mtmDSpect_ADVnuv_137min.mat')

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/VELOCITY
%load('VECTwin_8min.mat')
%load('VECTwin_17min.mat')
load('VECTwin_68min.mat')
%load('VECTwin_137min.mat')
dt = stVECTwin.dt;
dtLandADV = repmat(datetime('25-Aug-2017 22:00:00', 'TimeZone', ...
                    'America/Chicago'), size(dt(1,:),2), 1);             
dtLandADV = hours(dt(1,:)' - dtLandADV)';

% indices for integration over the IG band
iIGlow  = 13%4;%13;%7;%4;
iIGhigh = 165%42;%165;%83;%42;
nG = 9.8; 
rH = repmat(h, size(f,1),1);

% calculate the variance for u, v, and n (single sided spectra)
rM0n = trapz(f(iIGlow:iIGhigh,1), snn(iIGlow:iIGhigh,:));  
rM0u = trapz(f(iIGlow:iIGhigh,1), suu(iIGlow:iIGhigh,:));  
rM0v = trapz(f(iIGlow:iIGhigh,1), svv(iIGlow:iIGhigh,:));  

% R by Lippman et al., (1999) - normalized ratio of total velocity to 
% pressure (n) variance
rR = ( (rM0u + rM0v) ./ rM0n) .* (h / nG);

% fraction of alongshore to cross-shore velocity at IG frequencies
r2D = rM0v ./ rM0u;

% fraction of the velocity variance contributed by shear waves
rAlpha = 1 - (1 ./ rR);

% plot
figure; 
plot(dtLandADV, r2D)
% subplot(2,1,1)
% plot(dtLandADV, stVECTwin.rmsU, dtLandADV, stVECTwin.rmsV)
% xlabel('Hours relative to hurricane landfall')
% ylabel('RMS velocity (m/s)')
% legend('u', 'v')
% plotFancyAxis

figure; 
subplot(2,1,1)
%yyaxis left
plot(dtLandADV, rR); 
ylabel('R')
xlabel('Hours relative to hurricane landfall')
%yyaxis right
subplot(2,1,2)
plot(dtLandADV, rAlpha)
xlabel('Hours relative to hurricane landfall')
ylabel('\alpha')
plotFancyAxis

% equivalent velocity spectra
figure;
for i =  10 % OPTIONS ARE 3,4,10,12,21 for *significant* VLF wave motions
    
    subplot(2,2,[1 2])
    plot(stVECTwin.dt(:,i)+hours(5), stVECTwin.seapress(:,i)) % UTC
    
    rVarUV = suu(:,i) + svv(:,i);
    rVarEta = snn(:,i) .* (nG / h(:,i));  % equivalent velocity
    nR = trapz(f(iIGlow:iIGhigh,1), rVarUV(iIGlow:iIGhigh)) / trapz(f(iIGlow:iIGhigh,1), rVarEta(iIGlow:iIGhigh));
    
    subplot(2,2,3)
    loglog(f(:,i), rVarUV, f(:,i), rVarEta); hold on
    loglog([0.003 0.003], [0.01 100]);
    xlim([0 0.04])
    legend('Suu + Svv', 'S\eta\eta')
    plotFancyAxis 
    
end    

%% sensitivity figure for choice of NW and K -- ADV --
%%% NOTE, none of these spectra are zero-padded and they are not corrected
%%% for attenuation of the pressure signal by linear wave theory

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV

%%% --------------------------- 8 min --------------------------------- %%%
% plot Burst 141-168 (landfall, energy in low-frequency IG mode)
figure;
iBurst = 141;
nNW = [3 4 5 6 7 8]; 
n2W = nNW/(8192*(1/16))*2; % W = NW/(N*deltaT);

subplot(1,4,1)
load('OUT-mtmADVsp_nPLUS_8minNW3K5.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmADVsp_nPLUS_8minNW4K7.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmADVsp_nPLUS_8minNW5K9.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_8minNW6K11.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_8minNW7K13.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_8minNW8K15.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
loglog([0.04 0.04], [0.01 10]); hold on
xlim([0 0.25])
ylim([0.01 10])
xticks([0 0.01 0.1])
xlabel('Frequency (Hz)')
ylabel('PSD (m^2/Hz)')
legend('NW=3,K=5', ...
       'NW=4,K=7', ...
       'NW=5,K=9', ...
       'NW=6,K=11', ...
       'NW=7,K=13', ...
       'NW=8,K=15')
%legend('NW=3,K=5,dof=10,2W=0.0117', ...
%       'NW=4,K=7,dof=14,2W=0.0156', ...
%       'NW=5,K=9,dof=18,2W=0.0195', ...
%       'NW=6,K=11,dof=22,2W=0.0234', ...
%       'NW=7,K=13,dof=26,2W=0.0273', ...
%       'NW=8,K=15,dof=30,2W=0.0312')
% formatting   
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.03 .03] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 14        , ...
  'FontName'    , 'Arial'   );  
   
%%% -------------------------- 17 min --------------------------------- %%%
% plot Burst 70 (landfall, energy in low-frequency IG mode)

iBurst = 71;
nNW = [3 4 5 6 7 8]; 
n2W = nNW/(16384*(1/16))*2; % W = NW/(N*deltaT);

subplot(1,4,2)
load('OUT-mtmADVsp_nPLUS_17minNW3K5.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmADVsp_nPLUS_17minNW4K7.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmADVsp_nPLUS_17minNW5K9.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_17minNW6K11.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_17minNW7K13.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_17minNW8K15.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
loglog([0.04 0.04], [0.01 10]); hold on
xlim([0 0.25])
ylim([0.01 10])
xticks([0 0.001 0.01 0.1])
xlabel('Frequency (Hz)')
ylabel('PSD (m^2/Hz)')
%legend('NW=3,K=5,dof=10,2W=0.0059', ...
%       'NW=4,K=7,dof=14,2W=0.0078', ...
%       'NW=5,K=9,dof=18,2W=0.0098', ...
%       'NW=6,K=11,dof=22,2W=0.0117', ...
%       'NW=7,K=13,dof=26,2W=0.0137', ...
%       'NW=8,K=15,dof=30,2W=0.0156')
legend('NW=3,K=5', ...
       'NW=4,K=7', ...
       'NW=5,K=9', ...
       'NW=6,K=11', ...
       'NW=7,K=13', ...
       'NW=8,K=15')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.03 .03] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 14        , ...
  'FontName'    , 'Arial'   );

%%% -------------------------- 34 min --------------------------------- %%%
% plot Burst 30ish? (landfall, energy in low-frequency IG mode)

iBurst = 28;
nNW = [3 4 5 6 7 8]; 
n2W = nNW/(32768*(1/16))*2; % W = NW/(N*deltaT);

subplot(1,4,3)
load('OUT-mtmADVsp_nPLUS_34minNW3K5.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmADVsp_nPLUS_34minNW4K7.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmADVsp_nPLUS_34minNW5K9.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_34minNW6K11.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_34minNW7K13.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_34minNW8K15.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
loglog([0.04 0.04], [0.01 10]); hold on
xlim([0 0.25])
ylim([0.01 10])
xticks([0 0.001 0.01 0.1])
xlabel('Frequency (Hz)')
ylabel('PSD (m^2/Hz)')
%legend('NW=3,K=5,dof=10,2W=0.0059', ...
%       'NW=4,K=7,dof=14,2W=0.0078', ...
%       'NW=5,K=9,dof=18,2W=0.0098', ...
%       'NW=6,K=11,dof=22,2W=0.0117', ...
%       'NW=7,K=13,dof=26,2W=0.0137', ...
%       'NW=8,K=15,dof=30,2W=0.0156')
legend('NW=3,K=5', ...
       'NW=4,K=7', ...
       'NW=5,K=9', ...
       'NW=6,K=11', ...
       'NW=7,K=13', ...
       'NW=8,K=15')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.03 .03] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 14        , ...
  'FontName'    , 'Arial'   ); 

%%% -------------------------- 137 min -------------------------------- %%%   
% plot Burst 10 (many instances of VLF signals)
iBurst = 10;
nNW = [3 4 5 6 7 8]; 
n2W = nNW/(131072*(1/16))*2; % W = NW/(N*deltaT);

subplot(1,4,4)
load('OUT-mtmADVsp_nPLUS_137minNW3K5.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmADVsp_nPLUS_137minNW4K7.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmADVsp_nPLUS_137minNW5K9.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_137minNW6K11.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_137minNW7K13.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmADVsp_nPLUS_137minNW8K15.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
loglog([0.04 0.04], [0.01 10]); hold on
loglog([0.003 0.003], [0.01 10]); hold off
xlim([0 0.25])
ylim([0.01 10])
xticks([0 0.001 0.01 0.1])
xlabel('Frequency (Hz)')
ylabel('PSD (m^2/Hz)')
%legend('NW=3,K=5,dof=10,2W=0.0007', ...
%       'NW=4,K=7,dof=14,2W=0.0010', ...
%       'NW=5,K=9,dof=18,2W=0.0012', ...
%       'NW=6,K=11,dof=22,2W=0.0015', ...
%       'NW=7,K=13,dof=26,2W=0.0017', ...
%       'NW=8,K=15,dof=30,2W=0.0020')
legend('NW=3,K=5', ...
       'NW=4,K=7', ...
       'NW=5,K=9', ...
       'NW=6,K=11', ...
       'NW=7,K=13', ...
       'NW=8,K=15')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.03 .03] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 14        , ...
  'FontName'    , 'Arial'   ); 

%% sensitivity figure for choice of NW and K -- RRUs --
%%% NOTE, none of these spectra are zero-padded and they are not corrected
%%% for attenuation of the pressure signal by linear wave theory
%%% ------------------------ 17 min - RRU1 ---------------------------- %%%

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/SensitivityTests/RRU1

% plot Burst 55 (during overwash)

figure;
iBurst = 59; %59
nNW = [3 4 5 6 7 8]; 
n2W = nNW/(16384*(1/16))*2; % W = NW/(N*deltaT);

subplot(2,3,1)
load('OUT-mtmRRU1_n_17minZP_NW3K5.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU1_n_17minZP_NW4K7.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU1_n_17minZP_NW5K9.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_17minZP_NW6K11.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_17minZP_NW7K13.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_17minZP_NW8K15.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
xlim([0 0.04])
xlabel('Frequency (Hz)')
ylabel('PSD (m^2/Hz)')
legend('NW=3,K=5', ...
       'NW=4,K=7', ...
       'NW=5,K=9', ...
       'NW=6,K=11', ...
       'NW=7,K=13', ...
       'NW=8,K=15')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.03 .03] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 14        , ...
  'FontName'    , 'Arial'   );

%%% ------------------------ 34 min - RRU1 ---------------------------- %%%

iBurst = 26; 
nNW = [3 4 5 6 7 8]; 
n2W = nNW/(32768*(1/16))*2; % W = NW/(N*deltaT);

subplot(2,3,2)
load('OUT-mtmRRU1_n_34minNW3K5.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU1_n_34minNW4K7.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU1_n_34minNW5K9.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_34minNW6K11.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_34minNW7K13.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_34minNW8K15.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
xlim([0 0.04])
xlabel('Frequency (Hz)')
ylabel('PSD (m^2/Hz)')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.03 .03] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 14        , ...
  'FontName'    , 'Arial'   );

%%% ------------------------ 68 min - RRU1 ---------------------------- %%%
% plot Burst 13 (during overwash)

iBurst = 13;
nNW = [3 4 5 6 7 8]; 
n2W = nNW/(16384*(1/16))*2; % W = NW/(N*deltaT);

subplot(2,3,3)
load('OUT-mtmRRU1_n_68min_NW3K5.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU1_n_68min_NW4K7.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU1_n_68min_NW5K9.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_68min_NW6K11.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_68min_NW7K13.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_68min_NW8K15.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
xlim([0 0.04])
xlabel('Frequency (Hz)')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.03 .03] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 14        , ...
  'FontName'    , 'Arial'   );

%%% ------------------------ 17 min - RRU2 ---------------------------- %%%
% plot Burst 55 (during overwash)

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/SensitivityTests/RRU2

iBurst = 59;
nNW = [3 4 5 6 7 8]; 
n2W = nNW/(2048*(1/2))*2; % W = NW/(N*deltaT);

subplot(2,3,4)
load('OUT-mtmRRU2_n_17min_NW3K5.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU2_n_17min_NW4K7.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU2_n_17min_NW5K9.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU2_n_17min_NW6K11.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU2_n_17min_NW7K13.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU2_n_17min_NW8K15.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
xlim([0 0.04])
xlabel('Frequency (Hz)')
ylabel('PSD (m^2/Hz)')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.03 .03] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 14        , ...
  'FontName'    , 'Arial'   );

%%% ------------------------ 34 min - RRU2 ---------------------------- %%%

iBurst = 26; 
nNW = [3 4 5 6 7 8]; 
n2W = nNW/(4096*(1/2))*2; % W = NW/(N*deltaT);

subplot(2,3,5)
load('OUT-mtmRRU1_n_34min_NW3K5.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU1_n_34min_NW4K7.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU1_n_34min_NW5K9.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_34min_NW6K11.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_34min_NW7K13.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU1_n_34min_NW8K15.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
xlim([0 0.04])
xlabel('Frequency (Hz)')
ylabel('PSD (m^2/Hz)')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.03 .03] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 14        , ...
  'FontName'    , 'Arial'   );

%%% ------------------------ 68 min - RRU2 ---------------------------- %%%
% plot Burst 13 (during overwash)

iBurst = 13;
nNW = [3 4 5 6 7 8]; 
n2W = nNW/(8192*(1/2))*2; % W = NW/(N*deltaT);

subplot(2,3,6)
load('OUT-mtmRRU2_n_68min_NW3K5.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU2_n_68min_NW4K7.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold  on
load('OUT-mtmRRU2_n_68min_NW5K9.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU2_n_68min_NW6K11.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU2_n_68min_NW7K13.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
load('OUT-mtmRRU2_n_68min_NW8K15.mat')
loglog(f(:,iBurst), snn(:,iBurst)); hold on
xlim([0 0.04])
xlabel('Frequency (Hz)')
ylabel('PSD (m^2/Hz)')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.03 .03] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 14        , ...
  'FontName'    , 'Arial'   );

%% prewhiten pressure data using Bob Parker's algorithm

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/INPUTS
%cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV/INPUTS

%rP = load('ADVsp_128min_psin.txt');
rP = load('RRU2_34min_psin.txt');
sPath = '/Users/KatherineAnardeWheels/Research/Software/BDW_Scripts/psd';
sFile = 'temp';
bPwht = true;
nFsamp = 2;

for i = 29%1: size(rP,2)
    [rSpp, rF, rPwht] = runPsd(rP(:,i), nFsamp, sFile, sPath, 5, bPwht, false, []);
end

figure; loglog(rF, rSpp); xlim([0 0.04])

dlmwrite('RRU1_68min_psin_PW.txt', rPwht, 'delimiter', ...
         '\t', 'precision', 12)
     
%% try creating a 68 min spectrogram from the hpf data 

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV
load('dtADV_68min.mat');
load('OUT-etaIN_68min.mat'); % this was created adhoc using the code in R

% try a 3 pole butterworth filter
nFsamp = 16;
nNy = nFsamp/2;
[b, a] = butter(3, 0.002/nNy, 'high'); % get the transfer function
rHPF = filtfilt(b, a, etaIN);

% let's see what the spectrum looks like
sPath = '/Users/KatherineAnardeWheels/Research/Software/BDW_Scripts/psd';
sFile = 'temp';
bPwht = true;

for i = 1 : size(rHPF,2)
    [rSnn(:,i), rF(:,i)] = runPsd(rHPF(:,i), nFsamp, sFile, sPath, 8, bPwht, false, []);
end
    
figure; loglog(rF(:,4), rSnn(:,4)); xlim([0.003 0.25])
figure; loglog(rF(:,14), rSnn(:,14)); xlim([0.003 0.25])

figure;
imagesc(dtLandADV, rF(1:165,1), rSnn(1:165,:)) 
% formatting
    axis tight
    colormap('jet')
    c = colorbar;
    ylabel(c, 'S\eta\eta (m^2/Hz)', 'FontSize', 13)
    ylabel('Frequency (Hz)')
    xlabel('Hours relative to hurricane landfall')
    set(gca, 'YGrid', 'off', 'YLim', [0 0.04])
    set(gca, 'YDir', 'normal')
    set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'LineWidth'   , 1         , ...
      'FontSize'    , 12        , ...
      'FontName'    , 'Arial'   ); 
  
%% high-pass filter the entire ADV etaIN time record with butterworth filter
% save 8 min windows for bispectral analysis and 68 min windows for stats
% NOTE: modified on 6/12/19 to hpf other ADV and RRU pressure records

% cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV
% load('OUT-etaIN_68min.mat');
% etaIN = etaIN(:);  

% also HPF the sea pressure time series for the ADV
% cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV/INPUTS
% rSP = load('ADVsp_68min_psin.txt'); % already in Pascals
% rSP = rSP(:);

% and RRU-1/RRU-2/RRU-6
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/INPUTS
rSP = load('RRU1_68min_psin.txt'); % in Pascals
rSP = rSP(:);

% and QC'ed velocity data
% cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/VELOCITY/
% load('stVECTwin')
% rVelU = stVECTwin.velU(:);
% rVelV = stVECTwin.velV(:);

% apply a 3 pole butterworth HP filter
nFsamp = 16;
nNy = nFsamp/2;
[b, a] = butter(3, 0.002/nNy, 'high'); % get the transfer function
%rHPF = filtfilt(b, a, etaIN);
rHPFsp = filtfilt(b, a, rSP);
% rHPFvelu = filtfilt(b, a, rVelU);
% rHPFvelv = filtfilt(b, a, rVelV);
    
% create windows
% inputs
n2nRRU = 13;             
nWin   = 2^n2nRRU;
%nTotal = length(rHPFvelu);    % total # of indices 
nTotal = length(rHPFsp);     % total # of indices 
%nTotal = length(rHPF);       % total # of indices 
nBurst = floor(nTotal/nWin);  % number of bursts (should be >75)
cntWin = 0;                   % initialize counter 

for iBurst = 1 : nBurst

    %rEtaIN(:,iBurst)  = rHPF( (1 + cntWin) : (1 + cntWin + nWin - 1) );
    rSPhpf(:,iBurst)  = rHPFsp( (1 + cntWin) : (1 + cntWin + nWin - 1) );
    %rVELUhpf(:,iBurst)  = rHPFvelu( (1 + cntWin) : (1 + cntWin + nWin - 1) );
    %rVELVhpf(:,iBurst)  = rHPFvelv( (1 + cntWin) : (1 + cntWin + nWin - 1) );
    
    % increase window counter
    cntWin = cntWin + nWin;
end

% dlmwrite('ADVetaIN_68min_hpf_0pt002.txt', rEtaIN, 'delimiter', ...
%          '\t', 'precision', 12)    
% dlmwrite('ADVsp_68min_hpf_0pt002.txt', rSPhpf, 'delimiter', ...
%          '\t', 'precision', 12)     
% dlmwrite('RRU2_8min_hpf_0pt002.txt', rSPhpf, 'delimiter', ...
%          '\t', 'precision', 12) 

nRho = 1025;   % seawater density [kg/m3]
nG   = 9.81;   % gravitational constant [m/s2]
rEtaHPF = rSPhpf / (nRho*nG);
dlmwrite('RRU1_8min_eta_hpf_0pt002.txt', rEtaHPF, 'delimiter', ...
         '\t', 'precision', 12) 

% dlmwrite('IN-VECTwinHPF_0pt002Hz_U_68min.txt', rVELUhpf, 'delimiter', ...
%          '\t', 'precision', 12) 
% dlmwrite('IN-VECTwinHPF_0pt002Hz_V_68min.txt', rVELVhpf, 'delimiter', ...
%          '\t', 'precision', 12) 

%% quantify assymetry for the HPF IG wave time series at Matagorda
%%% calculate assymetry using the hilbert transform and 34 minute time
%%% records

nRho = 1025;   % seawater density [kg/m3]
nG   = 9.81;   % gravitational constant [m/s2]

% load the HPF time series of pressure in Pascals
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/INPUTS
load('dtRRU1_34min.mat')

% full time record in the backshore
rSP = load('RRU1_34min_hpf_0pt002.txt'); % already in Pascals
%rSP = load('RRU1_68min_hpf_0pt002.txt');
rEtaHPF = rSP / (nRho * nG);
rEtaHPF = detrend(rEtaHPF);
rHt = imag(hilbert(rEtaHPF));
rA1 = mean(rHt.^3) ./ (mean(rEtaHPF.^2)).^(3/2);
figure; plot(dtLandRRU1, rA1)
mean(rA1(1:25)) % shallow flooding
mean(rA1(26:35)) % overwash

% low-pass filtered time record in the backshore
% apply a 3 pole butterworth LP filter
nFsamp = 16;
nNy = nFsamp/2;
[b, a] = butter(3, 0.04/nNy, 'low'); % get the transfer function
rEtaLPF = filtfilt(b, a, rEtaHPF);
rEtaLPF = detrend(rEtaLPF);
rHt = imag(hilbert(rEtaLPF));
rA1lpf = mean(rHt.^3) ./ (mean(rEtaLPF.^2)).^(3/2);
figure; plot(dtLandRRU1, rA1, dtLandRRU1, rA1lpf)
mean(rA1lpf(1:25)) % shallow flooding
mean(rA1lpf(26:35)) % overwash

% full time record in the back-barrier
rSP = load('RRU2_34min_hpf_0pt002.txt'); % already in Pascals
rEtaHPF = rSP / (nRho * nG);
rEtaHPF = detrend(rEtaHPF);
rHt = imag(hilbert(rEtaHPF));
rA2 = mean(rHt.^3) ./ (mean(rEtaHPF.^2)).^(3/2);
figure; plot(dtLandRRU1, rA2)

%% check amplitude modulation of IG waves by VLF waves in the backshore/back-barrier

nRho = 1025;   % seawater density [kg/m3]
nG   = 9.81;   % gravitational constant [m/s2]

% load the HPF time series of pressure in Pascals
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/INPUTS
%rSP = load('RRU1_34min_hpf_0pt002.txt'); % already in Pascals
rSP = load('RRU2_34min_hpf_0pt002.txt'); % already in Pascals
%load('dtRRU1_34min.mat')
rEtaHPF = rSP / (nRho * nG);

% now let's just look at overwash
iBlow = 26; % was 25 but distracting
iBup  = 28;
rEtaHPF = rEtaHPF(:,iBlow:iBup);
rEtaHPF = rEtaHPF(:);
%dt = dt(:,iBlow:iBup); dt = dt(:);

% get the wave envelope for each record
rEnv = abs(hilbert(rEtaHPF));

% now low pass filter to check for envelope modulations at VLF
nFsamp = 2;   % apply a 3 pole butterworth HP filter
nNy = nFsamp/2;
[b, a] = butter(3, 0.001/nNy, 'low'); % get the transfer function
rLPF = filtfilt(b, a, rEnv);

% plot to check
%figure; plot(dt, rEtaHPF); hold on; plot(dt, rEnv); hold on; plot(dt, rLPF)
figure; plot(rEtaHPF); hold on; plot(rEnv); hold on; plot(rLPF)

% now find the max amplitude modulation (absolute value)
max(abs(detrend(rLPF)))
%figure; plot(dt, detrend(rLPF))
figure; plot(detrend(rLPF))


%% check amplitude modulation of IG waves by VLF waves in the surf zone

% load the HPF, etaIN time series
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV/INPUTS
rEtaHPF = load('ADVetaIN_68min_hpf_0pt002.txt'); 
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV
load('dtADV_68min'); 

% bandpass filter to just get the IG waves
nFsamp = 16;   % apply a 3 pole butterworth filter
nNy = nFsamp/2;
high_corner = 0.045;
low_corner  = 0.003;
[b, a] = butter(3, [low_corner high_corner]/nNy); % get the transfer function
rBPF = filtfilt(b, a, rEtaHPF);

% get the wave envelope for each record
rEnv = abs(hilbert(rBPF));

% now low pass filter to check for envelope modulations at VLF
[b, a] = butter(3, 0.002/nNy, 'low'); % get the transfer function
rLPF = filtfilt(b, a, rEnv);

% plot to check
iBurst = 20;
% inputs for windowing the PT time series LST/LDT
iMin = find(stVLF.HPF.dt == dt(1,iBurst));
iMax = find(stVLF.HPF.dt == dt(end,iBurst));
figure; 
subplot(3,1,1); plot(stVLF.HPF.dt(iMin:iMax), stVLF.HPF.eta(iMin:iMax))
subplot(3,1,2); plot(dt(:,iBurst), rBPF(:,iBurst), dt(:,iBurst), rEnv(:,iBurst), dt(:,iBurst), rLPF(:,iBurst))

% now find the max amplitude modulation (absolute value)
mean(max(abs(detrend(rLPF(:,17:22))))) % this spans landfall

figure; plot(dt(:,iBurst), detrend(rLPF(:,iBurst)))

%% Test on methodology for obtaining the HPF time series

nRho = 1025;   % seawater density [kg/m3]
nG   = 9.81;   % gravitational constant [m/s2]

% load the HPF time series of pressure in Pascals
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/INPUTS
load('dtRRU1_34min.mat')

% full time record in the backshore
rSPhpf = load('RRU1_34min_hpf_0pt002.txt'); % already in Pascals
rSP = load('RRU1_34min_psin.txt');

% method 1: using a hpf
rEtaHPF = rSPhpf / (nRho * nG);
rEtaHPF = detrend(rEtaHPF);

% method 2: subtracting the LPF time series
rEta = rSP / (nRho * nG);
nFsamp = 16;   % apply a 3 pole butterworth filter
nNy = nFsamp/2;
[b, a] = butter(3, 0.002/nNy, 'low'); % get the transfer function
rLPF = filtfilt(b, a, rEta);
rEtaHPF2 = rEta-rLPF;
rEtaHPF2 = detrend(rEtaHPF2);

% plot
i = 28;
figure; plot(dt(:,i), rEtaHPF(:,i), dt(:,i), rEtaHPF2(:,i))

% now test the ADV data
% load the HPF time series of pressure in Pascals
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV
load('dtADV_68min.mat')
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV/INPUTS
rSPhpf = load('ADVsp_68min_hpf_0pt002.txt'); % already in Pascals
rSP = load('ADVsp_68min_psin.txt');

% method 1: using a hpf
rEtaHPF = rSPhpf / (nRho * nG);
rEtaHPF = detrend(rEtaHPF);

% method 2: subtracting the LPF time series
rEta = rSP / (nRho * nG);
nFsamp = 16;   % apply a 3 pole butterworth filter
nNy = nFsamp/2;
[b, a] = butter(3, 0.002/nNy, 'low'); % get the transfer function
rLPF = filtfilt(b, a, rEta);
rEtaHPF2 = rEta-rLPF;
rEtaHPF2 = detrend(rEtaHPF2);

% plot
i = 14; %14
figure; plot(dt(:,i), rEtaHPF(:,i), dt(:,i), rEtaHPF2(:,i))
figure; plot(dt(:,i), rEtaHPF(:,i)-rEtaHPF2(:,i))

% ok, so the HPF messes with the edges of the signal (this we know), but doesn't really
% appear to screw with the amplitudes of the IG waves; let's check the
% spectra just to be sure
sPath = '/Users/KatherineAnardeWheels/Research/Software/BDW_Scripts/psd';
sFile = 'temp';
bPwht = true;
nFsamp = 16;
NW = 8;

% note, this is the 1 std error

for iBurst = 1 : size(rEtaHPF,2)
    [rSnnHPF(:,iBurst), rF] = runPsd(rEtaHPF(:,iBurst), nFsamp, sFile, sPath, ...
                                                     NW, bPwht, false, []);
end

for iBurst = 1 : size(rEtaHPF2,2)
    [rSnnHPF2(:,iBurst), rF] = runPsd(rEtaHPF2(:,iBurst), nFsamp, sFile, sPath, ...
                                                     NW, bPwht, false, []);
end

figure;
imagesc(dtLandADV, rF(:,1), rSnnHPF) % only the IG band
    % formatting
    axis tight
    %colormap('jet')
    %cmocean('thermal') 
    c = colorbar;
    ylabel(c, 'S\eta\eta (m^2/Hz)', 'FontSize', 13)
    ylabel('Frequency (Hz)')
    xlabel('Hours relative to hurricane landfall')
    set(gca, 'YGrid', 'off', 'YLim', [0 0.04], 'XLim', [-23 20])
    set(gca, 'YDir', 'normal')
    
figure;
imagesc(dtLandADV, rF(:,1), rSnnHPF2) % only the IG band
    % formatting
    axis tight
    %colormap('jet')
    %cmocean('thermal') 
    c = colorbar;
    ylabel(c, 'S\eta\eta (m^2/Hz)', 'FontSize', 13)
    ylabel('Frequency (Hz)')
    xlabel('Hours relative to hurricane landfall')
    set(gca, 'YGrid', 'off', 'YLim', [0 0.04], 'XLim', [-23 20])
    set(gca, 'YDir', 'normal')    
    
figure;
imagesc(dtLandADV, rF(:,1), rSnnHPF2-rSnnHPF) % only the IG band
    % formatting
    axis tight
    %colormap('jet')
    %cmocean('thermal') 
    c = colorbar;
    ylabel(c, 'S\eta\eta (m^2/Hz)', 'FontSize', 13)
    ylabel('Frequency (Hz)')
    xlabel('Hours relative to hurricane landfall')
    set(gca, 'YGrid', 'off', 'YLim', [0 0.04], 'XLim', [-23 20])
    set(gca, 'YDir', 'normal')    
    
%% Make new high-pass filtered time series using method of subtraction

% RRUs (in Pascals)
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/INPUTS
%rSP = load('RRU1_34min_psin.txt');
rSP = load('RRU2_34min_psin.txt');

% subtracting the LPF time series
%nFsamp = 16;     % apply a 3 pole butterworth filter
nFsamp = 2;     % apply a 3 pole butterworth filter
nNy = nFsamp/2;
[b, a] = butter(3, 0.002/nNy, 'low'); % get the transfer function
rLPF = filtfilt(b, a, rSP);
rHPF = rSP-rLPF;
figure; plot(rHPF(:,1))

% save 34 min windows for stats (in Pascals)
dlmwrite('RRU2_34min_hpf_v2_0pt002.txt', rHPF, 'delimiter', ...
         '\t', 'precision', 12)

% save 68 and 8 min windows 
% ADV
% cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV
% load('OUT-etaIN_68min.mat'); 

% and RRU-1/RRU-2/RRU-6
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/INPUTS
%rSP = load('RRU1_68min_psin.txt'); % in Pascals
rSP = load('RRU2_68min_psin.txt'); % in Pascals
%rSP = load('RRU6_68min_psin.txt');
%rSP = rSP(:);

% subtracting the LPF time series
%nFsamp = 16;     % apply a 3 pole butterworth filter
nFsamp = 2;     % apply a 3 pole butterworth filter
nNy = nFsamp/2;
[b, a] = butter(3, 0.002/nNy, 'low'); % get the transfer function
rLPF = filtfilt(b, a, rSP);
%rLPF = filtfilt(b, a, etaIN);
rHPF = rSP-rLPF;
%rHPF = etaIN-rLPF;
figure; plot(rHPF(:,1))
rHPF = rHPF(:);
    
% create windows
% inputs
%n2nRRU = 14; % 13 for 8 min windows for 16 Hz data  
n2nRRU = 11; % 10 for 8 min windows for 2 Hz data
nWin   = 2^n2nRRU;
nTotal = length(rHPF);        % total # of indices 
nBurst = floor(nTotal/nWin);  % number of bursts (should be >75)
cntWin = 0;                   % initialize counter 

for iBurst = 1 : nBurst

    %rEtaIN(:,iBurst)  = rHPF( (1 + cntWin) : (1 + cntWin + nWin - 1) );
    rSPhpf(:,iBurst)  = rHPF( (1 + cntWin) : (1 + cntWin + nWin - 1) );
    
    % increase window counter
    cntWin = cntWin + nWin;
end

%dlmwrite('ADVetaIN_8min_hpf_v2_0pt002.txt', rEtaIN, 'delimiter', ...
%         '\t', 'precision', 12)    
% dlmwrite('ADVetaIN_68min_hpf_v2_0pt002.txt', rHPF, 'delimiter', ...
%          '\t', 'precision', 12)    

nRho = 1025;   % seawater density [kg/m3]
nG   = 9.81;   % gravitational constant [m/s2]
rEtaHPF = rSPhpf / (nRho*nG);

%dlmwrite('RRU1_17min_eta_hpf_v2_0pt002.txt', rEtaHPF, 'delimiter', ...
%          '\t', 'precision', 12) 
dlmwrite('RRU2_17min_eta_hpf_v2_0pt002.txt', rEtaHPF, 'delimiter', ...
         '\t', 'precision', 12) 
     
%% check to see if the short wave group structure disappears

% load time series with VLF removed
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV/INPUTS
rEtaInHPF = load('ADVetaIN_68min_hpf_v2_0pt002.txt');

cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV
load('dtADV_68min.mat')

% load time series without VLF removed
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/ADV/INPUTS
%load('OUT-etaIN_68min.mat')
rEta = detrend(load('ADVsp_68min_psin-tide.txt'));
nRho = 1025;   % seawater density [kg/m3]
nG   = 9.81;   % gravitational constant [m/s2]
rEta = rEta / (nRho * nG);

% bandpass filter to just get the sea-swell waves
nFsamp = 16;   % apply a 3 pole butterworth filter
nNy = nFsamp/2;
high_corner = 0.33;
low_corner  = 0.04;
[b, a] = butter(3, [low_corner high_corner]/nNy); % get the transfer function
rSS = filtfilt(b, a, rEtaInHPF);

% get the SS wave envelope and filter at IG frequencies
%rEnvSS = abs(hilbert(detrend(rSS)));
rEnvSS = abs(detrend(rSS)); 
[b, a] = butter(3, 0.04/nNy, 'low'); % get the transfer function
rEnvSS = (pi/2) * filtfilt(b, a, rEnvSS); % I just like this version better
%rEnvSS = filtfilt(b, a, rEnv);
rIG = filtfilt(b, a, rEtaInHPF);

% now only the lowest frequency IG waves
[b, a] = butter(3, 0.015/nNy, 'low'); % get the transfer function
rIGlow = filtfilt(b, a, rEtaInHPF);

% now the vlf envelope of the IG waves and VLF waves
rEnvIG = abs(detrend(rIG)); 
[b, a] = butter(3, 0.003/nNy, 'low'); % get the transfer function
rEnvIG = (pi/2) * filtfilt(b, a, rEnvIG);
%rVLF = filtfilt(b, a, etaIN); % total, not just incoming only
rVLF = filtfilt(b, a, rEta); % total, not just incoming only

% plot the short waves + envelope, IG waves + envelope, VLF waves
figure; 
i = 14;
ax(1) = subplot(3,1,1); plot(dt(:,i)+hours(5), rSS(:,i), dt(:,i)+hours(5), rEnvSS(:,i));
ylabel('\eta+ (m)'); legend('sea/swell', 'SS wave envelope')
plotFancyAxis
ax(2) = subplot(3,1,2); plot(dt(:,i)+hours(5), rIG(:,i), dt(:,i)+hours(5), rIGlow(:,i), dt(:,i)+hours(5), rEnvIG(:,i));
ylabel('\eta+ (m)'); legend('IG', 'IG<0.015 Hz', 'IG wave envelope')
plotFancyAxis
ax(3) = subplot(3,1,3); plot(dt(:,i)+hours(5), rEta(:,i), dt(:,i)+hours(5), rVLF(:,i));
ylabel('\eta (m)'); legend('total', 'VLF')
plotFancyAxis
linkaxes(ax, 'x')

