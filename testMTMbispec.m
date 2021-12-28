 % --------------------------testMTMbispec---------------------------------%
% Purpose: This run script uses synthetic QPC data and Gaussian white noise
% (zero mean unit variance) to analyze the statistics of the multitaper
% bispectral estimator
%
% SEE ALSO: Elgar & Guza (1988), Birkelund et al. (2000)
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       03/15/18        KAnarde             Original code
%
%% preamble
bQPC = false;           % test mtmBispecd vs. bispecd/bicoher for QPC case
bWhiteNoise = false;    % test variance reduction of the MBE for WGN
bConfidence = true;     % evaluate the use of 6/dof for MBE error

%% QPC
% MBE ---------------------------------------------------------------
% N  = 64
% NW = 2
% K  = 4
% nFFT = "default", 128

if bQPC
    % run the adaptive case for all segments (ensemble average)  
    load('qpc3_zp_nw2k4.mat')   
    [zBspecA, rBic1A, rBic2A, rFA] = mtmBispecd(Yk, vk, dk, f, snn, 'false');
    % run only for the first realization
    [zBspecA1, rBic1A1, rBic2A1, rFA1] = mtmBispecd(Yk(:,:,1), vk(:,:,1), ...
                                           dk(:,:,1), f(:,1), ...
                                           snn(:,1), 'false');

    % run the full time series - adaptive (TAKES TOO LONG)
    % load('qpcFull_zp_nw160k320.mat')
    % [zBspecAF, rBic1AF, rBic2AF, rFAF] = mtmBispecd(Yk, vk, dk, f, snn, 'false');

    % % run the non-adaptive case for all segments
    % load('qpc_NAzp_nw3pt5k7.mat') 
    % [zBspecNA, rBic1NA, rBic2NA, rFNA] = mtmBispecdNA(Yk, vk, f, snn, 'false');
    % % run only for the first realization
    % [zBspecNA1, rBic1NA1, rBic2NA1, rFNA1] = mtmBispecdNA(Yk(:,:,1), ...
    %                                             vk(:,:,1), f(:,1), snn(:,1),...
    %                                             'false');                                   

    % bispecd -----------------------------------------------------------
    %load('qpc.mat')
    zdat = qpcgen('default');
    % only the first realization
    %[bd,w,xf] = bispecd(zmat(:,1),128,3,64,0);
    [bd,w,xf] = bispecd(zdat(:,1),128,3,64,0);

    % % this is the raw biperiodogram (very high variance)
    % [~,w,xf] = bispecd(zmat(:,1),128,0,64,0);
    % mask = hankel([1:128],[128,1:128-1]);
    % cxf = conj(xf);
    % bd = fftshift(xf*xf.'.*cxf(mask))*length(zmat(:,1))^2;

    % periodogram, here multiplied by N because xf is normalized by N
    %per = fftshift(xf.*conj(xf)) * length(zmat(:,1)); 
    per = fftshift(xf.*conj(xf)) * length(zdat(:,1)); 

    % the segment average
    %[bdAve,wAve] = bispecd(zmat,128,3,64,0);
    [bdAve,wAve] = bispecd(zdat,128,3,64,0);

    % bicoher -----------------------------------------------------------
    % default hanning window
    %[bic1,bic2,wb1,bspec1] = bicoher(zmat(:,1),128,[],64,0);
    [bic1,bic2,wb1,bspec1] = bicoher(zdat(:,1),128,[],64,0);
    %[bic1Ave,bic2Ave,wbAve,bspecAve,pyyAve] = bicoher(zmat,128,[],64,0);
    [bic1Ave,bic2Ave,wbAve,bspecAve,pyyAve] = bicoher(zdat,128,[],64,0);

    % shift the periodogram
    %pyyAve = fftshift(pyyAve) * length(zmat(:,1));
    pyyAve = fftshift(pyyAve) * length(zdat(:,1)); 

    % PLOTS
    % BURST ----------------------------------------------------------------- %
    % ADAPTIVE
    figure;
    subplot(1,4,1)
       contour(rFA1, rFA1, abs(zBspecA1), 10), grid on 
       title('mtmBispecd - $|{B}|$ Burst 1', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]);
    subplot(1,4,2)
       contour(rFA1, rFA1, imag(zBspecA1), 10), grid on 
       title('$imag(B)$', 'interpreter', 'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]);
    subplot(1,4,3)
       surf(rFA1, rFA1, rBic1A1); view(2); shading interp
       %contour(rFA1, rFA1, rBic1A1, 10), grid on 
       title('$b$', 'interpreter', 'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]); caxis([0 1]);
    subplot(1,4,4)
       surf(rFA1, rFA1, rBic2A1); view(2); shading interp
       %contour(rFA1, rFA1, rBic2A1, 10), grid on 
       title('$b^2$', 'interpreter', 'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]);

    % BISPECD
    figure;
    subplot(1,4,1)
       contour(w, w, abs(bd), 10), grid on 
       title('bispecd - $|B|$ Burst 1', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]);
    subplot(1,4,2)
       contour(w, w, imag(bd), 10), grid on 
       title('$imag(B)$ Burst 1', 'interpreter', 'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]);
    subplot(1,4,3)
       surf(wb1, wb1, bic1); view(2); shading interp
       %contour(wb1, wb1, bic1, 10), grid on 
       title('$b$', 'interpreter', 'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]); caxis([0 1]);
    subplot(1,4,4)
       surf(wb1, wb1, bic2); view(2); shading interp
       %contour(wb1, wb1, bic2, 10), grid on 
       title('$b^2$', 'interpreter','Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]);

    % AVERAGE --------------------------------------------------------------- %
    % ADAPTIVE
    figure;
    subplot(1,4,1)
       contour(rFA, rFA, abs(zBspecA), 10), grid on 
       title('mtmBispecd - $|\bar{B}|$', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]);
    subplot(1,4,2)
       contour(rFA, rFA, imag(zBspecA), 10), grid on 
       title('$imag(\bar{B})$', 'interpreter', 'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]);
    subplot(1,4,3)
       surf(rFA, rFA, rBic1A); view(2); shading interp
       %contour(rFA, rFA, rBic1A, 10), grid on 
       title('$\bar{b}$', 'interpreter', 'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]); caxis([0 1]);
    subplot(1,4,4)
       surf(rFA, rFA, rBic2A); view(2); shading interp
       %contour(rFA, rFA, rBic2A, 10), grid on 
       title('$\bar{b}^2$', 'interpreter', 'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]);

    % BISPECD
    figure;
    subplot(1,4,1)
       contour(wAve, wAve, abs(bdAve), 10), grid on 
       title('bispecd - $|\bar{B}|$', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]);
    subplot(1,4,2)
       contour(wAve, wAve, imag(bdAve), 10), grid on 
       title('bispecd - $imag(\bar{B})$', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]);
    subplot(1,4,3)
       surf(wbAve, wbAve, bic1Ave); view(2); shading interp
       %contour(wbAve, wbAve, bic1Ave, 10), grid on 
       title('$\bar{b}$', 'interpreter', 'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]); caxis([0 1]);
    subplot(1,4,4)
       surf(wbAve, wbAve, bic2Ave); view(2); shading interp
       %contour(wbAve, wbAve, bic2Ave, 10), grid on 
       title('$\bar{b}^2$', 'interpreter', 'Latex', 'FontSize', 14)
       xlabel('$f_1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f_2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       xlim([0 0.5]); ylim([0 0.5]); 

    % POWER SPECTRUMS 
    figure;
    subplot(1,2,1)
    loglog(rFA, snn(:,1), w(65:end), per(65:end))
    title('mean power spectra', 'interpreter', 'Latex', 'FontSize', 14)
    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
    ylabel('arb/Hz')
    legend('mtm-adaptive', 'periodogram')

    subplot(1,2,2)
    loglog(rFA, mean(snn,2), w(65:end), pyyAve(65:end))
    title('mean power spectra', 'interpreter', 'Latex', 'FontSize', 14)
    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
    ylabel('arb/Hz')
    legend('mtm-adaptive', 'periodogram')
end

%% white noise variance reduction 

if bWhiteNoise
    
    figure; 
    % create white gaussian data - iid (independent & identically distributed) 
    n = 64;
    rWNoise = randn(n, 1);
    % remove the mean and make unity variance
    rWNoise = (rWNoise - mean(rWNoise))/std(rWNoise); 
    var(rWNoise) % check that the variance is 1
    plot(rWNoise)
    % 
    % % raw biperiodogram, should have variance of 64 (because power spectrum
    % % should be) BUT not quite (more like high 50s)
    % [~,w,xf] = bispecd(rWNoise,128,0,n,0);
    % mask = hankel([1:128],[128,1:128-1]);
    % cxf = conj(xf);
    % bd = fftshift(xf * xf.' .* cxf(mask)) * length(rWNoise)^2;
    % figure;
    % contour(w, w, abs(bd), 10), grid on
    % rVarWN = mean(var(bd(3:120, 3:120))); % don't include edges b/c of edge effects
    % 
    % % check periodogram
    % per = fftshift(xf.*conj(xf)) * length(rWNoise); 
    % figure; loglog(w,per)

    % theoretical approximation 
    nNW = [1.5 2 2.5 3 3.5 4 4.5 5];
    rTVar = 3*nNW.^2;
    figure; 
    plot(nNW, rTVar, '-o')
    xlabel('NW', 'interpreter', 'Latex', 'FontSize', 16)
    ylabel('variance reduction factor [$var(B^{per})/var(\hat{B})$]', ...
           'interpreter', 'Latex', 'FontSize', 16)
    hold on 

    % load mtm white noise for each and run bispectrum, calculate variance
    fid = ["1pt5", "2", "2pt5", "3", "3pt5", "4", "4pt5", "5"];
    rVarMTM = NaN(size(fid));
    for i = 1 : length(fid)
        load(sprintf('wn_zp_NW%s.mat', fid(i))) 
        bPlot = false;
        [zBspecWN, ~, ~, rFWN] = mtmBispecd(Yk, vk, dk, f, snn, bPlot);
        %rVarMTM(i) = mean(var(zBspecWN(2:63, 2:63)));
        rVarMTM(i) = mean(var(zBspecWN));
    end
    plot( nNW, (n ./ rVarMTM) , '-o')
    legend('theoretical', 'MBE, K=2NW')
end

%% Confidence limits for MBE ----------------------------------------------
if bConfidence
    
    nSeg = 64;  % 64 segments
    nTS  = 512; % 512 data points
    rWNoise = zeros(nTS, nSeg);
    bPlot = false;
    
%     for i=1 : nSeg
%         % create white gaussian data - iid 
%         rWNoise(:,i) = randn(nTS, 1);
% 
%         % remove the mean and make unity variance
%         rWNoise(:,i) = (rWNoise(:,i) - mean(rWNoise(:,i)))/std(rWNoise(:,i)); 
%         
%         % debugging
%         var(rWNoise(:,i)) % check that the variance is 1
%         figure; plot(rWNoise(:,i))
%     end    
%        
%     dlmwrite('Confidence_GWN.txt',rWNoise,'delimiter','\t','precision',12)

    % run the adaptive case for all segments (ensemble average)  
    load('CL-GWN_NW3_K4.mat') 
    load('CL-GWN_NW3_K4_NA.mat') 
    NW = 3;
    K  = 4;
    
    % make symmetric filter for reducing size of bispec matrices
    nCol = size(f,1);
    rSym = NaN(nCol,nCol);
    cntRow = 1;
    for i = 1 : nCol
        rSym(1:cntRow,i) = 1;
        cntRow = cntRow + 1;
    end
    
    % run and filter
    bPlot = false;
    %[zBspecA, rBic1A, rBic2A, rFA] = mtmBispecd(Yk, vk, dk, f, snn, bPlot);
    [zBspecA, rBic1A, rBic2A, rFA] = mtmBispecdNA(Yk, vk, f, snn, bPlot);
    rBic1 = rBic1A .* rSym;
    rBic2 = rBic2A .* rSym;
    nanmean(rBic2(:))
              
    % check the goodness of fit for X^2 with 2 degrees of freedom (95%
    % confidence)
    scale = K/(3*(NW)^2*(2*K*nSeg));
    dof   = 2;
    alpha = dof/2;
    beta  = 2*scale;
%     [h,p] = chi2gof(rBic2(:),'cdf',@(x)chi2cdf(x,dof),'nparams',1);%, 'Alpha', 0.95);
%     %h % print the goodness of fit (0 if they come from the same distribution)
%     p % pvalue 

    % this is the theoretical distribution fit to the data
    pd    = makedist('Gamma','a', alpha, 'b', beta)
    %pd    = makedist('Exponential','mu', beta)
    figure; qqplot(rBic2(:), pd); 
  
    % this is the maximum likelihood estimate of the parameters of a gamma 
    % distribution fit to the data (with confidence limits)
    [phat,pci] = mle(rBic2(:),'distribution','gamma')
    
    % fit gamma distribution to bicoherence and plot pdf
    pd = fitdist(rBic2(:),'Gamma'); % same as MLE distribution
    %mean(pd) % this is the same as nanmean(rBic2)
    %std(pd)  % not quite the same but close
    figure; histfit(rBic2(:), 257, 'Gamma') % same gamma as above
    %figure; qqplot(rBic2(:), pd); 
    
    % preallocate arrays
    [rAMean, rAStd, nAdof] = deal(NaN(nSeg,1));
    
    % calculate the mean degrees of freedom for each segment
    nMdof = mean(dofs);
    
    for i=2 : nSeg
        % calculate the bispectrum for various degrees of freedom, achieved
        % by changing the number of segments over which you average
        %[zBspecA, rBic1A, rBic2A, rFA] = mtmBispecd(Yk(:,:,1:i), ...
        %    vk(:,:,1:i), dk(:,:,1:i), f(:,1:i), snn(:,1:i), bPlot);
        [zBspecA, rBic1A, rBic2A, rFA] = mtmBispecdNA(Yk(:,:,1:i), ...
            vk(:,:,1:i), f(:,1:i), snn(:,1:i), bPlot);
        
        % only keep non symmetric elements
        rBic1 = rBic1A .* rSym;
        rBic2 = rBic2A .* rSym;
    
        % find the mean and std 
        rAMean(i) = nanmean(rBic1(:));
        rAStd(i)  = nanstd(rBic1(:));
        
        % find the dof as dof (from Snn) * 3 segments
        nAdof(i) = sum(nMdof(1:i));
    end
    
        % check the 95% CL against theory
        nTdof    = 6:1:800;                % theoretical if X_2^2 distributed
        rTdof    = sqrt(6./nTdof);  
        rTmean   = sqrt(2./nTdof);
        %nMTMdof  = 2*K*(2:1:64);
        nMTMdof  = 2*K*(2:1:64);
        %nMTMmean = sqrt( (2*K) ./ (3*(NW)^2*nMTMdof) );
        nMTMmean = sqrt( (2*K) ./ (3*(NW)^2*nMTMdof) );
        %rMTMdof  = sqrt( (6*K) ./ (3*(NW)^2*nMTMdof) );
        rMTMdof  = sqrt( (1) ./ ((NW)^2*(2:1:64)) ); % same as above
        rMTMdof  = sqrt( (1) ./ ((K/2)^2*(2:1:64)) ); % same as above
        rCL      = rAMean + 2*rAStd;       % actual
        
        % plots
        figure; semilogx(nTdof,rTdof) 
        hold on
        semilogx(nMTMdof,rMTMdof)
        hold on
        semilogx(nAdof, rCL, '*')
        hold on
        %semilogx(nMTMdof, rCL(2:end), '*')
        xlabel('degrees of freedom', 'Interpreter', 'LaTex', 'FontSize', 16)
        ylabel('95$\%$ significance $\bar{b}$', 'Interpreter', 'LaTex', 'FontSize', 16)
        legend({'theoretical for $\chi_2^2$: $\sqrt{6/v}$', ...
               'MBE theoretical: $\sqrt{6K/3(NW)^2v}$', 'MBE'}, ...
               'Interpreter', 'LaTex', 'FontSize', 14)
        xlim([0 1000])
        
        % check the expectation
        figure; semilogx(nTdof, rTmean); hold on
        semilogx(nMTMdof, nMTMmean); hold on
        semilogx(nAdof, rAMean, '*')
        xlabel('degrees of freedom', 'Interpreter', 'LaTex', 'FontSize', 16)
        ylabel('$E[\bar{b}]$', 'Interpreter', 'LaTex', 'FontSize', 16)
        legend({'theoretical for $\chi_2^2$: $\sqrt{2/v}$', ...
               'MBE theoretical: $\sqrt{2K/3(NW)^2v}$', 'MBE'}, ...
               'Interpreter', 'LaTex', 'FontSize', 14)
        xlim([0 1000])
                
end
%% OLD PLOTS

% % REAL ------------------------------------------------------------------ %
% figure;
% % ADAPTIVE ----------------------
% % plot just the first realization 
%    subplot(3,2,1)
%    contour(rFA1, rFA1, abs(zBspecA1), 10), grid on 
%    title('mtmBispecd - $|{B}|$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot the average
%    subplot(3,2,2)
%    contour(rFA, rFA, abs(zBspecA), 10), grid on 
%    title('mtmBispecd - $|\bar{B}|$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % NON-ADAPTIVE -------------------  
% % plot just the first realization 
%    subplot(3,2,3)
%    contour(rFNA1, rFNA1, abs(zBspecNA1), 10), grid on 
%    title('mtmBispecdNA - $|{B}|$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot the average
%    subplot(3,2,4)
%    contour(rFNA, rFNA, abs(zBspecNA), 10), grid on 
%    title('mtmBispecdNA - $|\bar{B}|$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % BISPECD -------------------   
% % plot just the first realization
%    subplot(3,2,5)
%    contour(w, w, abs(bd), 10), grid on 
%    title('bispecd - $|B|$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot just the average
%    subplot(3,2,6)
%    contour(wAve, wAve, abs(bdAve), 10), grid on 
%    title('bispecd - $|\bar{B}|$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
%    
% % IMAGINARY ------------------------------------------------------------- %
% figure;
% % ADAPTIVE ------------------- 
% % plot just the first realization
%    subplot(3,2,1)
%    contour(rFA1, rFA1, imag(zBspecA1), 10), grid on 
%    title('mtmBispecd - $imag(B)$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot the average
%    subplot(3,2,2)
%    contour(rFA, rFA, imag(zBspecA), 10), grid on 
%    title('mtmBispecd - $imag(\bar{B})$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % NON-ADAPTIVE ------------------- 
% % plot just the first realization
%    subplot(3,2,3)
%    contour(rFNA1, rFNA1, imag(zBspecNA1), 10), grid on 
%    title('mtmBispecdNA - $imag(B)$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot the average
%    subplot(3,2,4)
%    contour(rFNA, rFNA, imag(zBspecNA), 10), grid on 
%    title('mtmBispecdNA - $imag(\bar{B})$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % BISPECD ------------------- 
% % plot just the first realization
%    subplot(3,2,5)
%    contour(w, w, imag(bd), 10), grid on 
%    title('bispecd - $imag(B)$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot just the average
%    subplot(3,2,6)
%    contour(wAve, wAve, imag(bdAve), 10), grid on 
%    title('bispecd - $imag(\bar{B})$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
%    
% % BICOHERENCE ----------------------------------------------------------- %
% figure;
% % ADAPTIVE ------------------- 
% % plot just the first realization
%    subplot(3,2,1)
%    contour(rFA1, rFA1, rBic1A1, 10), grid on 
%    title('mtmBispecd - $b$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot the average
%    subplot(3,2,2)
%    contour(rFA, rFA, rBic1A, 10), grid on 
%    title('mtmBispecd - $\bar{b}$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % NON-ADAPTIVE ------------------- 
% % plot just the first realization
%    subplot(3,2,3)
%    contour(rFNA1, rFNA1, rBic1NA1, 10), grid on 
%    title('mtmBispecdNA - $b$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot the average
%    subplot(3,2,4)
%    contour(rFNA, rFNA, rBic1NA, 10), grid on 
%    title('mtmBispecdNA - $\bar{b}$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % BISPECD -------------------    
% % plot just the first realization
%    subplot(3,2,5)
%    contour(wb1, wb1, bic1, 10), grid on 
%    title('bispecd - $b$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot just the average
%    subplot(3,2,6)
%    contour(wbAve, wbAve, bic1Ave, 10), grid on 
%    title('bispecd - $\bar{b}$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
%    
% % BICOHERENCE^2 --------------------------------------------------------- %
% figure;
% % ADAPTIVE -------------------  
% % plot just the first realization
%    subplot(3,2,1)
%    contour(rFA1, rFA1, rBic2A1, 10), grid on 
%    title('mtmBispecd - $b^2$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot the average
%    subplot(3,2,2)
%    contour(rFA, rFA, rBic2A, 10), grid on 
%    title('mtmBispecd - $\bar{b}^2$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % NON-ADAPTIVE -------------------  
% % plot just the first realization
%    subplot(3,2,3)
%    contour(rFNA1, rFNA1, rBic2NA1, 10), grid on 
%    title('mtmBispecdNA - $b^2$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot the average
%    subplot(3,2,4)
%    contour(rFNA, rFNA, rBic2NA, 10), grid on 
%    title('mtmBispecdNA - $\bar{b}^2$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % BISPECD -------------------    
% % plot just the first realization
%    subplot(3,2,5)
%    contour(wb1, wb1, bic2, 10), grid on 
%    title('bispecd - $b^2$ Burst 1', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
% % plot just the average
%    subplot(3,2,6)
%    contour(wbAve, wbAve, bic2Ave, 10), grid on 
%    title('bispecd - $\bar{b}^2$', 'interpreter', ...
%          'Latex', 'FontSize', 14)
%    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
%    colorbar
%    xlim([0 0.5]); ylim([0 0.5]);
%    
% % FULL LENGTH ----------------------------------------------------------- %
% % figure;
% % % ADAPTIVE ------------------- 
% % % real
% %    subplot(2,2,1)
% %    contour(rFAF, rFAF, abs(zBspecAF), 10), grid on 
% %    title('mtmBispecd - $|{B}|$ Burst 1', 'interpreter', ...
% %          'Latex', 'FontSize', 14)
% %    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
% %    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
% %    colorbar
% %    xlim([0 0.5]); ylim([0 0.5]);
% % % imaginary
% %    subplot(2,2,2)
% %    contour(rFAF, rFAF, imag(zBspecAF), 10), grid on 
% %    title('mtmBispecd - $imag(B)$ Burst 1', 'interpreter', ...
% %          'Latex', 'FontSize', 14)
% %    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
% %    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
% %    colorbar
% %    xlim([0 0.5]); ylim([0 0.5]);
% % % bicoherence 
% %    subplot(2,2,3)
% %    contour(rFAF, rFAF, rBic1AF, 10), grid on 
% %    title('mtmBispecd - $b$ Burst 1', 'interpreter', ...
% %          'Latex', 'FontSize', 14)
% %    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
% %    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
% %    colorbar
% %    xlim([0 0.5]); ylim([0 0.5]);
% % % bicoherence^2  
% %    subplot(2,2,4)
% %    contour(rFAF, rFAF, rBic2AF, 10), grid on 
% %    title('mtmBispecd - $b^2$ Burst 1', 'interpreter', ...
% %          'Latex', 'FontSize', 14)
% %    xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
% %    ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
% %    colorbar
% %    xlim([0 0.5]); ylim([0 0.5]);
%    
% % POWER SPECTRUMS (burst 1)
% figure;
% subplot(1,2,1)
% plot( 1:length(ts(:,1)), ts(:,1), 1:length(zmat(:,1)), zmat(:,1)) 
% title('time series - burst 1', 'interpreter', 'Latex', 'FontSize', 14)
% 
% subplot(1,2,2)
% loglog(rFA(:,1), snn(:,1), w(65:end), per(65:end))
% title('power spectrum - burst 1', 'interpreter', 'Latex', 'FontSize', 14)
% xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
% ylabel('arb/Hz')
% legend('mtm-adaptive', 'periodogram')