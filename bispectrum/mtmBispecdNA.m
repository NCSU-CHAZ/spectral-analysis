function [zBspec, rBic1, rBic2, rF] = mtmBispecdNA(yk, vk, f, snn, bPlot) 
% -----------------------------mtmBispecdNA-------------------------------%
% Purpose: This function generates bispectra and calculates the bicoherence
% using outputs from Thomson's multitaper spectral estimator (with NO 
% adaptive weighting procedure), package "multitaper" in R (Rahim and 
% Thomson, 2017). Here we use the direct FFT approach following Thomson 
% (1989) and Birkelund (200) and assume power spectra were generated for 
% eta (water surface elevation) [m]. Two methods of normalization for
% bicoherence are provided. Inputs can be for multiple realizations or for 
% a single realization.
% 
% Inputs:
%       - yk:        eigencoefficients (complex, S(f)=|Yk(f)|^2), [m/Hz]
%       - vk:        Slepian sequences
%       - dk:        adaptive weighting coefficients
%       - f:         frequencies from fft
%       - snn:       double-sided power spectral density [m^2/Hz]
%       - bPlot:     boolean for plotting
%
% Outputs:
%       - zBspec:    bispectrum for the upper right quadrant (positive)
%       - rBic1:     bicoherence using the Haubrich (1965) normalization
%                    [b(f1,f2)] - upper bound of 1
%       - rBic2:     skewness function (often called bicoherence) 
%                    [b^2(f1,f2)] - no upper bound 
%       - rF:        vector of frequencies 
%
% SEE ALSO: bispecd from HOSA package, multitaper (R)
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       03/11/18        KAnarde             Original code
%
%% preamble ---------------------------------------------------------------

nLen  = size(vk,1);       % length of input time series
nFreq = size(f,1);        % number of frequencies
nK    = size(vk,2);       % number of tapers
nRz   = size(snn,2);      % number of realizations
rF    = f(:,1);           % only need one frequency vector

% bispectrum will be an nFreq x nFreq x nReal matrix (doesn't incl neg f)
zBspec = zeros(nFreq,nFreq,nRz);

% hankel mask quickens the triple products summation; returns a matrix
% whose first column is c and last r [hankel(c,r)]
rMask = hankel( (1:nFreq), [nFreq,1:nFreq-1]); 

%% bispectrum -------------------------------------------------------------
% calculate the bispectrum for each realization
for iRz = 1 : nRz
    
    % set local variables
    zYk = yk(:,:,iRz);
    rVk = vk(:,:,iRz);

    % ---- P(j,k,l) matrix ----    
    % calculate P(j,k,l) as in DJT's eq.(9) and gamma as in eq.(12) by 
    % summing over each discrete sample (n) in time series
    rP3d   = zeros(nK,nK,nK);
    nGamma = 0;

    for j = 1 : nK
        for k = 1 : nK
            for l = 1 : nK
                for n = 1 : nLen 
                    rP3d(j,k,l) = rP3d(j,k,l) + (rVk(n,j)*rVk(n,k)*rVk(n,l));
                end
                    % here Thomson's gamma = Birkelund's U3
                    nGamma = nGamma + rP3d(j,k,l)^2;
            end
        end
    end

    % ---- accumulate triple products ----
    % find the complex conjugate of the expansion coefficients
    zCjYk = conj(zYk); % for check on Birkelund's derivation

    % accumulate the triple product of the expansion coefficients & norm
    % constant for each taper; calculate the triple sum via outer product 
    % (aa') of \hat(x)(f) and toggling of indices with the complex conj
    % using the hankel mask
    for j = 1 : nK
        for k = 1 : nK
            for l = 1 : nK
                % using Birkelund's derivation (11), no adaptive weighting
                zBspec(:,:,iRz) = zBspec(:,:,iRz) + ...
                                ( (zYk(:,j) * zYk(:,k).') .* ...
                                  reshape(zCjYk(rMask,l),nFreq,nFreq) * ...
                                  rP3d(j,k,l) );
            end
        end
    end

    % normalize by gamma
    zBspec(:,:,iRz) = zBspec(:,:,iRz) / nGamma;
    
end

% sum the bispectra for each realization and average
rBspec = sum(real(zBspec),3) / nRz;
iBspec = sum(imag(zBspec),3) / nRz;
zBspec = complex(rBspec,iBspec);
%zBspec = mean(zBspec,3); % looks like this does the same thing as above

% same for the power spectra
rSnn = mean(snn,2);

%% bicoherence ------------------------------------------------------------
% convert process power spectrum and bispectrum to moments
rM2  = rSnn / nLen;      % S(f) = N*M2(f)
zM3  = zBspec / nLen^2;  % B(f1,f2) = N^2*M3(f1,f2) 

% Haubrich normalization which has no upper bound
% [see Birkelund et al. (2001) and Haubrich (1965)]                   
rBic1 = abs(zM3) ./ sqrt(rM2 * rM2.' .* rM2(rMask));  

% skewness function [Collis et al. (1998)]
rBic2 = abs(zM3).^2 ./ (rM2 * rM2.' .* rM2(rMask)); 

%rBic = abs(zBspec).^2 ./ (rSnn * rSnn.' .* rSnn(rMask)); % debugging 

%% plot -------------------------------------------------------------------

if bPlot
    figure;
       subplot(2,2,1)
       contour(rF, rF, real(zBspec), 10), grid on 
       title('mtmBispecdNA - real $\hat{B}(f1,f2)$', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       %xlim([0 0.04]); ylim([0 0.04]);

       subplot(2,2,2)
       contour(rF, rF, imag(zBspec), 10), grid on 
       title('mtmBispecdNA - imaginary $\hat{B}(f1,f2)$', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       %xlim([0 0.04]); ylim([0 0.04]);

       subplot(2,2,3)
       %surf(rF, rF, abs(rBic1)); view(2); shading interp
       contour(rF, rF, abs(rBic1), 10), grid on 
       title('mtmBispecdNA - $|\hat{b}(f1,f2)|$', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       %xlim([0 0.04]); ylim([0 0.04]);

       subplot(2,2,4)
       %surf(rF, rF, abs(rBic2)); view(2); shading interp
       contour(rF, rF, abs(rBic2).^2, 10), grid on 
       title('mtmBispecdNA - $|\hat{b}(f1,f2)|^2$', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       %xlim([0 0.04]); ylim([0 0.04]);
end
return