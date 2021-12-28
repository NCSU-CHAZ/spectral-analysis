function [zBspec, rBic1, rBic2, rF] = mtmBispecd(yk, vk, dk, f, snn, bPlot) 
% -----------------------------mtmBispecd---------------------------------%
% Purpose: This function generates bispectra and calculates the bicoherence
% using outputs from Thomson's multitaper spectral estimator, package 
% "multitaper" in R (Rahim and Thomson, 2017). Here we use the direct FFT 
% approach following Thomson (1989) and Birkelund (200) and assume power 
% spectra were generated for \eta (water surface elevation) [m]. Two 
% methods of normalization for bicoherence are provided. Inputs can be for
% multiple realizations or for a single realization.
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
    rDk = dk(:,:,iRz);

    % ---- P(j,k,l) matrix ----    
    % calculate P(j,k,l) as in DJT's eq.(9) and gamma as in eq.(12) by 
    % summing over each discrete sample (n) in time series
    rP3d   = zeros(nK,nK,nK);
    nGamma = zeros(nFreq,nFreq);

    for j = 1 : nK
        for k = 1 : nK
            for l = 1 : nK
                for n = 1 : nLen 
                    rP3d(j,k,l) = rP3d(j,k,l) + ...
                                  ( rVk(n,j) * rVk(n,k) * rVk(n,l) );
                end

                % use Birkelund's normalization (U3') 
                nGamma = nGamma + ( (rDk(:,j) * rDk(:,k).') .* ...
                                    reshape(rDk(rMask,l), nFreq, nFreq) ...
                                    * rP3d(j,k,l)^2 );
            end
        end
    end

    % ---- accumulate triple products ----
    % calculate the expansion coefficients, DJT's \hat(x)(f) = dk(f)Yk(f)  
    zHatXf = rDk .* zYk;

    % find the complex conjugate of the expansion coefficients
    zCjHatXf = conj(zHatXf);
    %zCjYk   = conj(zYk);    % for check on Birkelund's derivation

    % accumulate the triple product of the expansion coefficients & norm
    % constant for each taper; calculate the triple sum via outer product 
    % (aa') of \hat(x)(f) and toggling of indices with the complex conj
    % using the hankel mask
    for j = 1 : nK
        for k = 1 : nK
            for l = 1 : nK
                % using Thomson's equation (9)
                zBspec(:,:,iRz) = zBspec(:,:,iRz) + ...
                           ( (zHatXf(:,j) * zHatXf(:,k).') .* ...
                             reshape(zCjHatXf(rMask,l), nFreq, nFreq) * ...
                             rP3d(j,k,l) );

                % using Birkelund's derivation (11 and 24) --> equivalent!
                %zBspec = zBspec + ( (zYk(:,j) * zYk(:,k).') .* ...
                %    reshape(zCjYk(rMask,l), nFreq, nFreq) ...
                %    .* ( (rDk(:,j) * rDk(:,k).') .* ...
                %    reshape(rDk(rMask,l), nFreq, nFreq) ) * rP3d(j,k,l) );
            end
        end
    end

    % normalize by gamma
    zBspec(:,:,iRz) = zBspec(:,:,iRz) ./ nGamma;
    
end

% mean bispectrum (separate real and imaginary)
rBspec = sum(real(zBspec),3) / nRz;
iBspec = sum(imag(zBspec),3) / nRz;
zBspec = complex(rBspec,iBspec);
%zBspec = mean(zBspec,3); % looks like this does the same thing as above

% mean power spectra
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

%% plot -------------------------------------------------------------------

if bPlot
    figure;
       subplot(2,2,1)
       %surf(rF, rF, abs(zBspec)); view(2); shading interp
       contour(rF, rF, abs(zBspec), 10), grid on 
       title('mtmBispecd - real $\hat{B}(f1,f2)$', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       %xlim([0 0.04]); ylim([0 0.04]);

       subplot(2,2,2)
       %surf(rF, rF, imag(zBspec)); view(2); shading interp
       contour(rF, rF, imag(zBspec), 10), grid on 
       title('mtmBispecd - imaginary $\hat{B}(f1,f2)$', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       %xlim([0 0.04]); ylim([0 0.04]);

       subplot(2,2,3)
       %surf(rF, rF, abs(rBic1)); view(2); shading interp
       contour(rF, rF, abs(rBic1), 10), grid on 
       title('mtmBispecd - $|\hat{b}(f1,f2)|$', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       %xlim([0 0.04]); ylim([0 0.04]);

       subplot(2,2,4)
       %surf(rF, rF, abs(rBic2)); view(2); shading interp
       contour(rF, rF, abs(rBic2).^2, 10), grid on 
       title('mtmBispecd - $|\hat{b}(f1,f2)|^2$', 'interpreter', ...
             'Latex', 'FontSize', 14)
       xlabel('$f1 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       ylabel('$f2 \ (Hz)$', 'interpreter', 'Latex', 'FontSize', 14)
       colorbar
       %xlim([0 0.04]); ylim([0 0.04]);
end
return