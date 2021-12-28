#mtmBispecdIN <- function(p, u, PTele, dt, NW, K) {
# -----------------------------mtmBispecdIN----------------------------
# Purpose: This function runs package multitaper to generate multitaper
# parameters for bispectral analysis of sea-surface elevations (incoming
# wave signal) in Matlab. Also generates power spectra (double sided)
# with confidence limits and F-test. Requires multitaper & R.matlab 
# packages. By default, it only uses the incoming wave signal.
#
# Inputs:
#       - p:         evenly sampled array of pressure time series (Pa)
#       - u:         " " of cross-shore velocity time series (m/s)
#       - PTele:     elevation of pressure transducer above the bed (m)
#       - dt:        sampling period (s)
#       - NW:        multitaper time-bandwidth product (typically 3-5)
#       - K:         number of Slepian tapers (NW*2-1 good starting pt)
#
# SEE ALSO: mtmBispecd.m
#
# Record of revisions:
#       Date            Programmer          Description of Change
#       =========================================================
#       09/23/18        K. Anarde           Modified mtmBispecdExport.R
#                                           for incoming wave signal
#
# ------------------------------user input------------------------------

# set working directory and load R packages
setwd("/Users/KatherineAnardeWheels/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/WAVES/RRUs/INPUTS")
library(multitaper)
library(R.matlab)

# user inputs
#etaIN <- read.table("ADVetaIN_8min_hpf_v2_0pt002.txt") # eta, incoming only
etaIN <- read.table("RRU1_17min_eta_hpf_v2_0pt002.txt")  # just eta
dt    <- 1/16;        # sampling period (s)
NW    <- 5;           # MTM: time-bandwidth product (for 8 min=3, 17 min=5)
K     <- 9;           # MTM: number of tapers (for 8 min=5, 17 min=9)

# ------------------------------preamble-------------------------------

# additional calculated input variables
N    <- length(etaIN[,1]); # number of samples 
b    <- NW/(N*dt);     # bandwidth
C    <- length(etaIN[1,]); # number of bursts
nFFT <- (N)/2+1;     # length of the fft (don't ZERO-PAD and include zero)

# compute the time vector for detrending
time  <- seq((1), N*dt+1, dt) #adds one data point...too tired to fix
time  <- time[1:N]
ttbar <- time - (time[N] + time[1])/2

# preallocate arrays
snn    <- matrix(data=NA, nrow=nFFT, ncol=C)
f      <- matrix(data=NA, nrow=nFFT, ncol=C)
dofs   <- matrix(data=NA, nrow=nFFT, ncol=C)
vk     <- array(data=NA, dim = c(N, K, C))
dk     <- array(data=NA, dim = c(nFFT, K, C))
Yk     <- array(data=NA, dim = c(nFFT, K, C))
nFreq  <- matrix(data=NA, nrow=1, ncol=C)
nfft   <- matrix(data=NA, nrow=1, ncol=C)
CIup   <- matrix(data=NA, nrow=nFFT, ncol=C)
CIlow  <- matrix(data=NA, nrow=nFFT, ncol=C)
fTest  <- matrix(data=NA, nrow=nFFT, ncol=C)

# ----------------------------generate MTM inputs--------------------------

for (i in 1:C){
    
    ## calculate the power spectra for chosen NW and K; zero-pad and return eigencoefs
    res <- spec.mtm(etaIN[,i], NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=TRUE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    Ftest=TRUE, jackknife=TRUE, centre=c("Slepian"))
    
    ## plot spectra (for debugging and testing NW,K), only in infragravity range
    #plot(dropFreqs(res, 0, 0.20), log='yes', Ftest=FALSE, jackknife = TRUE)
    #plot(dropFreqs(res, 0, 0.25), Ftest=TRUE, siglines=c(0.90, 0.99))
    
    ## save internals for power spectra
    snn[,i]   <- res$spec      # don't multiply by two here (keep as a double-sided)
    f[,i]     <- res$freq      # actual frequencies
    dofs[,i]  <- res$mtm$dofs  # degrees of freedom (2*K, except for zero/nyquist which have half)
    CIup[,i]  <- res$mtm$jk$upperCI # don't multiply by two here (keep as a double-sided)
	 CIlow[,i] <- res$mtm$jk$lowerCI # don't multiply by two here (keep as a double-sided)
    fTest[,i] <- res$mtm$Ftest
    
    ## save internals 
    vk[ , , i] <- res$mtm$dpss$v       # slepian sequences
    dk[ , , i] <- res$mtm$eigenCoefWt  # eigenweights (here the sqrt(|dk|^2))
    Yk[ , , i] <- res$mtm$eigenCoefs   # eigencoefficients (fourier transform of tapered data)
    nFreq[1,i] <- res$mtm$nfreqs       # number of frequencies
    nfft[1,i]  <- res$mtm$nFFT         # size of nfft
}

# write outputs to a mat file
writeMat('IN-mtmBisp_RRU1_HPF_v2_17min_NW5K9.mat', snn=snn, f=f, dofs=dofs, vk = vk, Yk = Yk, nFreq = nFreq, 
         nfft = nfft, NW=NW, K=K, dk=dk, CIup=CIup, CIlow=CIlow, fTest=fTest, fixNames=TRUE, matVersion="5", 
         onWrite=NULL, verbose=FALSE)