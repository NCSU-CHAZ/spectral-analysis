## preamble
# set working directory and load R packages
setwd("/Users/KatherineAnardeWheels/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM")
library(multitaper)
library(R.matlab)

data <- read.table("Confidence_GWN.txt")

# set variables
N  <- 512;       # number of samples
dt <- 1;         # period
NW <- 5;         # time-bandwidth
K  <- 10;        # tapers
b  <- NW/(N*dt)  # bandwidth
C  <- length(data[1,])
nFFT <- N/2+1;  # length of the nfft (don't need to zero pad, include zero)

# preallocate arrays
snn    <- matrix(data=NA, nrow=nFFT, ncol=C)
f      <- matrix(data=NA, nrow=nFFT, ncol=C)
dofs   <- matrix(data=NA, nrow=nFFT, ncol=C)
vk     <- array(data=NA, dim = c(N, K, C))
dk     <- array(data=NA, dim = c(nFFT, K, C))
Yk     <- array(data=NA, dim = c(nFFT, K, C))
nFreq  <- matrix(data=NA, nrow=1, ncol=C)
nfft   <- matrix(data=NA, nrow=1, ncol=C)

## adaptive (mean already removed)
for (i in 1:C){
    
    ## calculate the power spectra for chosen NW and K; zero-pad and return eigencoefs
    #res <- spec.mtm(data[,i], NW, K, deltat=dt, plot=FALSE, nFFT=N, returnInternals=TRUE, adaptiveWeighting=FALSE, returnZeroFreq=TRUE)
    res <- spec.mtm(data[,i], NW, K, deltat=dt, plot=FALSE, nFFT=N, returnInternals=TRUE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE)
    
    #plot(res, log='yes', Ftest=FALSE, jackknife = FALSE)
    
    ## save internals for power spectra
    snn[,i]   <- res$spec                # don't multiply by two here (keep as a double-sided power spectra)
    f[,i]     <- res$freq                # actual frequencies
    dofs[,i]  <- res$mtm$dofs            # degrees of freedom (2*K, except for zero and nyquist which have half)
    
    ## save internals for bispectrum/bicoherence
    vk[ , , i] <- res$mtm$dpss$v       # slepian sequences
    dk[ , , i] <- res$mtm$eigenCoefWt  # eigenweights (here the sqrt(|dk|^2))
    Yk[ , , i] <- res$mtm$eigenCoefs   # eigencoefficients (fourier transform of tapered data)
    nFreq[1,i] <- res$mtm$nfreqs       # number of frequencies
    nfft[1,i]  <- res$mtm$nFFT         # size of nfft
}
writeMat('CL-GWN_NW5_K10.mat', snn=snn, f = f, vk = vk, Yk = Yk, dk=dk, dofs=dofs, fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)
#writeMat('CL-GWN_NW5_K10_NA.mat', snn=snn, f = f, vk = vk, Yk = Yk, dofs=dofs, fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)