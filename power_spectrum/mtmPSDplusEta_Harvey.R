#mtmPSD_Harvey <- function(p, NW, K, dt) {
# -----------------------------mtmPSD_Harvey----------------------------
# Purpose: This function runs package multitaper to generate multitaper
# power spectra (single sided) from pressure data with confidence limits 
# and F-test. Requires multitaper & R.matlab packages. Has options for 
# producing Spp or Snn (WSE).
#
# Inputs:
#       - p:         evenly sampled array of pressure time series (Pa)
#       - dt:        sampling period (s)
#       - NW:        multitaper time-bandwidth product (typically 3-5)
#       - K:         number of Slepian tapers (NW*2-1 good starting pt)
#
# SEE ALSO: none
#
# Record of revisions:
#       Date            Programmer          Description of Change
#       =========================================================
#       11/18/18        K. Anarde           Condensed and organized
#
# ------------------------------user input------------------------------

# set working directory and load R packages
setwd("/Users/KatherineAnardeWheels/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM")
library(multitaper)
library(R.matlab)

# user inputs
etaIN <- read.table("ADVetaIN_68min_hpf_v2_0pt002.txt")  # pressure
dt    <- 1/16;        # sampling period (s)
NW    <- 8;           # MTM: time-bandwidth product
K     <- 15;          # MTM: number of tapers 

# ------------------------------preamble-------------------------------

# additional calculated input variables
N    <- length(etaIN[,1]); # number of samples 
b    <- NW/(N*dt);     # bandwidth
C    <- length(etaIN[1,]); # number of bursts
nFFT <- (N)/2+1;       # length of the fft (DON'T ZERO-PAD and include zero)

# compute the time vector for detrending
time  <- seq((1), N*dt+1, dt) #adds one data point...too tired to fix
time  <- time[1:N]
ttbar <- time - (time[N] + time[1])/2

# preallocate arrays
snn    <- matrix(data=NA, nrow=nFFT, ncol=C)
f      <- matrix(data=NA, nrow=nFFT, ncol=C)
dofs   <- matrix(data=NA, nrow=nFFT, ncol=C)
CIup   <- matrix(data=NA, nrow=nFFT, ncol=C)
CIlow  <- matrix(data=NA, nrow=nFFT, ncol=C)
fTest  <- matrix(data=NA, nrow=nFFT, ncol=C)

# ----------------------------generate MTM inputs--------------------------

for (i in 1:C){
	    
    ## calculate the power spectra for chosen NW and K; don't zero-pad and return eigencoefs
    psd <- spec.mtm(etaIN[,i], NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=FALSE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    Ftest=TRUE, jackknife=TRUE, centre=c("Slepian"))
                            
    ## plot spectra (for debugging and testing NW,K)
    #plot(dropFreqs(psd, 0, 0.04), log='yes', Ftest=FALSE, jackknife = TRUE)
    #plot(dropFreqs(psd, 0, 0.04), Ftest=TRUE, siglines=c(0.90, 0.99))

    ## save internals
    snn[,i]   <- psd$spec * 2     # multiply by two here (make single-sided)
    f[,i]     <- psd$freq         # actual frequencies
    dofs[,i]  <- psd$mtm$dofs     # degrees of freedom (2*K, except for zero/nyquist which have half)
    CIup[,i]  <- psd$mtm$jk$upperCI * 2 # multiply by two here (make single-sided)
    CIlow[,i] <- psd$mtm$jk$lowerCI * 2 # multiply by two here (make single-sided)
    fTest[,i] <- psd$mtm$Ftest
}

# write to .mat file
writeMat('OUT-mtmADVsp_nPlusHPF_v2_68minNW8K15.mat', dofs=dofs, fTest = fTest, CIlow = CIlow, CIup = CIup, f = f, snn = snn, 
    	fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)