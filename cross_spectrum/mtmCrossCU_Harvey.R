#mtmCrossReflect <- function(c, u, v, dt, NW, K) {
# -----------------------------mtmCrossCU----------------------------
# Purpose: This function runs package multitaper to generate multitaper
# estimates of the power and co-spectrum for suspended sediment conce.
# and cross-shore infragravity velocity. We include the alongshore velocity
# to check the magniture of ross-shore vs alongshore velocity variance 
# (rule of thumb, alongshore >75% cross-shore discarded). Exports 
# outputs to Matlab. 
#
# Inputs:
#       - c:         evenly sampled array of sediment concentration (kg/m^3)
#       - u:         " " of cross-shore velocity time series (m/s)
#       - v:         " " of alongshore velocity time series (m/s)
#       - dt:        sampling period (m)
#       - NW:        multitaper time-bandwidth product (typically 3-5)
#       - K:         number of Slepian tapers (NW*2-1 good starting pt)
#
# Record of revisions:
#       Date            Programmer          Description of Change
#       =========================================================
#       03/25/19        K. Anarde           New code
#
# ------------------------------user input------------------------------

# set working directory and load R packages
setwd("/Users/KatherineAnardeWheels/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/ANALYSIS/VELOCITY/INPUTS")
library(multitaper)
library(R.matlab)

# user inputs
ss <- read.table("IN-VECTwin_SL_137min.txt")
u <- read.table("IN-VECTwin_U_137min.txt")
v <- read.table("IN-VECTwin_V_137min.txt")
dt    <- 1/16;        # sampling period (s)
NW    <- 3;           # MTM: time-bandwidth product
K     <- 4;           # MTM: number of tapers 

# ------------------------------preamble-------------------------------

# additional calculated input variables
N    <- length(ss[,1]); # number of samples 
b    <- NW/(N*dt);     # bandwidth
C    <- length(ss[1,]); # number of bursts
nFFT <- (N)/2+1;       # length of the fft (DON'T ZERO-PAD and include zero)

# compute the time vector for detrending
time  <- seq((1), N*dt+1, dt) #adds one data point...too tired to fix
time  <- time[1:N]
ttbar <- time - (time[N] + time[1])/2

# preallocate arrays
scc    <- matrix(data=NA, nrow=nFFT, ncol=C)
suu    <- matrix(data=NA, nrow=nFFT, ncol=C)
svv    <- matrix(data=NA, nrow=nFFT, ncol=C) 
mscCU  <- matrix(data=NA, nrow=nFFT, ncol=C)
mscCV  <- matrix(data=NA, nrow=nFFT, ncol=C)
f      <- matrix(data=NA, nrow=nFFT, ncol=C)
Ykc    <- array(data=NA, dim = c(nFFT, K, C))
Yku    <- array(data=NA, dim = c(nFFT, K, C))
Ykv    <- array(data=NA, dim = c(nFFT, K, C))

# ----------------------------generate MTM inputs--------------------------

for (i in 1:C){
	    
    ## find a polynomial trend using the Slepian tapers and calculate the
    ## residuals
	velXtrend <- multitaperTrend(u[,i], B = b, deltat = dt, t.in = time)
	velX.resid <- u[,i] - velXtrend[[1]] - velXtrend[[2]] * ttbar
	
    velYtrend <- multitaperTrend(v[,i], B = b, deltat = dt, t.in = time)
	velY.resid <- v[,i] - velYtrend[[1]] - velYtrend[[2]] * ttbar
	
	ctrend <- multitaperTrend(ss[,i], B = b, deltat = dt, t.in = time)
	c.resid <- ss[,i] - ctrend[[1]] - ctrend[[2]] * ttbar
    
    # compute power spectra and return internals (for mtm.coh to work, 
    # both spectra must have the same frequency resolution and return the zero 
    # frequency
                    
    Suu <- spec.mtm(velX.resid, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=TRUE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    centre=c("Slepian"))
    
    Svv <- spec.mtm(velY.resid, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=TRUE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    centre=c("Slepian")) 
                    
    Scc <- spec.mtm(c.resid, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=TRUE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    centre=c("Slepian"))                   
                    
    # compute the magnitude squared coherence between c and u
    Scu <- mtm.coh(Scc, Suu, phcorr=FALSE, plot=FALSE)  
    Scv <- mtm.coh(Scc, Svv, phcorr=FALSE, plot=FALSE)     
    
    ## save internals for power spectra
    scc[,i]    <- Scc$spec    # double-sided spectra
    suu[,i]    <- Suu$spec    # double-sided spectra 
    svv[,i]    <- Svv$spec    # double-sided spectra 
    Ykc[,,i]   <- Scc$mtm$eigenCoefs # eigencoefficients (fourier transform of tapered data)
    Yku[,,i]   <- Suu$mtm$eigenCoefs 
    Ykv[,,i]   <- Svv$mtm$eigenCoefs 
    f[,i]      <- Scc$freq    # actual frequencies
    mscCU[,i]  <- Scu$msc
    mscCV[,i]  <- Scv$msc    
}

# write outputs to a mat file
writeMat('IN-mtmCrossSL-U_ADV_137min.mat', scc=scc, suu=suu, svv=svv, Ykc = Ykc, Yku = Yku, Ykv = Ykv, f = f, mscCU = mscCU, mscCV = mscCV, fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)