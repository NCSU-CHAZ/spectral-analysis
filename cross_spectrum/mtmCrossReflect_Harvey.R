#mtmCrossReflect <- function(p, u, PTele, dt, NW, K) {
# -----------------------------mtmCrossReflect----------------------------
# Purpose: This function runs package multitaper to generate multitaper
# estimates of the power and co-spectrum for calculation of incoming
# and outgoing infragravity energy fluxes. Exports outputs to Matlab. 
#
# Inputs:
#       - p:         evenly sampled array of pressure time series (Pa)
#       - u:         " " of cross-shore velocity time series (m/s)
#       - PTele:     elevation of pressure transducer above the bed (m)
#       - dt:        sampling period (m)
#       - NW:        multitaper time-bandwidth product (typically 3-5)
#       - K:         number of Slepian tapers (NW*2-1 good starting pt)
#
# NOTE: to edit the source code in an R package interactively, use 
# e.g. trace("mtm.bispectrum", edit=TRUE), save, and then continue in editor
#
# Record of revisions:
#       Date            Programmer          Description of Change
#       =========================================================
#       09/28/18        K. Anarde           Original code
#
# ------------------------------user input------------------------------

# set working directory and load R packages
setwd("/Users/KatherineAnardeWheels/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM")
library(multitaper)
library(R.matlab)

# user inputs
p <- read.table("ADVsp_68min_hpf.txt")
u <- read.table("IN-VECTwinHPF_U_68min.txt")
PTele <- 0.7;  	      # PT elevation
dt    <- 1/16;        # sampling period (s)
NW    <- 8;           # MTM: time-bandwidth product
K     <- 15;          # MTM: number of tapers 

# ------------------------------preamble-------------------------------

# additional calculated input variables
N    <- length(p[,1]); # number of samples 
b    <- NW/(N*dt);     # bandwidth
C    <- length(p[1,]); # number of bursts
nFFT <- (N)/2+1;     # length of the fft (ZERO-PAD and include zero)

# convert pressure to water depth
rho  <- 1025;          # seawater density [kg/m3]
g    <- 9.81;          # gravitational constant [m/s2]
wse  <- p / (rho*g);   # water depth [m]

# compute the time vector for detrending
time  <- seq((1), N*dt+1, dt) #adds one data point...too tired to fix
time  <- time[1:N]
ttbar <- time - (time[N] + time[1])/2

# preallocate arrays
snn    <- matrix(data=NA, nrow=nFFT, ncol=C)
suu    <- matrix(data=NA, nrow=nFFT, ncol=C)
snnNA  <- matrix(data=NA, nrow=nFFT, ncol=C)
mscUN  <- matrix(data=NA, nrow=nFFT, ncol=C)
f      <- matrix(data=NA, nrow=nFFT, ncol=C)
h      <- matrix(data=NA, nrow=1, ncol=C)
Ykn    <- array(data=NA, dim = c(nFFT, K, C))
Yku    <- array(data=NA, dim = c(nFFT, K, C))
YknNA  <- array(data=NA, dim = c(nFFT, K, C))

# ----------------------------generate MTM inputs--------------------------

for (i in 1:C){
	    
    ## find a polynomial trend using the Slepian tapers and calculate the
    ## sea-surface residual (m)
    wsetrend <- multitaperTrend(wse[,i], B = b, deltat = dt, t.in = time)
    eta   <- wse[,i] - wsetrend[[1]] - wsetrend[[2]] * ttbar   
    h[,i] <- mean(wse[,i]+PTele)	  # mean water depth
    
    # find a polynomial trend using the Slepian tapers and calculate the 
    # velocity residual
	velXtrend <- multitaperTrend(u[,i], B = b, deltat = dt, t.in = time)
	velX.resid <- u[,i] - velXtrend[[1]] - velXtrend[[2]] * ttbar
    
    # compute power spectra and return internals (for mtm.coh to work, 
    # both spectra must have the same frequency resolution and return the zero 
    # frequency
    	
    SnnNA <- spec.mtm(eta, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=TRUE, adaptiveWeighting=FALSE, returnZeroFreq=TRUE, 
                    centre=c("Slepian"))
                    
    Suu <- spec.mtm(velX.resid, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=TRUE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    centre=c("Slepian"))
    
    Snn <- spec.mtm(eta, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=TRUE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    centre=c("Slepian"))   
                    
    # compute the magnitude squared coherence between eta and u
    Snu <- mtm.coh(Snn, Suu, phcorr=FALSE, plot=FALSE)     
    
    ## save internals for power spectra
    snn[,i]    <- Snn$spec    # double-sided spectra
    suu[,i]    <- Suu$spec    # double-sided spectra  
    snnNA[,i]  <- SnnNA$spec  # double-sided spectra 
    Ykn[,,i]   <- Snn$mtm$eigenCoefs # eigencoefficients (fourier transform of tapered data)
    Yku[,,i]   <- Suu$mtm$eigenCoefs 
    YknNA[,,i] <- SnnNA$mtm$eigenCoefs  
    f[,i]      <- Snn$freq      # actual frequencies
    mscUN[,i]  <- Snu$msc
}

# write outputs to a mat file
writeMat('IN-mtmCrossReflect_ADVnHPF_68min_NW8K15.mat', snn=snn, suu=suu, snnNA=snnNA, Ykn = Ykn, Yku = Yku, f = f, 
         mscUN = mscUN, YknNA=YknNA, h=h, fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)