#mtmCrossADVRRU6 <- function(p1, p2, u1, dt, NW, K) {
# --------------------------mtmCrossADVRRU6-----------------------------
# Purpose: This function runs package multitaper to generate multitaper
# estimates of the power and cross-spectrum between atmospheric pressure
# and seapressure. Exports outputs to Matlab. 
#
# Inputs:
#       - p1:        evenly sampled array of ADV seapressure (Pa)
#       - p2:        " " of RRU6 seapressure (Pa)
#       - dt:        sampling period (s)
#       - NW:        multitaper time-bandwidth product (typically 3-5)
#       - K:         number of Slepian tapers (NW*2-1 good starting pt)
#
# SEE ALSO: mtmCrossReflect_Harvey.R
#
# Record of revisions:
#       Date            Programmer          Description of Change
#       =========================================================
#       12/4/18         K. Anarde           Modified from mtmCrossVel_Harvey
#											for n+ at ADV and upsampling
#
# ------------------------------user input------------------------------

# set working directory and load R packages
setwd("/Users/KatherineAnardeWheels/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM")
library(multitaper)
library(R.matlab)

# user inputs
spADV <- read.table("ADVsp_128min_psin.txt") # ADV
spRRU <- read.table("RRU6_137min_psin_16Hz.txt") # RRU6
uADV  <- read.table("IN-VECTwin_U_137min.txt")   # ADV cross-shore velocity
PTele <- 0.7;  	      # ADV PT elevation
dt    <- 1/16;        # sampling period (s)
NW    <- 3;           # MTM: time-bandwidth product
K     <- 4;           # MTM: number of tapers 

# ------------------------------preamble-------------------------------

# additional calculated input variables
N    <- length(spADV[,1]); # number of samples 
b    <- NW/(N*dt);         # bandwidth
C    <- length(spADV[1,]); # number of bursts
nFFT <- (N)/2+1;     # length of the fft (DON'T ZERO-PAD and include zero)

# convert pressure to water depth
rho  <- 1025;          # seawater density [kg/m3]
g    <- 9.81;          # gravitational constant [m/s2]
wseADV  <- spADV / (rho*g)
wseRRU  <- spRRU / (rho*g)

# compute the time vector for detrending
time  <- seq((1), N*dt+1, dt) #adds one data point...too tired to fix
time  <- time[1:N]
ttbar <- time - (time[N] + time[1])/2

# preallocate arrays
snnADV    <- matrix(data=NA, nrow=nFFT, ncol=C)
snnRRU    <- matrix(data=NA, nrow=nFFT, ncol=C)
f         <- matrix(data=NA, nrow=nFFT, ncol=C)
CIupADV   <- matrix(data=NA, nrow=nFFT, ncol=C)
CIlowADV  <- matrix(data=NA, nrow=nFFT, ncol=C)
CIupRRU   <- matrix(data=NA, nrow=nFFT, ncol=C)
CIlowRRU  <- matrix(data=NA, nrow=nFFT, ncol=C)
ph    <- matrix(data=NA, nrow=nFFT, ncol=C)
msc   <- matrix(data=NA, nrow=nFFT, ncol=C)
nRRU   <- matrix(data=NA, nrow=N, ncol=C)
nADVin <- matrix(data=NA, nrow=N, ncol=C)
ntmsc <- matrix(data=NA, nrow=nFFT, ncol=C)
phVar <- matrix(data=NA, nrow=nFFT, ncol=C)
ntVar <- matrix(data=NA, nrow=nFFT, ncol=C)
ntSig <- matrix(data=NA, nrow=3, ncol=C)

# ----------------------------generate MSC --------------------------

# compute magnitude squared coherence between ADV and RRU wse over hour-long bursts
for (i in 1:C){

    ## find a linear trend using the Slepian tapers and calculate the
    ## sea-surface residual (m)
    trend <- multitaperTrend(wseADV[,i], B = b, deltat = dt, t.in = time)
    etaADV<- wseADV[,i] - trend[[1]] - trend[[2]] * ttbar   
    
    trend <- multitaperTrend(wseRRU[,i], B = b, deltat = dt, t.in = time)
    etaRRU<- wseRRU[,i] - trend[[1]] - trend[[2]] * ttbar   
    
    ## separate ADV eta to incoming component only
    h     <- (wseADV[,i] - etaADV) + PTele	 # hydrostatic water depth (not mean)
	etaIN <- (etaADV + ( sqrt(h/g) * uADV[,i] )) * (1/2)
    etaOUT<- (etaADV - ( sqrt(h/g) * uADV[,i] )) * (1/2)
    
    ## check that the incoming and outgoing wave signals look correct
	#plot(time, etaADV, type="l", xlab="time (s)", ylab = "water depth (m)")
    #lines(time, etaIN-mean(etaIN), type="l", col="red")
    #lines(time, etaOUT-mean(etaOUT), type="l", col="green")
    
    ## check that the incoming and outgoing wave signals look correct
	#plot(time, etaRRU, type="l", xlab="time (s)", ylab = "water depth (m)")
    
    # compute power spectra and return internals (for mtm.coh to work, 
    # both spectra must have the same frequency resolution and return the zero 
    # frequency) 
    adv.psd <- spec.mtm(etaIN, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=TRUE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    jackknife = TRUE, centre=c("Slepian"))#, Ftest=TRUE)
                    
    rru.psd <- spec.mtm(etaRRU, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=TRUE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    jackknife = TRUE, centre=c("Slepian"))#, Ftest=TRUE)

    # compute the magnitude squared coherence
    coh <- mtm.coh(adv.psd, rru.psd, phcorr=FALSE, plot=FALSE)
    
    # plot coherence and save to get confidence intervals for uncorrelated data from the transformed msc 	
    # (nehlim=10, nehc=4 for smoothing)
    pcoh <- plot(dropFreqs(coh, 0, 0.01), percentGreater=NULL, nehlim=0, nehc=0, cdfQuantilesTicks=NULL, 	
    			 drawPercentLines=TRUE, percentG=c(.5,.95,.99))
    
         # save variables to matrix
         snnADV[,i]   <- adv.psd$spec * 2  # have to multiply by two for one-sided spectra (normalization)
         snnRRU[,i]   <- rru.psd$spec * 2  # have to multiply by two for one-sided spectra (normalization)
         f[,i]        <- adv.psd$freq      # same for all
         CIupADV[,i]  <- adv.psd$mtm$jk$upperCI * 2
         CIlowADV[,i] <- adv.psd$mtm$jk$lowerCI * 2
         CIupRRU[,i]  <- rru.psd$mtm$jk$upperCI * 2
         CIlowRRU[,i] <- rru.psd$mtm$jk$lowerCI * 2
         ph[,i]       <- coh$ph
         msc[,i]      <- coh$msc
         nADVin[,i]   <- etaIN
         nRRU[,i]     <- etaRRU
         ntmsc[,i]    <- coh$NTmsc
         phVar[,i]    <- coh$phvar
         ntVar[,i]    <- coh$NTvar
         ntSig[,i]    <- pcoh$sigNT
  }

# write to .mat file
writeMat('OUT-mtmADVvRRU6_16Hz_nPlus_137min.mat', snnADV=snnADV, snnRRU=snnRRU, f=f, CIupADV=CIupADV, CIupRRU=CIupRRU, CIlowADV=CIlowADV, CIlowRRU=CIlowRRU, ph=ph, msc=msc, nADVin=nADVin, nRRU=nRRU, ntmsc=ntmsc, phVar = phVar, ntVar = ntVar, ntSig = ntSig, fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)