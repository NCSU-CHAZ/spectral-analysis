#mtmCrossADVRRU6 <- function(p1, p2, u1, dt, NW, K) {
# --------------------------mtmCrossATMnPlus----------------------------
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
spATM <- read.table("ATMsub_128min_16Hz.txt")     # ATM
#tmpATM <- c(spATM[65536:131072,11], spATM[1:65536,12])
spADV <- read.table("ADVsp_128min_psin.txt")   # ADV
#tmpADV <- c(spADV[65536:131072,11], spADV[1:65536,12])
uADV  <- read.table("IN-VECTwin_U_137min.txt")    # ADV cross-shore velocity
#tmpU  <- c(uADV[65536:131072,11], uADV[1:65536,12])
PTele <- 0.7;  	      # ADV PT elevation
dt    <- 1/16;        # sampling period (s)
NW    <- 3;           # MTM: time-bandwidth product
K     <- 4;           # MTM: number of tapers 

# ------------------------------preamble-------------------------------

# additional calculated input variables
N    <- length(spADV[,1]); # number of samples 
#N <- length(tmpADV)
b    <- NW/(N*dt);         # bandwidth
C    <- length(spADV[1,]); # number of bursts
nFFT <- (N)/2+1;     # length of the fft (DON'T ZERO-PAD and include zero)

# convert pressure to water depth
rho  <- 1025;          # seawater density [kg/m3]
g    <- 9.81;          # gravitational constant [m/s2]
wseADV  <- spADV / (rho*g)
#wseADV  <- tmpADV / (rho*g)

# compute the time vector for detrending
time  <- seq((1), N*dt+1, dt) #adds one data point...too tired to fix
time  <- time[1:N]
ttbar <- time - (time[N] + time[1])/2

# preallocate arrays
snnADV    <- matrix(data=NA, nrow=nFFT, ncol=C)
sppATM    <- matrix(data=NA, nrow=nFFT, ncol=C)
f         <- matrix(data=NA, nrow=nFFT, ncol=C)
CIupADV   <- matrix(data=NA, nrow=nFFT, ncol=C)
CIlowADV  <- matrix(data=NA, nrow=nFFT, ncol=C)
CIupATM   <- matrix(data=NA, nrow=nFFT, ncol=C)
CIlowATM  <- matrix(data=NA, nrow=nFFT, ncol=C)
ph    <- matrix(data=NA, nrow=nFFT, ncol=C)
msc   <- matrix(data=NA, nrow=nFFT, ncol=C)
pATM  <- matrix(data=NA, nrow=N, ncol=C)
nADVin <- matrix(data=NA, nrow=N, ncol=C)
ntmsc <- matrix(data=NA, nrow=nFFT, ncol=C)
phVar <- matrix(data=NA, nrow=nFFT, ncol=C)
ntVar <- matrix(data=NA, nrow=nFFT, ncol=C)
ntSig <- matrix(data=NA, nrow=3, ncol=C)

# ----------------------------generate MSC --------------------------

# compute magnitude squared coherence between ADV and RRU wse over hour-long bursts (15-20, 30)
for (i in 1:5){

    ## find a linear trend using the Slepian tapers and calculate the
    ## sea-surface residual (m)
    trend <- multitaperTrend(wseADV[,i], B = b, deltat = dt, t.in = time)
    etaADV<- wseADV[,i] - trend[[1]] - trend[[2]] * ttbar   
    #trend <- multitaperTrend(wseADV, B = b, deltat = dt, t.in = time)
    #etaADV<- wseADV - trend[[1]] - trend[[2]] * ttbar 
    
    trend <- multitaperTrend(spATM[,i], B = b, deltat = dt, t.in = time)
    resid <- spATM[,i] - trend[[1]] - trend[[2]] * ttbar   
    #trend <- multitaperTrend(tmpATM, B = b, deltat = dt, t.in = time)
	#resid <- tmpATM - trend[[1]] - trend[[2]] * ttbar 
    
    ## separate ADV eta to incoming component only
    h     <- (wseADV[,i] - etaADV) + PTele	 # hydrostatic water depth (not mean)
	etaIN <- (etaADV + ( sqrt(h/g) * uADV[,i] )) * (1/2)
    etaOUT<- (etaADV - ( sqrt(h/g) * uADV[,i] )) * (1/2)
    #h     <- (wseADV - etaADV) + PTele	 # hydrostatic water depth (not mean)
	#etaIN <- (etaADV + ( sqrt(h/g) * tmpU )) * (1/2)
    #etaOUT<- (etaADV - ( sqrt(h/g) * tmpU )) * (1/2)
    
    ## check that the incoming and outgoing wave signals look correct
	plot(time, etaADV, type="l", xlab="time (s)", ylab = "water depth (m)")
    lines(time, etaIN-mean(etaIN), type="l", col="red")
    #lines(time, etaOUT-mean(etaOUT), type="l", col="green")
    
    ## check that the incoming and outgoing wave signals look correct
	#plot(time, resid, type="l", xlab="time (s)", ylab = "atmospheric pressure (Pa)")
    
    # compute power spectra and return internals (for mtm.coh to work, 
    # both spectra must have the same frequency resolution and return the zero 
    # frequency) 
    adv.psd <- spec.mtm(etaIN, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=TRUE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    jackknife = TRUE, centre=c("Slepian"))
                    
    atm.psd <- spec.mtm(resid, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=TRUE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    jackknife = TRUE, centre=c("Slepian"))
                    
    # plot the PSDs
    #plot(dropFreqs(atm.psd, 0, 0.01))
    #plot(dropFreqs(adv.psd, 0, 0.01))

    # compute the magnitude squared coherence
    coh <- mtm.coh(atm.psd, adv.psd, phcorr=FALSE, plot=FALSE)
    
    # plot coherence and save to get confidence intervals for uncorrelated data from the transformed msc 	
    # (nehlim=10, nehc=4 for smoothing)
    pcoh <- plot(dropFreqs(coh, 0, 0.01), percentGreater=NULL, nehlim=0, nehc=0, cdfQuantilesTicks=NULL, 	
    			 drawPercentLines=TRUE,percentG=c(.5,.95,.99))
    
         # save variables to matrix
         snnADV[,i]   <- adv.psd$spec * 2  # have to multiply by two for one-sided spectra (normalization)
         sppATM[,i]   <- atm.psd$spec * 2  # have to multiply by two for one-sided spectra (normalization)
         f[,i]        <- adv.psd$freq      # same for all
         CIupADV[,i]  <- adv.psd$mtm$jk$upperCI * 2
         CIlowADV[,i] <- adv.psd$mtm$jk$lowerCI * 2
         CIupATM[,i]  <- atm.psd$mtm$jk$upperCI * 2
         CIlowATM[,i] <- atm.psd$mtm$jk$lowerCI * 2
         ph[,i]       <- coh$ph
         msc[,i]      <- coh$msc
         nADVin[,i]   <- etaIN
         pATM[,i]     <- resid
         ntmsc[,i]    <- coh$NTmsc
         phVar[,i]    <- coh$phvar
         ntVar[,i]    <- coh$NTvar
         ntSig[1:3,i]    <- pcoh$sigNT
  }

# write to .mat file
writeMat('OUT-mtmCrossATMnPLUS_137minUS.mat', snnADV=snnADV, sppATM=sppATM, f=f, CIupADV=CIupADV, CIupATM=CIupATM, CIlowADV=CIlowADV, CIlowATM=CIlowATM, ph=ph, msc=msc, nADVin=nADVin, pATM=pATM, ntmsc=ntmsc, phVar = phVar, ntVar = ntVar, ntSig = ntSig, fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)