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
p     <- read.table("RRU1_34min_hpf_v2_0pt002.txt")  # pressure
dt    <- 1/16;       # sampling period (s)
NW    <- 5;         # MTM: time-bandwidth product
K     <- 9;         # MTM: number of tapers 
bSnn  <- TRUE;      # true for Snn (WSE PSD)

# ------------------------------preamble-------------------------------

# additional calculated input variables
N    <- length(p[,1]); # number of samples 
b    <- NW/(N*dt);     # bandwidth
C    <- length(p[1,]); # number of bursts
nFFT <- (N)/2+1;       # length of the fft (DON'T ZERO-PAD and include zero)

# compute the time vector for detrending
time  <- seq((1), N*dt+1, dt) #adds one data point...too tired to fix
time  <- time[1:N]
ttbar <- time - (time[N] + time[1])/2

# preallocate arrays
if (bSnn){
	snn    <- matrix(data=NA, nrow=nFFT, ncol=C)
} else{
	spp    <- matrix(data=NA, nrow=nFFT, ncol=C)
}
f      <- matrix(data=NA, nrow=nFFT, ncol=C)
dofs   <- matrix(data=NA, nrow=nFFT, ncol=C)
CIup   <- matrix(data=NA, nrow=nFFT, ncol=C)
CIlow  <- matrix(data=NA, nrow=nFFT, ncol=C)
fTest  <- matrix(data=NA, nrow=nFFT, ncol=C)

# ----------------------------generate MTM inputs--------------------------

for (i in 1:C){
	
	## amplitude function
    #asig <- HilbertTransform(p[,i])
    #env <- HilbertEnvelope(asig)
	    
    ## find a linear trend using the Slepian tapers and calculate the
    ## [sea-surface] residual 
    trend <- multitaperTrend(p[,i], B = b, deltat = dt, t.in = time)
    resid <- p[,i] - trend[[1]] - trend[[2]] * ttbar   
    
    # check that the detrended time series look right
    #plot(time, p[,i], type="l", xlab="time (s)", ylab = "pressure [Pa]")
    #lines(time, p[,i] - resid, type="l", col="red")
        
	# convert pressure to water depth (if user specified above)
	if (bSnn){
		rho  <- 1025;              # seawater density [kg/m3]
		g    <- 9.81;              # gravitational constant [m/s2]
		eta  <- resid / (rho*g);   # residual water surface elevation [m]
    
    	## calculate the power spectra for chosen NW and K; don't zero-pad and return eigencoefs
    	psd <- spec.mtm(eta, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=FALSE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    Ftest=TRUE, jackknife=TRUE, centre=c("Slepian"))
        } else {
      	## for some reason the program was triple padding here, so chose to default pad by 2
        psd <- spec.mtm(resid, NW, K, deltat=dt, plot=FALSE, nFFT=N, 
                    returnInternals=FALSE, adaptiveWeighting=TRUE, returnZeroFreq=TRUE, 
                    Ftest=TRUE, jackknife=TRUE, centre=c("Slepian"))	
        }
                    
    	## plot spectra (for debugging and testing NW,K)
    #	plot(dropFreqs(psd, 0, 0.04), log='yes', Ftest=FALSE, jackknife = TRUE)
    #	plot(dropFreqs(psd, 0, 0.04), Ftest=TRUE, siglines=c(0.90, 0.99))

    ## save internals
    if (bSnn){
    	snn[,i]   <- psd$spec * 2  # multiply by two here (make single-sided)
    } else {
    	spp[,i]   <- psd$spec * 2  # multiply by two here (make single-sided)
    }
    f[,i]     <- psd$freq      # actual frequencies
    dofs[,i]  <- psd$mtm$dofs  # degrees of freedom (2*K, except for zero/nyquist which have half)
    CIup[,i]  <- psd$mtm$jk$upperCI * 2 # multiply by two here (make single-sided)
    CIlow[,i] <- psd$mtm$jk$lowerCI * 2 # multiply by two here (make single-sided)
    fTest[,i] <- psd$mtm$Ftest
}

# write to .mat file
    if (bSnn){
    	writeMat('OUT-mtmRRU1_HPF_v2_34minNW5K9.mat', dofs=dofs, fTest = fTest, CIlow = CIlow, CIup = CIup, f = f, snn = snn, 
    	fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)
    } else {
    	writeMat('OUT-mtmG5_WG_68min.mat', dofs=dofs, fTest = fTest, CIlow = CIlow, CIup = CIup, f = f, snn = spp, 
    	fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)
    }