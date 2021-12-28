## preamble
setwd("/Users/KatherineAnardeWheels/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM")

data <- read.table("WhiteNoise.txt")
N <- 64
NW <- c(1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
dt <- 1
K <- c(3, 4, 5, 6, 7, 8, 9, 10)
nFFT <- (N*2)/2+1;

## adaptive (mean already removed)
for (i in 1:8){

	# zero-pad
	res <- spec.mtm(data, NW[i], K[i], deltat=1, plot=FALSE, nFFT="default", returnInternals=TRUE,
					adaptiveWeighting=TRUE, returnZeroFreq=TRUE)
	snn <- res$spec             # don't multiply by two here (keep as a double-sided power spectra)
	f   <- res$freq             # actual frequencies
	vk  <- res$mtm$dpss$v       # slepian sequences
	dk  <- res$mtm$eigenCoefWt  # eigenweights (here the sqrt(|dk|^2))
	Yk  <- res$mtm$eigenCoefs   # eigencoefficients (fourier transform of tapered data)
	ts  <- data                 # zero-mean time series
		
	plot(res, log='yes', Ftest=FALSE, jackknife = FALSE)
	
	}
	
writeMat('wn_zp_NW5.mat', snn=snn, f = f, vk = vk, Yk = Yk, ts = ts, dk=dk, fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)