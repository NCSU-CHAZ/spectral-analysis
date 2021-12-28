## preamble
setwd("/Users/KatherineAnardeWheels/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM")

#data <- read.table("qpc_mat.txt")
data <- read.table("QPC3.txt")
N <- 64
NW <- 2
dt <- 1
K <- 4
C <- 64
b  <- NW/(N*dt)  # bandwidth
nFFT <- (N*2)/2+1;

snn    <- matrix(data=NA, nrow=nFFT, ncol=C)
f      <- matrix(data=NA, nrow=nFFT, ncol=C)
dofs   <- matrix(data=NA, nrow=nFFT, ncol=C)
vk     <- array(data=NA, dim = c(N, K, C))
dk     <- array(data=NA, dim = c(nFFT, K, C))
Yk     <- array(data=NA, dim = c(nFFT, K, C))
nFreq  <- matrix(data=NA, nrow=1, ncol=C)
nfft   <- matrix(data=NA, nrow=1, ncol=C)
ts     <- matrix(data=NA, nrow=N, ncol=C)

# compute the time vector for detrending
#time  <- seq((1), N*dt+1, dt) #adds one data point...too tired to fix
#time  <- time[1:N]
#ttbar <- time - (time[N] + time[1])/2

## non-adpative
for (i in 1:C){
	
	# find a polynomial trend using the Slepian tapers and calculate the residual
    #trend <- multitaperTrend(data[,i], B = b, deltat = dt, t.in = time)
    #data.resid <- data[,i] - trend[[1]] - trend[[2]] * ttbar
    	data.resid <- data[,i] - mean(data[,i])
    	
    # check that the detrended time series look right
    plot(time, data[,i], type="l", xlab="sample", ylab = "arb")
    lines(time, data[,i] - data.resid, type="l", col="red")
    
	# zero-pad
	res <- spec.mtm(data.resid, NW, K, deltat=1, plot=FALSE, nFFT="default", returnInternals=TRUE,
					adaptiveWeighting=FALSE, returnZeroFreq=TRUE)
	snn[,i]    <- res$spec             # don't multiply by two here (keep as a double-sided power spectra)
	f[,i]      <- res$freq             # actual frequencies
	vk[ , , i] <- res$mtm$dpss$v       # slepian sequences
	Yk[ , , i] <- res$mtm$eigenCoefs   # eigencoefficients (fourier transform of tapered data)
	nFreq[1,i] <- res$mtm$nfreqs       # number of frequencies
	nfft[1,i]  <- res$mtm$nFFT         # size of nfft
	ts[,i]     <- data.resid           # zero-mean time series

	}
	
writeMat('qpc_NAzp_nw3pt5k7.mat', snn=snn, f = f, vk = vk, Yk = Yk, nFreq = nFreq, nfft = nfft, ts = ts, fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)

## adaptive
for (i in 1:C){
	
	# find a polynomial trend using the Slepian tapers and calculate the residual
    #trend <- multitaperTrend(data[,i], B = b, deltat = dt, t.in = time)
    #data.resid <- data[,i] - trend[[1]] - trend[[2]] * ttbar
	data.resid <- data[,i] - mean(data[,i])
	    
    # check that the detrended time series look right
    plot(time, data[,i], type="l", xlab="sample", ylab = "arb")
    lines(time, data[,i] - data.resid, type="l", col="red")

	# zero-pad
	res <- spec.mtm(data.resid, NW, K, deltat=1, plot=FALSE, nFFT="default", returnInternals=TRUE,
					adaptiveWeighting=TRUE, returnZeroFreq=TRUE, Ftest=TRUE, jackknife=TRUE)
	snn[,i]    <- res$spec             # don't multiply by two here (keep as a double-sided power spectra)
	f[,i]      <- res$freq             # actual frequencies
	vk[ , , i] <- res$mtm$dpss$v       # slepian sequences
	dk[ , , i] <- res$mtm$eigenCoefWt  # eigenweights (here the sqrt(|dk|^2))
	Yk[ , , i] <- res$mtm$eigenCoefs   # eigencoefficients (fourier transform of tapered data)
	nFreq[1,i] <- res$mtm$nfreqs       # number of frequencies
	nfft[1,i]  <- res$mtm$nFFT         # size of nfft
	ts[,i]     <- data.resid           # zero-mean time series
		
	plot(res, log='yes', Ftest=FALSE, jackknife = TRUE)
	plot(res, Ftest=TRUE, siglines=c(0.90, 0.99))
	
	}
	
writeMat('qpc3_zp_nw2k4.mat', snn=snn, f = f, vk = vk, Yk = Yk, nFreq = nFreq, nfft = nfft, ts = ts, dk=dk, fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)

## full time-series (concatenate)
dataFull <- as.vector(as.matrix(data))
N <- length(dataFull)
NW <- 160
dt <- 1
K <- 320
C <- 1
b  <- NW/(N*dt)  # bandwidth
nFFT <- (N*2)/2+1;

snn    <- matrix(data=NA, nrow=nFFT, ncol=C)
f      <- matrix(data=NA, nrow=nFFT, ncol=C)
dofs   <- matrix(data=NA, nrow=nFFT, ncol=C)
vk     <- array(data=NA, dim = c(N, K, C))
dk     <- array(data=NA, dim = c(nFFT, K, C))
Yk     <- array(data=NA, dim = c(nFFT, K, C))
nFreq  <- matrix(data=NA, nrow=1, ncol=C)
nfft   <- matrix(data=NA, nrow=1, ncol=C)
ts     <- matrix(data=NA, nrow=N, ncol=C)

# compute the time vector for detrending
time  <- seq((1), N*dt+1, dt) #adds one data point...too tired to fix
time  <- time[1:N]
ttbar <- time - (time[N] + time[1])/2

# find a polynomial trend using the Slepian tapers and calculate the residual
trend <- multitaperTrend(dataFull, B = b, deltat = dt, t.in = time)
data.resid <- dataFull - trend[[1]] - trend[[2]] * ttbar
    
# check that the detrended time series look right
plot(time, dataFull, type="l", xlab="sample", ylab = "arb")
lines(time, dataFull - data.resid, type="l", col="red")
    
##data.resid <- data[,i] - mean(data[,i])
# zero-pad
res <- spec.mtm(data.resid, NW, K, deltat=1, plot=FALSE, nFFT="default", returnInternals=TRUE,
					adaptiveWeighting=TRUE, returnZeroFreq=TRUE, Ftest=TRUE, jackknife=TRUE)
snn   <- res$spec             # don't multiply by two here (keep as a double-sided power spectra)
f     <- res$freq             # actual frequencies
vk    <- res$mtm$dpss$v       # slepian sequences
dk    <- res$mtm$eigenCoefWt  # eigenweights (here the sqrt(|dk|^2))
Yk 	  <- res$mtm$eigenCoefs   # eigencoefficients (fourier transform of tapered data)
nFreq <- res$mtm$nfreqs       # number of frequencies
nfft  <- res$mtm$nFFT         # size of nfft
ts    <- data.resid           # zero-mean time series
		
plot(dropFreqs(res, 0, 0.5), Ftest=FALSE, jackknife = TRUE)
plot(res, Ftest=TRUE, siglines=c(0.90, 0.99))
	
writeMat('qpcFull_zp_nw160k320.mat', snn=snn, f = f, vk = vk, Yk = Yk, nFreq = nFreq, nfft = nfft, ts = ts, dk=dk, fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)
