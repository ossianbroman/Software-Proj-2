
# 1 Modify the the R bootstrapping code to be more efficient. 
#   It should also be parallelised at some level. 
#   You should profile both the original and final versions as well as determining the overall speed increase. 
#   Include the profile file in your repository.
# 2 Your new bootstrap function in R should be altered to accept an arbitrary number of covariates.
# 3 Modifications to the bootstrap function should be version controlled via your project GitHub repository.
# 4 Micro-benchmark (package microbenchmark) your R bootstrap against bootstraps via the package boot.
# 5 NB Your improved function will be used for your individual submissions

# Import the data
fitness <- read.csv("~/Desktop/5763- A2/Software-Proj-2/data/fitness.csv")

# Load the package
library(boot)

# Create data and set seed to ensure results are reproducible
x <- fitness$Age
y <- fitness$Oxygen
regData <- data.frame(x, y)

################

# Carls function to pass boot
lmBoot <- function(inputData, nBoot) {
  
  for(i in 1:nBoot){
    
    # resample our data with replacement
    bootData <- inputData[sample(1:nrow(inputData), nrow(inputData), replace = T),]
    
    # fit the model under this alternative reality
    bootLM <- lm(y ~ x, data = bootData)
    
    # store the coefs
    if(i == 1){
      bootResults <- matrix(coef(bootLM), ncol = 2)
    } else {
      bootResults<- rbind(bootResults, matrix(coef(bootLM), ncol = 2))
    }
    
  } 
  
  # end of i loop
  
  bootResults
}

system.time(test <- lmBoot(regData, nBoot = 100000))

##############

# Parallelising

library(parallel)

detectCores()

nCores <- detectCores()
nCores

myClust <- makeCluster(nCores-1, type = "FORK")

clusterEvalQ(myClust, library(boot))

sampleData <- fitness

clusterExport(myClust, "fitness")
clusterEvalQ(myClust, dataSum <- sum(fitness))

dataSum <- clusterEvalQ(myClust, sum(fitness))
dataSum

################ Carls improved code

lmBoot1 <- function(inputData, nBoot) {
  
  x <- cbind(1, scale(inputData$x, scale = F))
  y <- scale(inputData$y, scale = F)
  scaleData <- as.matrix(cbind(x, y))
  
  bootResults <- array(dim=c(nBoot, 2))
  
  for(i in 1:nBoot){
    
    # resample our data with replacement
    bootData <- scaleData[sample(1:nrow(inputData), nrow(inputData), replace = T),]
    xmat <- bootData[,1:2]
    ymat <- bootData[,3]
    
    # fit the model under this alternative reality
    # Changed the lm part to matrix form
    
    beta <- solve(t(xmat) %*% xmat) %*% t(xmat) %*% ymat
    
    bootResults[i,] <- beta
    
  } # end of i loop
  
  bootResults
  
}

system.time(lmBoot1(regData, nBoot = 100000))


################## Using foreach

library('foreach')

lmBoot2 <- function(inputData, nBoot) {
  
  x <- cbind(1, scale(inputData$x, scale = F))
  y <- scale(inputData$y, scale = F)
  scaleData <- as.matrix(cbind(x, y))
  
  bootResults <- array(dim=c(nBoot, 2))
  
  foreach(i = 1:nBoot)  %dopar% {
    
    # resample our data with replacement
    bootData <- scaleData[sample(1:nrow(inputData), nrow(inputData), replace = T),]
    xmat <- bootData[,1:2]
    ymat <- bootData[,3]
    
    # fit the model under this alternative reality
    # Changed the lm part to matrix form
    
    beta <- solve(t(xmat) %*% xmat) %*% t(xmat) %*% ymat
    
    bootResults[i,] <- beta
    
  } # end of i loop
  
  bootResults
  
}

system.time(lmBoot2(regData, nBoot = 100000))

################# Efficient

parBadBoot <- function( index, scaleData, N )
{
  bootData <- scaleData[sample( 1:N, N, replace = T ),]
  Xmat <- bootData[,1:2]
  Ymat <- bootData[,3]
  
  # fit the model under this alternative reality
  # Changed the lm part to matrix form
  beta <- solve( t( Xmat ) %*% Xmat ) %*% t( Xmat ) %*% Ymat
  bootResults <- beta
}

# final bootstrap function that utilises parLapply function to speed
# up process
bestBadBootBrother <- function( inputData, nBoot )
{
  x <- cbind( 1, scale( inputData$x, scale = F ) )
  y <- scale( inputData$y, scale = F )
  scaleData <- as.matrix( cbind( x, y ) )
  
  bootResults <- array( dim = c( nBoot, 2 ) )
  
  N <- nrow( inputData )
  
  tempbootResults <- parLapply( myClust, 1:nBoot, parBadBoot, 
                                scaleData = scaleData, N = N ) 
  
  bootResults <- tempbootResults 
  
  r <- bootResults[[1]]
  c <- bootResults[[2]]
  m <- matrix( c( r, c ), ncol = 2, nrow = N )
  
  return( m ) 
}

system.time( test2 <- bestBadBootBrother( inputData = regData, 100000 ) )

##################
# Exteding above to any number of covars

parBadBootAnyCovars <- function( index, scaleData, N )
{
  bootData <- scaleData[sample( 1:N, N, replace = T ),]
  Xmat <- bootData[,1:(ncol(bootData)-1)]
  Ymat <- bootData[,ncol(bootData)]
  
  # fit the model under this alternative reality
  # Changed the lm part to matrix form
  beta <- solve( t( Xmat ) %*% Xmat ) %*% t( Xmat ) %*% Ymat
  bootResults <- beta
}

bestBadBootBrotherAnyCovars <- function( nBoot, yDat, ... ) 
{
  xDat <- list(...)
  x <- cbind(1, scale(xDat[[1]], scale = F))
  if(length(xDat) > 1) {
    for(i in 2:length(xDat)) {
      x <- cbind(x, scale(xDat[[i]], scale = F))
    }
  }
  y <- scale( yDat, scale = F )
  scaleData <- as.matrix( cbind( x, y ) )
  bootResults <- array( dim = c( nBoot, 2 ) )
  
  N <- length( yDat ) 
  
  tempbootResults <- parLapply( myClust, 1:nBoot, parBadBoot, 
                                scaleData = scaleData, N = N ) 
  
  bootResults <- tempbootResults 
  
  r <- bootResults[[1]]
  c <- bootResults[[2]]
  m <- matrix( c( r, c ), ncol = 2, nrow = N )
  
  return( m ) 
}

system.time( test3 <- bestBadBootBrotherAnyCovars( 100, regData$y, regData$x ) )

# Compare our bootstrap to initial via microbenchmark
library(microbenchmark)
microbenchmark(lmBoot(regData, 100), bestBadBootBrotherAnyCovars( 100, regData$y, regData$x))















