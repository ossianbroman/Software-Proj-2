
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
set.seed(1435)

# Load the package
install.packages("boot")
library("boot")

# Create data and set seed to ensure results are reproducible
x <- fitness$Age
y <- fitness$Oxygen
regData <- data.frame(x, y)

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

# Parallelising

library("parallel")

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

# Carls improved code

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


# Using foreach
install.packages("foreach")
library("foreach")

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

# Make the function more efficient by using a helper function
parBadBoot <- function( index, scaleData, N )
{
  # sample data with replacement
  bootData <- scaleData[sample( 1:N, N, replace = T ),]
  Xmat <- bootData[,1:2]
  Ymat <- bootData[,3]
  
  # fit the model under this alternative reality
  # Changed the lm part to matrix form
  beta <- solve( t( Xmat ) %*% Xmat ) %*% t( Xmat ) %*% Ymat
  bootResults <- beta
}

# Bootstrap function that utilises parLapply function to speed up the process
bestBadBootBrother <- function( inputData, nBoot )
{
  # parallelize
  library(parallel)
  nCores <- detectCores()
  myClust <- makeCluster(nCores-1, type = "FORK")
  # scale inputData and pass to helper function in form of a matrix
  x <- cbind( 1, scale( inputData$x, scale = F ) )
  y <- scale( inputData$y, scale = F )
  scaleData <- as.matrix( cbind( x, y ) )
  
  bootResults <- array( dim = c( nBoot, 2 ) )
  
  N <- nrow( inputData )
  
  tempbootResults <- parLapply( myClust, 1:nBoot, parBadBoot, 
                                scaleData = scaleData, N = N ) 
  
  bootResults <- tempbootResults 
  # save all bootstap data into matrix and output
  r <- bootResults[[1]]
  c <- bootResults[[2]]
  m <- matrix( c( r, c ), ncol = 2, nrow = N )
  
  return( m ) 
}

system.time( test2 <- bestBadBootBrother( inputData = regData, 100000 ) )

# Extending the above helper function to accept an arbitraty number of covariates
parBootAnyCovars <- function( index, scaleData, N )
{
  # sample data with replacement
  bootData <- scaleData[sample( 1:N, N, replace = T ),]
  Xmat <- bootData[,1:(ncol(bootData)-1)]
  Ymat <- bootData[,ncol(bootData)]
  
  # fit the model under this alternative reality
  # Changed the lm part to matrix form
  beta <- solve( t( Xmat ) %*% Xmat ) %*% t( Xmat ) %*% Ymat
  return(beta)
}

# Include covariates by looping through
bestBootCovars1 <- function( nBoot, yDat, ... ) 
{
  # parallelize
  library(parallel)
  nCores <- detectCores()
  myClust <- makeCluster(nCores-1, type = "FORK")
  # save x covariates into list and loop through them to constuct covariate matrix
  xDat <- list(...)
  x <- cbind(1, scale(xDat[[1]], scale = F))
  if(length(xDat) > 1) {
    for(i in 2:length(xDat)) {
       x <- cbind(x, scale(xDat[[i]], scale = F))
    }
  }
  
  y <- scale( yDat, scale = F )
  
  # unite x covariates and y into matrix
  scaleData <- cbind( x, y )
  
  N <- length( yDat ) 
  
  # pass scaled data to helper function to bootstrap
  tempbootResults <- parLapply( myClust, 1:nBoot, parBootAnyCovars, 
                                scaleData = scaleData, N = N ) 
  # save all bootstap data into matrix and output
  m <- matrix( unlist(tempbootResults), ncol = length(xDat)+1, nrow = N )  
  return( m ) 
}

# Changing the for loop on xDat into something more dynamic
bestBootCovars2 <- function( nBoot, yDat, ... ) 
{
  # parallelize
  library(parallel)
  nCores <- detectCores()
  myClust <- makeCluster(nCores-1, type = "FORK")
  # save x covariates into list
  xDat <- list(...)
  # create matrix with x covariates and y and scale it
  mdata <- unlist(xDat)
  mcol <- length(xDat)
  mrow <- length( mdata ) / mcol
   
  baseone <- rep( 1 , mrow)
   
  mnew <- matrix( data = c(mdata, yDat), nrow = mrow, ncol = mcol+1)
  mscale <- scale( mnew, scale = F )
  scaleData <- matrix( data = c( baseone, mscale ), nrow = mrow, ncol = mcol + 2 )
  
  N <- length( yDat ) 
  # pass scaled data to helper function to bootstrap
  tempbootResults <- parLapply( myClust, 1:nBoot, parBootAnyCovars, 
                                scaleData = scaleData, N = N ) 
  # save all bootstap data into matrix and output 
  m <- matrix( unlist(tempbootResults), ncol = length(xDat)+1, nrow = N )
  
  return( m ) 
}

# Include covariates via matrix
bestBootCovars3 <- function( nBoot, yDat, ... ) 
{
  # parallelize
  library(parallel)
  nCores <- detectCores()
  myClust <- makeCluster(nCores-1, type = "FORK")
  
  # save x covariates into list, and turn them into a scaled matrix in one step
  xDat <- list(...)
  x <- cbind(1, scale(do.call(cbind, xDat), scale =F ))
  y <- scale( yDat, scale = F )

  scaleData <- cbind(x, y)
  
  N <- length( yDat ) 
  
  # pass scaled data to helper function to bootstrap
  tempbootResults <- parLapply( myClust, 1:nBoot, parBootAnyCovars, 
                                scaleData = scaleData, N = N ) 
  # save all bootstap data into matrix and output
  m <- matrix( unlist(tempbootResults), ncol = length(xDat)+1, nrow = N )
  
  return( m ) 
}

# Fastest version by scaling data in one step
oneBootToRuleThemAll <- function( nBoot, yDat, ... ) 
{
  # parallelize
  library(parallel)
  nCores <- detectCores()
  myClust <- makeCluster(nCores-1, type = "FORK")
  
  # save x covariates into list
  xDat <- list(...)
  # create data matrix in one step by using scale function only once
  scaleData <- cbind(1, scale(do.call(cbind, xDat, yDat), scale =F ))

  N <- length( yDat ) 
  
  # pass scaled data to helper function to bootstrap
  tempbootResults <- parLapply( myClust, 1:nBoot, parBadBootAnyCovars, 
                                scaleData = scaleData, N = N ) 
  # save all bootstap data into matrix and output
  m <- matrix( unlist(tempbootResults), ncol = length(xDat)+1, nrow = N )
  
  return( m ) 
}

## In all functions above, we return a matrix, where the first column is the intercept
## and the rest of the columns are best estimate for each covar (in slope terms)
## To have final bootstrap results, take column means and calculate 95% CI.

## Testing for 1 as well as more x covariates
system.time( test3 <- bestBootCovars1( 100000, regData$y, regData$x) )
system.time( test4 <- bestBootCovars2( 100000, regData$y, regData$x) )
system.time( test5 <- bestBootCovars3( 100000, regData$y, regData$x) )
system.time( test6 <- oneBootToRuleThemAll( 100000, regData$y, regData$x) )

system.time( test7 <- bestBootCovars1( 100000, regData$y, regData$x, sampleData$Weight, sampleData$RestPulse, sampleData$RunPulse, sampleData$MaxPulse) )
system.time( test8 <- bestBootCovars2( 100000, regData$y, regData$x, sampleData$Weight, sampleData$RestPulse, sampleData$RunPulse, sampleData$MaxPulse) )
system.time( test9 <- oneBootToRuleThemAll( 100000, regData$y, regData$x, sampleData$Weight, sampleData$RestPulse, sampleData$RunPulse, sampleData$MaxPulse) )


# Compare our efficient bootstrap to initial via microbenchmark
library(microbenchmark)
microbenchmark(lmBoot(regData, 100), oneBootToRuleThemAll( 100, regData$y, regData$x))















