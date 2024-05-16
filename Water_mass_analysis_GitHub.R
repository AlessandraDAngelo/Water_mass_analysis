# R script of simple 4-endmember mass balance analysis (Östlund and Hut, 1984) for the detection of the water masses.
# D'Angelo A. and Loose B. (2021) - University of Rhode Island, Graduate School of Oceanography

# Load necessary libraries
library(dplyr)

# Load dataframe "data"

# Standardize the end member data by subtracting the mean and dividing by the standard deviation
  #end members: Absolute salinity (SA), Arctic Nitrate-Phosphate (ANP), Oxygen isotope delta-18O (delta.O)
# 1. Calculate the mean
mean.SA <- mean(data$SA)
mean.ANP <- mean(data$ANP)
mean.delta.O <- mean(data$delta.O)

# 2. Calculate the standard deviation
sd.SA <- sd(data$SA)
sd.ANP <- sd(data$ANP)
sd.delta.O <- sd(data$delta.O)

# 3. Standardize data
data$SA.std <- (data$SA - mean.SA) / sd.SA
data$ANP.std <- (data$ANP - mean.ANP) / sd.ANP
data$delta.O.std <- (data$delta.O - mean.delta.O) / sd.delta.O

# Define endmembers
Sal <- data$SA.std
ANP <- data$ANP.std
delta.delta.O <- data$delta.O.std

# Newton et al 2013: End-Member Parameter Values for SIM: Surface + 2.6 for d18O, Surface for ANP
S.df <- subset(data, Depth >= -50)

# Pick surface average by CTD
S.df <- S.df %>%
  group_by(CTD) %>%
  summarise(s.ANP = mean(ANP), s.delta.O = mean(delta.O) + 2.6)

# Join averaged values with the main data frame
DF <- left_join(data, S.df)

# Create a loop
for (i in DF$CTD) {
  AW.test <- c(34.5, 0, 0, 1)
  PW.test <- c(32.5, 1, -2.5, 1)
  MW.test <- c(0, 0, -20, 1)
  SIM.test <- c(4, DF$s.ANP[i], DF$s.delta.O[i], 1)
  
  A1 <- cbind(AW.test, PW.test)
  A2 <- cbind(A1, MW.test)
  A <- cbind(A2, SIM.test)
  
  A.M <- array(A, dim = c(4, 4))
  
  # Standardize the matrix
  A.M[1, ] <- (A.M[1, ] - mean.SA) / sd.SA
  A.M[2, ] <- (A.M[2, ] - mean.ANP) / sd.ANP
  A.M[3, ] <- (A.M[3, ] - mean.delta.O) / sd.delta.O
  
  X <- array(A.M, dim = c(4, 4))
}

# Create weighted diagonal matrix
W <- diag(4)
W[4, 4] <- 1 / 0.05

# Loop through the data
for (i in seq_len(nrow(data))) {
  # Vector of observations
  y <- c(Sal[i], ANP[i], delta.delta.O[i], 1)
  
  # Ordinary least squares regression
  OLS.df <- solve(t(X) %*% W %*% X, t(X) %*% W %*% y)
  
  # Update the data frame with the results
  data$faw[i] <- OLS.df[[1]]
  data$fpw[i] <- OLS.df[[2]]
  data$fmw[i] <- OLS.df[[3]]
  data$fsim[i] <- OLS.df[[4]]
}

# Check the reproducibility of the data
# Misfit > 5% indicates not good data
f_sum <- rowSums(data[, c("faw", "fpw", "fmw", "fsim")])
data$OLS_Residuals <- (f_sum - 1) * 100

print(data$OLS_Residuals) #check results

#References:
#Newton, R., Schlosser, P., Mortlock, R., Swift, J., and MacDonald, R. 2013. Canadian Basin freshwater sources and changes: Results from the 2005 Arctic Ocean Section. J Geophys Res Oceans, 118, 2133–2154. https://doi.org/10.1002/jgrc.20101.
#Östlund, H. G., and G. Hut 1984, Arctic Ocean water mass balance from isotope data, J. Geophys. Res., 89(C4), 6373–6381, doi:10.1029/JC089iC04p06373.
