#### Repeatability & Power
# This code 95%+ of the code is written by Matt Wolak :: https://matthewwolak.github.io/software#ICC

# Boring data 
# 29 August 2018

# Simulate boring data ----------------------------------------------------
boring <- rnorm(100, 0, 1)   #<-- mean=0  &  sd=sqrt(variance)=1
hist(boring, xlab = "y", xlim = c(-5, 5))
# sample mean
mean(boring)
# sample var
var(boring)


# Monte Carlo -------------------------------------------------------------
# Without set.seed()
boring1 <- rnorm(100, 0, 1)
boring2 <- rnorm(100, 0, 1)

# means
mean(boring1); mean(boring2) 
# variances
var(boring1); var(boring2)


# Monte Carlo set seed ----------------------------------------------------
set.seed(101)  #<-- choose any old value, I was warned against numbers <100
boring1 <- rnorm(100, 0, 1)
set.seed(101)
boring2 <- rnorm(100, 0, 1)

# means
mean(boring1); mean(boring2) 
# variances
var(boring1); var(boring2)

# predicatably random, same result

# Repeatability calculation: ICC -- Intraclass correlation coefficient -------------------------------
# Set some values
mu <- 35     # Mean trait value
n <- 20      # Number of individuals in sample
k <- 2       # Number of measures per individual

ICC <- 0.5   # repeatability
sigma2_p <- 100  # Total phenotypic variance (sigma2_among + sigma2_within)
### END user-defined values #####

sigma2_among <- ICC * sigma2_p
sigma2_within <- sigma2_p - sigma2_among

# note the use of `gl()` below to 'generate levels'
## Creating the sample data has been simplified/unpacked for illustration
ICCdata <- data.frame(id = gl(n, k))

# Create the individual expectations p_i
## Note, we are using `rnorm()` below, but we have already set the seed above
pi_tmp <- rnorm(n, 0, sd = sqrt(sigma2_among))
ICCdata$pi <- rep(pi_tmp, each = k)

# Create the jth observation deviations from each p_i
ICCdata$rij <- rnorm(n*k, 0, sd = sqrt(sigma2_within))

# Use eqn. 2 to construct the trait phenotypes
ICCdata$y <- mu + ICCdata$pi + ICCdata$rij

boxplot(y ~ id, ICCdata,
        xlab = "ID", ylab = "y")


# ANOVA -------------------------------------------------------------------

anova1 <- aov(y ~ id, ICCdata)
summary(anova1)
# mean Sq is just Sum Sq/df
# F ratio = mean sq among / mean sq within


# From MS to Variance -----------------------------------------------------

# Lessels & Boag (1987)

MS <- anova(anova1)[["Mean Sq"]]  #<-- extracts the Mean Squares
Vwithin <- MS[2L]
# Note, below is only correct with equal sample size per group/individual (`k`)
Vamong <- (MS[1L] - MS[2L]) / k  #<-- NOTE different notation here than L&B 1987

ICC_1 <- Vamong / (Vamong + Vwithin)
ICC_1

# ICC package -------------------------------------------------------------
library(ICC)
ICCest(x = id, y = y, data = ICCdata)

# ICC = 0.5616 which is different than what we specified "0.5" in the data


# Monte Carlo error or just error? ----------------------------------------

# expected variance among (V):
(Vhat_among <- var(pi_tmp))

# expected variance within:
(Vhat_within <- var(ICCdata$rij))

# expected ICC
round(Vhat_among / (Vhat_among + Vhat_within), 4)


# Other Monte Carlo sources of error --------------------------------------

# There are other sources of Monte Carlo sampling variation that could contribute to the differenc

# POWER = probability of rejecting the null hypothesis when it is false and an alternative hypothesis is correct 

# Function to run what we ran individual above called "ICCdataGen":

## ICC repeatability
## n Number of individuals in sample
## k Number of measures per individual
## mu Mean trait value
## sigma2_p Total phenotypic variance (sigma2_among + sigma2_within)

# I've slightly changed the order and decided not to give default values
ICCdataGen <- function(ICC, n, k, mu, sigma2_p){
  sigma2_among <- ICC * sigma2_p
  sigma2_within <- sigma2_p - sigma2_among
  
  # note the use of `gl()` below to 'generate levels'
  ## Creating the sample data has been condensed for shorter/optimal code
  ICCdata <- data.frame(id = gl(n, k),
                        pi = rep(rnorm(n, 0, sd = sqrt(sigma2_among)), each = k),  # individual expectations p_i
                        rij = rnorm(n*k, 0, sd = sqrt(sigma2_within)))  # jth observation deviations from each p_i
  
  # Use eqn. 2 to construct the trait phenotypes
  ICCdata$y <- mu + rowSums(ICCdata[, c("pi", "rij")])
  
  vhat_among = var(unique(ICCdata$pi))
  vhat_within = var(ICCdata$rij)
  
  return(list(ICCdata, vhat_among, vhat_within))  #<-- returns a list with object data.frame and two numeric objects
  # returning list like this makes going back to re-run the analysis with many fake data.frame
  
}  #<-- END function


tmp_df <- ICCdataGen(ICC = 0.5, n = 20, k = 2, mu = 35, sigma2_p = 100)
ICCest(x = id, y = y, data = tmp_df[[1]])$ICC  # call for tmp_df[[1]] because tmp_df has the dataframe


icc50 <- replicate(50, ICCest(x = id, y = y, data = ICCdataGen(ICC = 0.5, n = 20, k = 2, mu = 35, sigma2_p = 100))$ICC)
icc50
mean(icc50)
# Plot result
boxplot(icc50)

# POWER SIMULATION --------------------------------------------------------

## n Number of individuals in sample
## k Number of measures per individual
## alpha criteria for rejecting a Null Hypothesis
## ICC repeatability
## mu Mean trait value
## sigma2_p Total phenotypic variance (sigma2_among + sigma2_within)

# Again, I've slightly changed the order and decided not to give default values
anovaPfun <- function(n, k, alpha, ICC, mu, sigma2_p){
  sigma2_among <- ICC * sigma2_p
  sigma2_within <- sigma2_p - sigma2_among
  
  # note the use of `gl()` below to 'generate levels'
  ## Creating the sample data has been condensed for shorter/optimal code
  ICCdata <- data.frame(id = gl(n, k),
                        pi = rep(rnorm(n, 0, sd = sqrt(sigma2_among)), each = k),  # individual expectations p_i
                        rij = rnorm(n*k, 0, sd = sqrt(sigma2_within)))  # jth observation deviations from each p_i
  
  # Use eqn. 2 to construct the trait phenotypes
  ICCdata$y <- mu + rowSums(ICCdata[, c("pi", "rij")])
  
  # Now, run the one-way ANOVA and generate a summary table
  tmp_aov <- summary(aov(y ~ id, data = ICCdata))
  # extracts the p-value from the ANOVA table
  tmp_pval <- tmp_aov[[1L]][['Pr(>F)']][1L] 
  
  # assesses whether p-value < alpha and return result
  return(tmp_pval < alpha)  #<-- returns an object of class logical
  
}  #<-- END function

# Check output and replicate
# First, just run once to see the output
anovaPfun(n = 20, k = 2, alpha = 0.05, ICC = 0.5, mu = 35, sigma2_p = 100)

# Now use `replicate()`
## wrap it all in `system.time()` to report how long it took
system.time(pTF <- replicate(1000,
                             anovaPfun(n = 20, k = 2, alpha = 0.05, ICC = 0.5, mu = 35, sigma2_p = 100)))

## have a look-see at the first 100 values
head(pTF, 100)

# What is the POWER:
## Conveniently, R can do math on `logicals` (treats TRUE=1 and FALSE=0)
## The proportion of TRUE out of the total (e.g., `length(pTF)`) is simply
## the average of the 0s and 1s
mean(pTF)  #<-- POWER!
# general rule of thumb is to get 0.8 power
# This helps us determine how many n (number of indv in sample) and k (number of measures per indvl) is adequate



# POWER CURVE assignment --------------------------------------------------
# This is where students start writing...

# 1 Report your power for each sample size in a table.
# 2 Make a figure depicting the curve over the range of sample sizes.

# Condition
# ICC is 0.25
# n is at least 10

# sample size
n <- seq(10, 100, 10) #creates a vector 
vec_output <- rep(NA, length(n)) # creates an empty vector the same length as n

# Report power for each sample size | you can use 1:length(n) as well, here I use seq_along(n)

# This for loop does the trick
for (i in seq_along(n)) {
  
  vec_output[i] <- mean(replicate(1000, #replicate 1000 times the function anovaPfun
                                  anovaPfun(n = n[i], k = 2, alpha = 0.05, ICC = 0.5, 
                                            mu = 35, sigma2_p = 100)))
  
  # prints the output individually
  # print(vec_output[i]) 
}

# Report power for each sample size
# Power Curve
plot(vec_output ~ n)


# More fancy output -------------------------------------------------------
# store the output from the for loop in a data.frame 
df_output <- data.frame(vec_output)

# add column with values of n to df_output data frame. 
df_output$n <- data.frame(n)

# the data frame df_output now looks like this
df_output

# Fancier plot with ggplot
library(ggplot2)
ggplot(df_output, aes(x = n, y = vec_output)) +
  geom_point() +
  labs(title = "Power Curve", x = "Individual (n)", y = "Power") +
  geom_line(y = 0.8, aes(colour = "red"))


