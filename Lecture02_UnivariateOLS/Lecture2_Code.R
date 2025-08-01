########## Lecture2_Code.R
# Creator: Alex Hoagland, alcobe@bu.edu
# Created: 5/31/2022
# Last modified: 5/31/2022
#
# PURPOSE
#   1. Probability review
#   2. OLS Regression
#
# NOTES: 
#   - uses the Tidyverse package and Dplyr
################################################################################


##### Packages #####
# install.packages('tidyverse') # if needed, install the package
library(tidyverse) # call the relevant library
library(broom) # Useful package for cleaning up regression output
###################################


##### Probability Review #####
p_healthy <- 0.75 # Assign 3 main probabilities
p_acute <- 0.2
p_chronic <- 0.05

p_acute_chronic <- 0.02 # Assign joint probabilities
p_healthy_acute <- 0
p_healthy_chronic <- 0
p_healthy_acute_chronic <- 0

# Conditional probabilities
p_chronic_givenacute <- p_acute_chronic / p_acute

# Test of independence
p_chronic_givenacute == p_chronic * p_acute

# Total probability 
day1 <- p_chronic # On first day, probability of illness is the known probability
day2 <- p_chronic * (1-p_chronic) # The probability of being healthy on day1, then developing an illness on the second day 
totalprob <- day1+day2

# Question: how long until total prob > 0.15?
dayprob <- function(n) { # This function returns the probability of (i) being healthy until day n, then (ii) getting a chronic condition on day n 
  return(p_chronic * (1-p_chronic)^n)
} 
mydata <- tibble( x=seq(1,10,by=1) ) %>% mutate(dayprob=dayprob(x)) %>% mutate(totalprob = cumsum(dayprob))
  # This builds a data set examining 1 through 10 days, calculates the daily probability and total probability of getting a total illness in each day
  # Note: uses a new dplyr function "cumsum" to generate a cumulative sum
print(paste0("The number of days before totalprob > 0.15 is ", 
             mydata %>% mutate(tokeep = totalprob>0.15) %>% filter(tokeep==1) %>% mutate(first = min(x)) %>% 
               filter(row_number()==1) %>% pull(first),sep=""))
  # In order to show the answer in the report

# Bayes' Rule
# First, we define an update function for p_chronic
update_chronic <- function(oldprob) { 
  p_x1_you_given_x3_me <- .7 # Experiment with changing this probability to see how the graph below changes
    # What happens if this is .75? .1? What does changing this mean in terms of the example? 
  p_x1_you_given_x1_me <- .75
  num <- p_x1_you_given_x3_me * oldprob
  den <- p_x1_you_given_x3_me * oldprob + p_x1_you_given_x1_me * (1-oldprob)
  return(num/den)
}

# Now, loop through 100 points at which family members are healthy
mydata <- rep(NA, 100) # Set blanks for 100 days
mydata[1] <- .05 # Set first probability at .05
for (i in 2:100) { 
  mydata[i] <- update_chronic(mydata[i-1])
}
# Now plot
mydata <- tibble(mydata) %>% mutate(x=row_number())
ggplot(mydata,aes(x,mydata)) + geom_point() + theme_minimal() + labs(x = "Day", y = "P(x_3)")
########################################


##### Regression #####
set.seed(1) # Setting the seed helps to make random number generators give us the same numbers across machines
tb <- tibble(
  x = rnorm(10000),
  u = rnorm(10000),
  y = 5.5*x + 12*u
) # Here's our true model, with some randomness baked in

reg_tb <- tb %>% 
  lm(y ~ x, .) %>%
  print() # What happens if we just regress y on x?

# Looking at your model output
summary(reg_tb) # gives you a full summary
regdata <- tidy(reg_tb, conf.int=TRUE) # Stores regression output in a data frame
regouts <- glance(reg_tb) # Stores other regression features in a data frame

tb <- tb %>% 
  mutate(
    yhat1 = predict(lm(y ~ x, .)),
    yhat2 = 0.0732608 + 5.685033*x, 
    uhat1 = residuals(lm(y ~ x, .)),
    uhat2 = y - yhat2
  ) # How close are our predictions? 

summary(tb[-1:-3])

tb %>% # Let's talk about ggplot for a second
  lm(y ~ x, .) %>% 
  ggplot(aes(x=x, y=y)) + 
  ggtitle("OLS Regression Line") +
  geom_point(size = 0.05, color = "black", alpha = 0.5) +
  geom_smooth(method = lm, color = "black") +
  annotate("text", x = -1.5, y = 30, color = "red", 
           label = paste("Intercept = ", -0.0732608)) +
  annotate("text", x = 1.5, y = -30, color = "blue", 
           label = paste("Slope =", 5.685033)) + 
  theme_classic() + labs(x="X", y="Y")

### Calculate the R-squared
SST <- tb %>% mutate(ybar=mean(y), sst=(y-ybar)^2) %>% select(sst) %>% summarise(sum(sst))
SSE <- tb %>% mutate(ybar=mean(y), sse=(yhat1-ybar)^2) %>% select(sse) %>% summarise(sum(sse))
SSR <- tb %>% mutate(ssr=uhat1^2) %>% select(ssr) %>% summarise(sum(ssr))

R2 <- SSE/SST
  
# How do you pull R2 from a lm command? 
myreg <- lm(tb$y ~ tb$x)
summary(myreg)

### Show unbiasedness: 
# Start with a population model y = 3 + 2x + e, 
# Assume x is normally distributed N(0,9) and e is N(0,36)
# Assume mean independence 

lm <- lapply( # lapply is a useful command (see ?lapply)
  1:1000,
  function(x) tibble(
    x = 9*rnorm(10000), # Distribution for x (rnorm(N) creates N random draws from a N(0,1) distribution)
    u = 36*rnorm(10000), # Distribution for e
    y = 3 + 2*x + u
  ) %>% 
    lm(y ~ x, .)
) # This runs 1,000 regressions with new data every time

as_tibble(t(sapply(lm, coef))) %>%
  summary(x) # Summarize the regression coefficients (tabular)

as_tibble(t(sapply(lm, coef))) %>% 
  ggplot()+
  geom_histogram(aes(x), binwidth = 0.01) # Summarize the regression coefficients (graph)

# Clean up data viz
as_tibble(t(sapply(lm, coef))) %>% 
  ggplot()+
  geom_histogram(aes(x), binwidth = 0.01,fill="gray",color='black') + 
  theme_minimal() + 
  labs(x="Estimated Slope Coefficient", y="Count",title="Distribution of Slope Coefficients") + 
  geom_vline(xintercept=2,color='red')
#########################################