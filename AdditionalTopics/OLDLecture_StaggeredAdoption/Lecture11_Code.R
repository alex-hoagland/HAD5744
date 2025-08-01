########## Lecture11_Code.R
# Creator: Alex Hoagland, alcobe@bu.edu
# Created: 8/15/2022
# Last modified: 
#
# PURPOSE
#   Advances in Difference-in-differences literature
#
# NOTES: 
#   - uses the Tidyverse package and Dplyr
################################################################################


##### Packages #####
# install.packages('tidyverse') # if needed, install the package
library(tidyverse) # call the relevant library
library(faux) # Useful package for simulating data
library(modelsummary) 
library(causaldata) 
library(here)
library(foreign)
library(zoo)

# TWFE Packages
library(fixest)
library(TwoWayFEWeights)
library(bacondecomp)
library(DIDmultiplegt)
library(did)
library(did2s)

set.seed(03262020)
##########


##### 1. Testing Parallel Trends #####
mydata <- read.csv(here("us_state_vaccinations.csv"))
mydata$date <- as.yearmon(mydata$date)

# Trim out some "states" that we don't need
`%!in%` <- Negate(`%in%`)
mydata <- mydata %>% filter(location  %!in% c("American Samoa", "Bureau of Prisons", "Dept of Defense",
                                              "Federated States of Micronesia", "Guam", "Indian Health Svc",
                                              "Long Term Care", "Marshall Islands", "Northern Mariana Islands",
                                              "Puerto Rico", "Republic of Palau", "United States", "Veterans Health",
                                              "Virgin Islands"))

# Group to month level
mydata <- mydata %>% group_by(date, location) %>% 
  summarize(total_vaccinations_per_hundred = sum(total_vaccinations_per_hundred,na.rm=T),
            people_vaccinated_per_hundred = sum(people_vaccinated_per_hundred,na.rm=T),
            monthly_vaccinations_per_million = sum(daily_vaccinations_per_million,na.rm=T))

# Identify all treated states for our TWFE regression
regdata_dte <- mydata %>% mutate(treatdate = ifelse(location %in% c("Ohio", "New York State", "Maryland"), 2021.417,
                                                    ifelse(location %in% c("Massachusetts", "Michigan"), 2021.583,NA)), # Note the treatment date doesn't matter for non-treated states
                                 reltime = round((as.numeric(date) - treatdate)*12)) 
regdata_dte <- regdata_dte %>% mutate(state = ifelse(location %in% c("Ohio", "New York State", "Maryland", "Massachusetts", "Michigan"), 1, 0),
                                      logy = log(monthly_vaccinations_per_million+1))
# Make sure to reassign "treated" states

dte <- feols(logy ~ i(reltime, ref = -1) | 
               date + location, data = regdata_dte)

# And use coefplot() for a graph of effects
iplot(dte,ref.line=-0.2, 
      ref.line.par=list(col="red",lty=1,lwd=2)) # How do we interpret this? 
################################################################################


###### 3. Diagnostics: Weighting
# This uses TwoWayFEWeights
regdata_dte <- regdata_dte %>% ungroup()%>% mutate(treated = ifelse(reltime >= 0,1,0))
twowayfeweights(df=regdata_dte, # Dataframe
                Y="logy", # Dep var
                G="location", # group identifier
                T="date", # Date var
                D="treated", # Whether state was treated
                cmd_type="feTR") # type of estimation to perform -- most commonly feTR
weighttable <- twowayfeweights(df=regdata_dte, # Dataframe
                               Y="logy", # Dep var
                               G="location", # group identifier
                               T="date", # Date var
                               D="treated", # Whether state was treated
                               cmd_type="feTR") # type of estimation to perform -- most commonly feTR
# lots of negative weighting here!

# Another option: bacondecomp
regdata_dte <- regdata_dte %>% mutate(treated = ifelse(is.na(treated),0,treated))
regdata_dte <- regdata_dte %>% mutate(stateid = as.numeric(as.factor(location)),
                                      dateid=as.numeric(as.factor(date))) # Need a state id variable as a number, not string
bacon(formula = logy ~ treated, 
      data=regdata_dte,
      id_var="stateid",
      time_var="dateid",
      quietly=F)
################################################################################


##### 4. $DID_M$ Estimator #####
did_multiplegt(df=regdata_dte, Y="logy", G="stateid", T="date", D="treated")
# This gives us a point estimate of the treatment effect. What about inference? 

# Add boostraps! Obviously this takes time 
did_multiplegt(df=regdata_dte, Y="logy", G="stateid", T="date", D="treated", brep=100)

# Can also incorporate placebos to check this
didm_placebo <- did_multiplegt(df=regdata_dte, Y="logy", G="stateid", T="date", D="treated", 
                               placebo=17) # Max number of placebos: T-2 --> more placebos = more time
placebos <- rep(NA, 17)
for (i in 1:17) {
  placebos[i] <- didm_placebo[2*i+2] # Pull placebo estimates
}
df <- data.frame(placebos[!is.na(placebos)]) # Remove all NANs, turn into data frame
df <- data.frame(t(df), "Placebos")
names(df) <- c("Placebo", "test")
ggplot(df,aes(x=Placebo))+geom_histogram(bins=5,color='darkslategray4',fill='darkslategray3')+
  geom_vline(xintercept=as.numeric(didm_placebo[1]),size=2,color='red')+
  theme_classic()+labs(x="Estimated Effect SIze")
################################################################################


##### 5. Callaway and Sant'Anna #####
# Need a variable for first period of treatment (see below)
# Note: there has to be a cleaner way to make this work, but I haven't found it yet. 
regdata_dte <- regdata_dte %>%
  mutate(firsttreat = ifelse(treated==1,dateid,NA)) %>%
  group_by(stateid) %>% 
  mutate(firsttreat = min(firsttreat,na.rm=T)) %>% # first period where observation is treated
  mutate(firsttreat = ifelse(is.infinite(firsttreat),0,firsttreat)) # want untreated units to have a 0

csdid <- att_gt(yname='logy',
                tname='dateid',
                idname='stateid',
                gname='firsttreat', 
                # Note: needs to be a variable containing first period that each unit was treated 
                # (0 for untreated units)
                xformla = ~ total_vaccinations_per_hundred, # if there are covariates to include
                data=regdata_dte,
                control_group="notyettreated", # Default here is "nevertreated" only
                anticipation=2) # Can set # of periods anticipation might occur
summary(csdid) # how many of these are significant? 
df <- data.frame(csdid$group,csdid$t,csdid$att)

# Visualizations
ggplot(df,aes(csdid.att)) + 
  geom_histogram(fill='darkslategray3',color='darkslategray4',bins=15) +
  theme_classic()+labs(x="Estimated ATT",y="",title="ATT: All (g,t) cells")+
  geom_vline(xintercept = 0,size=1.5,color='red')
# Look at all that heterogeneity!

# Heterogeneity across groups
ggplot(df,aes(csdid.att,group=factor(csdid.group),fill=factor(csdid.group))) + 
  geom_histogram(bins=15) +
  theme_classic()+labs(x="Estimated ATT",y="",title="ATT By Groups",fill="Group")+
  geom_vline(xintercept = 0,size=1.5,color='red')
# Treated groups have different ATT distributions

# Heterogeneity across periods 
# Let's bin time periods
df <- df %>% mutate(bintime = ifelse(csdid.t <= 6, 1, ifelse(csdid.t <= 12, 2, 3)))
ggplot(df,aes(csdid.att,group=factor(bintime),fill=factor(bintime))) + 
  geom_histogram(bins=15) +
  theme_classic()+labs(x="Estimated ATT",y="",title="ATT By Periods",fill="Periods (binned)")+
  geom_vline(xintercept = 0,size=1.5,color='red')
# Treated groups have different ATT distributions
################################################################################


##### 6. Borusyak et al., Gardner #####
regdata_dte <- regdata_dte %>% mutate(reltime = ifelse(is.na(reltime),Inf,reltime))
# Works better if you replace NAs with Infs for untreated group

regdata_dte <- regdata_dte %>% mutate(reltime = factor(reltime,
                                                       levels=c(-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,Inf)))
# Need reltime to be a factor before being passed through feols().
# Just note that the levels are out of order here

es <- did2s(data=regdata_dte,
            yname = "logy", first_stage = ~ total_vaccinations_per_hundred | stateid + dateid, 
            second_stage = ~i(reltime,ref=c(-1,Inf)), treatment = "treated", 
            cluster_var = "dateid")

fixest::iplot(es, main = "Event study: Staggered treatment", xlab = "Relative time to treatment", col = "steelblue", ref.line = -0.5)

# Compare to the original event study plot (dte) 
dte <- feols(logy ~ i(reltime, ref = c(-1, Inf)) | 
               date + location, data = regdata_dte)

fixest::iplot(list(es, dte), sep = 0.2, ref.line = -0.5,
              col = c("steelblue", "#82b446"), pt.pch = c(20, 18), 
              xlab = "Relative time to treatment", 
              main = "Event study: Staggered treatment (comparison)")


# Legend
legend(x=12, y=-.45, col = c("steelblue", "#82b446"), pch = c(20, 18), 
       legend = c("Two-stage estimate", "TWFE"))
################################################################################


##### 7. All the Estimators! #####
# Note: this comes from the did2s package again 
out = event_study(
  data = regdata_dte,
  yname = "logy", 
  idname = "stateid",
  tname = "dateid",
  gname = "firsttreat", # Unit specific date of initial treatment
  estimator = "all" # Can ask for some of c("all", "TWFE", "did2s", "did", "impute", "sunab", "staggered")
)

plot_event_study(out)
################################################################################


##### 8. Wooldridge's ETWFE/POLS #####
# To run this, we just need interactions between all cohorts and periods. 
# We already have periods, let's add in cohorts
regdata_dte <- regdata_dte %>% mutate(cohortid = as.factor(firsttreat))
pols <- lm(logy ~ factor(cohortid) + factor(dateid) + factor(cohortid):factor(dateid), 
           data=regdata_dte)
summary(pols)
################################################################################
