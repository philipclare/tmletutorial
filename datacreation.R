################################################
### Data creation code                        ##
################################################

cloudstor <- "C:/Users/z3312911/Cloudstor/"
.libPaths(paste0(cloudstor,"R Library"))

library("tmle")
library("ltmle")
library("SuperLearner")
library("simcausal")
library("MASS")
library("ranger")
library("parallel")
library("doParallel")
library("foreach")
library("lme4")

RNGkind(kind="default",normal.kind="default")
set.seed(43236)

# Data creation using 'simcausal'
# Creates 3 z-distributed baseline (time-constant) variables (ba, bb, bc)
# Then create a normally distruted 'latent' variable, u_t
# Initialise the confounder 'l' based on baseline covariates and 'u'
# Initialise the exposure variable based on baseline covariates and initial 'l'
# For each t in 1-4, create a new version of u, l and a based on previous observations
# Create y variables for each t, based on ALL u,l,a 'prior' to that y, as well as y_t-1
D <- DAG.empty() + 
  node("ba", distr="rnorm", mean=0, sd = 1) +
  node("bb", distr="rnorm", mean=0, sd = 1) +
  node("bc", distr="rnorm", mean=0, sd = 1) +
  
  node("u", t=0, distr="rnorm", mean=0, sd = 1) +
  node("c", t=0, distr="rbern", prob=0) +
  node("l", t=0, distr="rbern", prob=plogis(-2 + 1.5*u[t] + 0.1*ba - 0.1*bb + 0.1*bc)) + 
  node("a", t=0, distr="rbern", prob=plogis(-2 + 1.5*l[t] + 0.2*ba - 0.2*bb + 0.2*bc)) +
  
  node("u", t=1:4, distr="rnorm", mean=0.7*u[t-1], sd = 1) +
  node("c", t=1:4, distr="rbern", prob=ifelse(c[t-1]==1,1,plogis(-4.75 + 2.0*a[t-1] + 2.0*l[t-1]))) +
  node("l", t=1:4, distr="rbern", prob=ifelse(c[t]==1,NA,plogis(-2 + 1.0*a[t-1] + 2.0*l[t-1] + 1.5*u[t] + 0.1*ba - 0.1*bb + 0.1*bc))) +
  node("a", t=1:4, distr="rbern", prob=ifelse(c[t]==1,NA,plogis(-2 + 2.0*a[t-1] + 1.5*l[t] + 0.2*ba - 0.2*bb + 0.2*bc))) +
  
  node("y", t=0, distr="rnorm", mean=(1.00*a[t]
                                      + 0.50*l[t]
                                      + 0.50*u[t]
                                      + 0.2*ba - 0.2*bb + 0.2*bc), sd=1) +
  node("y", t=1, distr="rnorm", mean=ifelse(c[t]==1,NA,(0.80*a[t-1] + 1.00*a[t]
                                                        + 0.50*l[t-1] + 0.50*l[t]
                                                        + 0.50*u[t-1] + 0.50*u[t]
                                                        + 0.2*ba - 0.2*bb + 0.2*bc)), sd=1) +
  node("y", t=2, distr="rnorm", mean=ifelse(c[t]==1,NA,(0.60*a[t-2] + 0.80*a[t-1] + 1.00*a[t]
                                                        + 0.50*l[t-2] + 0.50*l[t-1] + 0.50*l[t]
                                                        + 0.50*u[t-2] + 0.50*u[t-1] + 0.50*u[t]
                                                        + 0.2*ba - 0.2*bb + 0.2*bc)), sd=1) +
  node("y", t=3, distr="rnorm", mean=ifelse(c[t]==1,NA,(0.40*a[t-3] + 0.60*a[t-2] + 0.80*a[t-1] + 1.00*a[t]
                                                        + 0.50*l[t-3] + 0.50*l[t-2] + 0.50*l[t-1] + 0.50*l[t]
                                                        + 0.50*u[t-3] + 0.50*u[t-2] + 0.50*u[t-1] + 0.50*u[t]
                                                        + 0.2*ba - 0.2*bb + 0.2*bc)), sd=1) +
  node("y", t=4, distr="rnorm", mean=ifelse(c[t]==1,NA,(0.20*a[t-4] + 0.40*a[t-3] + 0.60*a[t-2] + 0.80*a[t-1] + 1.00*a[t]
                                                        + 0.50*l[t-4] + 0.50*l[t-3] + 0.50*l[t-2] + 0.50*l[t-1] + 0.50*l[t]
                                                        + 0.50*u[t-4] + 0.50*u[t-3] + 0.50*u[t-2] + 0.50*u[t-1] + 0.50*u[t]
                                                        + 0.2*ba - 0.2*bb + 0.2*bc)), sd=1)
# Set this causal structure, defining all 'u' variables as latent (so they will not be included in the data)
D <- suppressWarnings(set.DAG(D, latent.v = c("u_0","u_1","u_2","u_3","u_4")))
# Create final simulated dataset
ldata <- simcausal::sim(D,n=1000)

cov(ldata$a_0,ldata$l_1,use="complete.obs")
cov(ldata$a_1,ldata$l_2,use="complete.obs")
cov(ldata$a_2,ldata$l_3,use="complete.obs")
cov(ldata$a_3,ldata$l_4,use="complete.obs")
# Convert numeric censoring variables to 'censored' variable for ltmle
ldata$c_0 <- BinaryToCensoring(is.censored=ldata$c_0)
ldata$c_1 <- BinaryToCensoring(is.censored=ldata$c_1)
ldata$c_2 <- BinaryToCensoring(is.censored=ldata$c_2)
ldata$c_3 <- BinaryToCensoring(is.censored=ldata$c_3)
ldata$c_4 <- BinaryToCensoring(is.censored=ldata$c_4)
