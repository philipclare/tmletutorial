### NOTE ###
This code requires data created in the R-code 'datacreation.R', also stored on this repository

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

# Define SuperLearner libraries to be used by SL-TMLE
SLlib <- c("SL.glm")
SLlib2 <- c("SL.glm","SL.glm.interaction","SL.stepAIC","SL.ranger")
SLlib3 <- list(Q=c("SL.glm","SL.glm.interaction","SL.stepAIC"),
               g=c("SL.glm","SL.glm.interaction","SL.stepAIC","SL.ranger"))

# Define q and g forms for manually specified LTMLE models
qformb <- c(l_0="Q.kplus1 ~ ba + bb + bc",
            y_0="Q.kplus1 ~ ba + bb + bc + l_0 + a_0",
            l_1="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + y_0",
            y_1="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1",
            l_2="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1",
            y_2="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1 + l_2 + a_2",
            l_3="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1 + l_2 + a_2 + y_2",
            y_3="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1 + l_2 + a_2 + y_2 + l_3 + a_3",
            l_4="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1 + l_2 + a_2 + y_2 + l_3 + a_3 + y_3",
            y_4="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1 + l_2 + a_2 + y_2 + l_3 + a_3 + y_3 + l_4 + a_4")
mqformb <- c(l_0="Q.kplus1 ~ ba + bb + bc",
             y_0="Q.kplus1 ~ ba + bb + bc + a_0",
             l_1="Q.kplus1 ~ ba + bb + bc + a_0 + y_0",
             y_1="Q.kplus1 ~ ba + bb + bc + a_0 + y_0 + a_1",
             l_2="Q.kplus1 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1",
             y_2="Q.kplus1 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1 + a_2",
             l_3="Q.kplus1 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1 + a_2 + y_2",
             y_3="Q.kplus1 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1 + a_2 + y_2 + a_3",
             l_4="Q.kplus1 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1 + a_2 + y_2 + a_3 + y_3",
             y_4="Q.kplus1 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1 + a_2 + y_2 + a_3 + y_3 + a_4")

gformb <- c(c_0="c_0 ~ ba + bb + bc",
            a_0="a_0 ~ ba + bb + bc + l_0",
            c_1="c_1 ~ ba + bb + bc + l_0 + a_0 + y_0",
            a_1="a_1 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1",
            c_2="c_2 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1",
            a_2="a_2 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1 + l_2",
            c_3="c_3 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1 + l_2 + a_2 + y_2",
            a_3="a_3 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1 + l_2 + a_2 + y_2 + l_3",
            c_4="c_4 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1 + l_2 + a_2 + y_2 + l_3 + a_3 + y_3",
            a_4="a_4 ~ ba + bb + bc + l_0 + a_0 + y_0 + l_1 + a_1 + y_1 + l_2 + a_2 + y_2 + l_3 + a_3 + y_3 + l_4")
mgformb <- c(c_0="c_0 ~ ba + bb + bc",
             a_0="a_0 ~ ba + bb + bc",
             c_1="c_1 ~ ba + bb + bc + a_0 + y_0",
             a_1="a_1 ~ ba + bb + bc + a_0 + y_0",
             c_2="c_2 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1",
             a_2="a_2 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1",
             c_3="c_3 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1 + a_2 + y_2",
             a_3="a_3 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1 + a_2 + y_2",
             c_4="c_4 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1 + a_2 + y_2 + a_3 + y_3",
             a_4="a_4 ~ ba + bb + bc + a_0 + y_0 + a_1 + y_1 + a_2 + y_2 + a_3 + y_3")

# TMLE ANALYSIS #
# Estimation using just GLMs
# Correctly specified
rltmle2 <- suppressWarnings(ltmle(data[,-1],
                                  Anodes=c("a_0","a_1","a_2","a_3","a_4"),
                                  Lnodes=c("l_0","l_1","l_2","l_3","l_4"),
                                  Cnodes=c("c_0","c_1","c_2","c_3","c_4"),
                                  Ynodes=c("y_0","y_1","y_2","y_3","y_4"),
                                  survivalOutcome=FALSE,
                                  abar=list(c(1,1,1,1,1),c(0,0,0,0,0)),
                                  SL.library=SLlib,
                                  Qform=qformb,gform=gformb,
                                  estimate.time=FALSE))
# Outcome model misspecified
rltmle2m1 <- suppressWarnings(ltmle(data[,-1],
                                    Anodes=c("a_0","a_1","a_2","a_3","a_4"),
                                    Lnodes=c("l_0","l_1","l_2","l_3","l_4"),
                                    Cnodes=c("c_0","c_1","c_2","c_3","c_4"),
                                    Ynodes=c("y_0","y_1","y_2","y_3","y_4"),
                                    survivalOutcome=FALSE,
                                    abar=list(c(1,1,1,1,1),c(0,0,0,0,0)),
                                    SL.library=SLlib,
                                    Qform=mqformb,gform=gformb,
                                    estimate.time=FALSE))
# Both models misspecified
rltmle2m2 <- suppressWarnings(ltmle(data[,-1],
                                    Anodes=c("a_0","a_1","a_2","a_3","a_4"),
                                    Lnodes=c("l_0","l_1","l_2","l_3","l_4"),
                                    Cnodes=c("c_0","c_1","c_2","c_3","c_4"),
                                    Ynodes=c("y_0","y_1","y_2","y_3","y_4"),
                                    survivalOutcome=FALSE,
                                    abar=list(c(1,1,1,1,1),c(0,0,0,0,0)),
                                    SL.library=SLlib,
                                    Qform=mqformb,gform=mgformb,
                                    estimate.time=FALSE))

# Estimation using SuperLearner
slltmle2 <- suppressWarnings(ltmle(data[,-1],
                                   Anodes=c("a_0","a_1","a_2","a_3","a_4"),
                                   Lnodes=c("l_0","l_1","l_2","l_3","l_4"),
                                   Cnodes=c("c_0","c_1","c_2","c_3","c_4"),
                                   Ynodes=c("y_0","y_1","y_2","y_3","y_4"),
                                   survivalOutcome=FALSE,
                                   abar=list(c(1,1,1,1,1),c(0,0,0,0,0)),
                                   SL.library=SLlib3,
                                   estimate.time=FALSE))

# NAIVE ANALYSIS #
# Create cumulative exposure and confounder variables for naive analysis
data$acum_0 <- data$a_0
data$acum_1 <- data$acum_0 + data$a_1
data$acum_2 <- data$acum_1 + data$a_2
data$acum_3 <- data$acum_2 + data$a_3
data$acum_4 <- data$acum_3 + data$a_4
data$lcum_0 <- data$l_0
data$lcum_1 <- data$lcum_0 + data$l_1
data$lcum_2 <- data$lcum_1 + data$l_2
data$lcum_3 <- data$lcum_2 + data$l_3
data$lcum_4 <- data$lcum_3 + data$l_4

# Create long-form data for GLM analysis
longdata <- reshape(data=data,varying=list(c(8,12,16,20,24),c(7,11,15,19,23),c(25:29),c(6,10,14,18,22),c(30:34),c(5,9,13,17,21)),
                    v.names=c("y","a","acum","l","lcum","c"),
                    idvar="ID",timevar="obs",direction="long", sep = "_")
# Correctly specified cumulative exposure model using random intercept model
LRI <- lmer("y ~ acum + lcum + ba + bb + bc + (1|ID)",data=longdata)
# Incorrectly specified cumulative exposure model
LRIm <- lmer("y ~ acum + ba + bb + bc + (1|ID)",data=longdata)

# Create summary table of coefficients and standard errors
lresults2 <- matrix(c(5*coef(summary(LRI))[2,1],
                      5*coef(summary(LRI))[2,2],
                      5*coef(summary(LRIm))[2,1],
                      5*coef(summary(LRIm))[2,2],
                      summary(rltmle2)$effect.measures$ATE$estimate,
                      summary(rltmle2)$effect.measures$ATE$std.dev,
                      summary(rltmle2m1)$effect.measures$ATE$estimate,
                      summary(rltmle2m1)$effect.measures$ATE$std.dev,
                      summary(rltmle2m2)$effect.measures$ATE$estimate,
                      summary(rltmle2m2)$effect.measures$ATE$std.dev,
                      summary(slltmle2)$effect.measures$ATE$estimate,
                      summary(slltmle2)$effect.measures$ATE$std.dev),nrow=6,ncol=2,byrow=TRUE)
