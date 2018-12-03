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
qforma <- c(l_0="Q.kplus1 ~ ba + bb + bc",
            l_1="Q.kplus1 ~ ba + bb + bc + l_0 + a_0",
            l_2="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + l_1 + a_1",
            l_3="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + l_1 + a_1 + l_2 + a_2",
            l_4="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + l_1 + a_1 + l_2 + a_2 + l_3 + a_3",
            y_4="Q.kplus1 ~ ba + bb + bc + l_0 + a_0 + l_1 + a_1 + l_2 + a_2 + l_3 + a_3 + l_4 + a_4")
mqforma <- c(l_0="Q.kplus1 ~ ba + bb + bc",
             l_1="Q.kplus1 ~ ba + bb + bc + a_0",
             l_2="Q.kplus1 ~ ba + bb + bc + a_0 + a_1",
             l_3="Q.kplus1 ~ ba + bb + bc + a_0 + a_1 + a_2",
             l_4="Q.kplus1 ~ ba + bb + bc + a_0 + a_1 + a_2 + a_3",
             y_4="Q.kplus1 ~ ba + bb + bc + a_0 + a_1 + a_2 + a_3 + a_4")

gforma <- c(c_0="c_0 ~ ba + bb + bc",
            a_0="a_0 ~ ba + bb + bc + l_0",
            c_1="c_1 ~ ba + bb + bc + l_0 + a_0 ",
            a_1="a_1 ~ ba + bb + bc + l_0 + a_0 + l_1",
            c_2="c_2 ~ ba + bb + bc + l_0 + a_0 + l_1 + a_1",
            a_2="a_2 ~ ba + bb + bc + l_0 + a_0 + l_1 + a_1 + l_2",
            c_3="c_3 ~ ba + bb + bc + l_0 + a_0 + l_1 + a_1 + l_2 + a_2",
            a_3="a_3 ~ ba + bb + bc + l_0 + a_0 + l_1 + a_1 + l_2 + a_2 + l_3",
            c_4="c_4 ~ ba + bb + bc + l_0 + a_0 + l_1 + a_1 + l_2 + a_2 + l_3 + a_3",
            a_4="a_4 ~ ba + bb + bc + l_0 + a_0 + l_1 + a_1 + l_2 + a_2 + l_3 + a_3 + l_4")
mgforma <- c(c_0="c_0 ~ ba + bb + bc",
             a_0="a_0 ~ ba + bb + bc",
             c_1="c_1 ~ ba + bb + bc + a_0",
             a_1="a_1 ~ ba + bb + bc + a_0",
             c_2="c_2 ~ ba + bb + bc + a_0 + a_1",
             a_2="a_2 ~ ba + bb + bc + a_0 + a_1",
             c_3="c_3 ~ ba + bb + bc + a_0 + a_1 + a_2",
             a_3="a_3 ~ ba + bb + bc + a_0 + a_1 + a_2",
             c_4="c_4 ~ ba + bb + bc + a_0 + a_1 + a_2 + a_3",
             a_4="a_4 ~ ba + bb + bc + a_0 + a_1 + a_2 + a_3")

# Create data subset with all observations of exposure and confounders, but only final outcome Y4
ldata <- data[,c(-1,-8,-12,-16,-20)]

# NAIVE ANALYSIS #
# Correctly specified
LGLM <- glm(data=data,"y_4 ~ a_0 + a_1 + a_2 + a_3 + a_4 + l_0 + l_1 + l_2 + l_3 + l_4 + ba + bb + bc")
V1<-vcov(LGLM) # Save variance-covariance matrix to calculate joint standard error
# Incorrectly specified
LGLMm <- glm(data=data,"y_4 ~ a_0 + a_1 + a_2 + a_3 + a_4 + ba + bb + bc")
V2<-vcov(LGLMm) # Save variance-covariance matrix to calculate joint standard error

# TMLE ANALYSIS #
# Estimation using just GLMs
# Correctly specified
rltmle1 <- suppressWarnings(ltmle(ldata,
                                  Anodes=c("a_0","a_1","a_2","a_3","a_4"),
                                  Lnodes=c("l_0","l_1","l_2","l_3","l_4"),
                                  Cnodes=c("c_0","c_1","c_2","c_3","c_4"),
                                  Ynodes="y_4",
                                  abar=list(c(1,1,1,1,1),c(0,0,0,0,0)),
                                  SL.library=SLlib,
                                  Qform=qforma,gform=gforma,
                                  estimate.time=FALSE,
                                  survivalOutcome=FALSE))

# Outcome model misspecified
rltmle1m1 <- suppressWarnings(ltmle(ldata,
                                    Anodes=c("a_0","a_1","a_2","a_3","a_4"),
                                    Lnodes=c("l_0","l_1","l_2","l_3","l_4"),
                                    Cnodes=c("c_0","c_1","c_2","c_3","c_4"),
                                    Ynodes="y_4",
                                    abar=list(c(1,1,1,1,1),c(0,0,0,0,0)),
                                    SL.library=SLlib,
                                    Qform=mqforma,gform=gforma,
                                    estimate.time=FALSE,
                                    survivalOutcome=FALSE))

# Both models misspecified
rltmle1m2 <- suppressWarnings(ltmle(ldata,
                                    Anodes=c("a_0","a_1","a_2","a_3","a_4"),
                                    Lnodes=c("l_0","l_1","l_2","l_3","l_4"),
                                    Cnodes=c("c_0","c_1","c_2","c_3","c_4"),
                                    Ynodes="y_4",
                                    abar=list(c(1,1,1,1,1),c(0,0,0,0,0)),
                                    SL.library=SLlib,
                                    Qform=mqforma,gform=mgforma,
                                    estimate.time=FALSE,
                                    survivalOutcome=FALSE))

# Estimation via SuperLearner
slltmle1 <- suppressWarnings(ltmle(ldata,
                                   Anodes=c("a_0","a_1","a_2","a_3","a_4"),
                                   Lnodes=c("l_0","l_1","l_2","l_3","l_4"),
                                   Cnodes=c("c_0","c_1","c_2","c_3","c_4"),
                                   Ynodes="y_4",
                                   abar=list(c(1,1,1,1,1),c(0,0,0,0,0)),
                                   SL.library=SLlib3,
                                   estimate.time=FALSE,
                                   survivalOutcome=FALSE))

# Create summary table of coefficients and standard errors
lresults1 <- matrix(c(coef(LGLM)[2]+coef(LGLM)[3]+coef(LGLM)[4]+coef(LGLM)[5]+coef(LGLM)[6],
                      V1[2,2] + V1[3,3] + V1[4,4] + V1[5,5] + V1[6,6]
                      + 2*V1[2,3]+ 2*V1[2,4] + 2*V1[2,5] + 2*V1[2,6]
                      + 2*V1[3,4] + 2*V1[3,5] + 2*V1[3,6]
                      + 2*V1[4,5] + 2*V1[4,6]
                      + 2*V1[5,6],
                      coef(LGLMm)[2]+coef(LGLMm)[3]+coef(LGLMm)[4]+coef(LGLMm)[5]+coef(LGLMm)[6],
                      V2[2,2] + V2[3,3] + V2[4,4] + V2[5,5] + V2[6,6]
                      + 2*V2[2,3]+ 2*V2[2,4] + 2*V2[2,5] + 2*V2[2,6]
                      + 2*V2[3,4] + 2*V2[3,5] + 2*V2[3,6]
                      + 2*V2[4,5] + 2*V2[4,6]
                      + 2*V2[5,6],
                      summary(rltmle1)$effect.measures$ATE$estimate,
                      summary(rltmle1)$effect.measures$ATE$std.dev,
                      summary(rltmle1m1)$effect.measures$ATE$estimate,
                      summary(rltmle1m1)$effect.measures$ATE$std.dev,
                      summary(rltmle1m2)$effect.measures$ATE$estimate,
                      summary(rltmle1m2)$effect.measures$ATE$std.dev,
                      summary(slltmle1)$effect.measures$ATE$estimate,
                      summary(slltmle1)$effect.measures$ATE$std.dev),nrow=6,ncol=2,byrow=TRUE)
