################################################
### Cross-sectional example                   ##
################################################

### NOTE ###
# This code requires data created in the R-code 'datacreation.R', also stored on this repository #

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

# Create subset of data with only baseline variables BA, BB, BC, A0, L0, Y0
csdata <- ldata[,c(2,3,4,6,7,8)]

# Manual TMLE
# Correctly specified
Q0 <- glm(data=csdata,"y_0 ~ a_0 + l_0 + ba + bb + bc")
QAW <- data.frame(cbind(QA=predict(Q0,type="response"),
                        Q0=predict(Q0,type="response",newdata=cbind(csdata[,1:4],a_0=0)),
                        Q1=predict(Q0,type="response",newdata=cbind(csdata[,1:4],a_0=1))))
G <- glm(data=csdata,"a_0 ~ l_0 + ba + bb + bc",family=binomial)
GAW <- predict(G,type="response")
HA1 <- csdata[,5]/GAW
HA0 <- -(1-csdata[,5])/(1-GAW)
H <- HA1+HA0
Q1 <- glm(data=data.frame(cbind(Y=csdata[,6],HA1=HA1,HA0=-HA0,QAW)),"Y ~ -1 + HA1 + HA0 + offset(QA)")
muA1 <- QAW$Q1 + coef(Q1)[1]/GAW
muA0 <- QAW$Q0 + coef(Q1)[2]/(1-GAW)
TMLE <- c(coef=mean(muA1-muA0),
          se=var((HA1-HA0)*(csdata[,6]-QAW$QA) + QAW$Q1 - QAW$Q0 - (muA1-muA0))/length(data[,1]))

# Outcome model mispecified
Q0m1 <- glm(data=csdata,"y_0 ~ a_0 + ba + bb + bc")
QAWm1 <- data.frame(cbind(QA=predict(Q0m1,type="response"),
                          Q0=predict(Q0m1,type="response",newdata=cbind(csdata[,1:4],a_0=0)),
                          Q1=predict(Q0m1,type="response",newdata=cbind(csdata[,1:4],a_0=1))))
Gm1 <- glm(data=csdata,"a_0 ~ l_0 + ba + bb + bc",family=binomial)
GAWm1 <- predict(Gm1,type="response")
HA1m1 <- csdata[,5]/GAWm1
HA0m1 <- -(1-csdata[,5])/(1-GAWm1)
Hm1 <- HA1m1+HA0m1
Q1m1 <- glm(data=data.frame(cbind(Y=csdata[,6],HA1=HA1m1,HA0=-HA0m1,QAWm1)),"Y ~ -1 + HA1 + HA0 + offset(QA)")
muA1m1 <- QAWm1$Q1 + coef(Q1m1)[1]/GAWm1
muA0m1 <- QAWm1$Q0 + coef(Q1m1)[2]/(1-GAWm1)
TMLEm1 <- c(coef=mean(muA1m1-muA0m1),
            se=var((HA1m1-HA0m1)*(csdata[,6]-QAWm1$QA) + QAWm1$Q1 - QAWm1$Q0 - (muA1m1-muA0m1))/length(data[,1]))

# Both models mispecified
Q0m2 <- glm(data=csdata,"y_0 ~ a_0 + ba + bb + bc")
QAWm2 <- data.frame(cbind(QA=predict(Q0m2,type="response"),
                          Q0=predict(Q0m2,type="response",newdata=cbind(csdata[,1:4],a_0=0)),
                          Q1=predict(Q0m2,type="response",newdata=cbind(csdata[,1:4],a_0=1))))
Gm2 <- glm(data=csdata,"a_0 ~ ba + bb + bc",family=binomial)
GAWm2 <- predict(Gm2,type="response")
HA1m2 <- csdata[,5]/GAWm2
HA0m2 <- -(1-csdata[,5])/(1-GAWm2)
Hm2 <- HA1m2+HA0m2
Q1m2 <- glm(data=data.frame(cbind(Y=csdata[,6],HA1=HA1m2,HA0=-HA0m2,QAWm2)),"Y ~ -1 + HA1 + HA0 + offset(QA)")
muA1m2 <- QAWm2$Q1 + coef(Q1m2)[1]/GAWm2
muA0m2 <- QAWm2$Q0 + coef(Q1m2)[2]/(1-GAWm2)
TMLEm2 <- c(coef=mean(muA1m2-muA0m2),
            se=var((HA1m2-HA0m2)*(csdata[,6]-QAWm2$QA) + QAWm2$Q1 - QAWm2$Q0 - (muA1m2-muA0m2))/length(data[,1]))

# GLM-based R-TMLE
# Correctly specified
rtmle <- tmle(Y=csdata[,6],A=csdata[,5],W=csdata[,1:4],
              Q.SL.library=SLlib,
              g.SL.library=SLlib,
              Qform="Y~A+l_0+ba+bb+bc",
              gform="A~l_0+ba+bb+bc")
# Outcome model mispecified
rtmlem1 <- tmle(Y=csdata[,6],A=csdata[,5],W=csdata[,1:4],
                Q.SL.library=SLlib,
                g.SL.library=SLlib,
                Qform="Y~A+ba+bb+bc",
                gform="A~l_0+ba+bb+bc")
# Both models mispecified
rtmlem2 <- tmle(Y=csdata[,6],A=csdata[,5],W=csdata[,1:4],
                Q.SL.library=SLlib,
                g.SL.library=SLlib,
                Qform="Y~A+ba+bb+bc",
                gform="A~ba+bb+bc")

# SuperLearner-based LTMLE
sltmle <- tmle(Y=csdata[,6],A=csdata[,5],W=csdata[,1:4],
               Q.SL.library=SLlib2,
               g.SL.library=SLlib2)

# Create summary table of coefficients and standard errors
csresults <- matrix(c(coef(summary(Q0))[2,1],coef(summary(Q0))[2,2],
                      coef(summary(Q0m1))[2,1],coef(summary(Q0m1))[2,2],
                      mean(muA1-muA0),sqrt(var((HA1-HA0)*(csdata[,6]-QAW$QA) + QAW$Q1 - QAW$Q0 - (muA1-muA0))/length(data[,1])),
                      mean(muA1m1-muA0m1),sqrt(var((HA1m1-HA0m1)*(csdata[,6]-QAWm1$QA) + QAWm1$Q1 - QAWm1$Q0 - (muA1m1-muA0m1))/length(data[,1])),
                      mean(muA1m2-muA0m2),sqrt(var((HA1m2-HA0m2)*(csdata[,6]-QAWm2$QA) + QAWm2$Q1 - QAWm2$Q0 - (muA1m2-muA0m2))/length(data[,1])),
                      rtmle$estimates$ATE$psi,sqrt(rtmle$estimates$ATE$var.psi),
                      rtmlem1$estimates$ATE$psi,sqrt(rtmlem1$estimates$ATE$var.psi),
                      rtmlem2$estimates$ATE$psi,sqrt(rtmlem2$estimates$ATE$var.psi),
                      sltmle$estimates$ATE$psi,sqrt(sltmle$estimates$ATE$var.psi)),nrow=9,ncol=2,byrow=TRUE)
