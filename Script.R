library(BWGS)
data(inra)

#GPLUB with RMR (reduct.marker.size = 5000)

testGBLUPRMR5000 <- bwgs.cv(TRAIN47K,
                            YieldBLUE,
                            # random.pop.size="NULL",
                            # geno.reduct.method="RMR",
                            # reduct.size=5000,
                            geno.impute.method="MNI",
                            predict.method="GBLUP",
                            nFolds=10,
                            nTimes=50)

#RKHS with marker selection by ANOVA (pval=0.001)
testRKHSANO <- bwgs.cv(TRAIN47K,
                       YieldBLUE,
                       geno.reduct.method="ANO",
                       pval=0.001,
                       geno.impute.method="MNI",
                       predict.method="RKHS",
                       nFolds=10,
                       nTimes=50) #

#Boxplot to compare prediction methods
boxplot(cbind(testGBLUPRMR5000$cv, testRKHSANO$cv))

#  sampling training population 
testBWGSRPS300=bwgs.cv(geno_shrink0001,
                       pheno,
                       MAXNA=0.2,
                       MAF=0.05,
                       pop.reduct.method="RANDOM",
                       sample.pop.size=300,
                       geno.reduct.method="NULL",
                       reduct.marker.size="NULL",
                       pval="NULL",
                       r2="NULL",
                       MAP="NULL",
                       geno.impute.method="MNI",
                       predict.method="GBLUP",
                       nFolds=10,
                       nTimes=50)


# optimizing a subset for training
testBWGSOPT300=bwgs.cv(geno_shrink0001,
                       pheno,
                       MAXNA=0.2,
                       MAF=0.05,
                       pop.reduct.method="OPTI",
                       sample.pop.size=300,
                       geno.reduct.method="NULL",
                       reduct.marker.size="NULL",
                       pval="NULL",
                       r2="NULL",
                       MAP="NULL",
                       geno.impute.method="MNI",
                       predict.method="GBLUP",
                       nFolds=10,
                       nTimes=5)

library(tidyverse)

testGBLUPRMR5000$bv_table %>%
  as.data.frame() %>%
  dplyr::arrange(-bv_predict_mean)

testGBLUPRMR5000$summary

testGBLUPRMR5000

predict_rf <- bwgs.predict(geno_train = TRAIN47K,
                              pheno_train = YieldBLUE,
                              geno_target = TARGET47K,
                              MAXNA = 0.2, 
                              MAF = 0.05,
                              geno.reduct.method = "NULL",
                              reduct.size = "NULL",
                              r2 = "NULL",
                              pval = "NULL",
                              MAP = "NULL",
                              geno.impute.method = "MNI",
                              predict.method = "RF")

predict_gblup %>%
  as.data.frame() %>%
  dplyr::arrange(-gpred)

library(rrBLUP)

#random population of 200 lines with 1000 markers
M <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  M[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

#random phenotypes
u <- rnorm(1000)
g <- as.vector(crossprod(t(M),u))
h2 <- 0.5  #heritability
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))

#predict marker effects
ans <- mixed.solve(y,Z=M)  #By default K = I
accuracy <- cor(u,ans$u)

#predict breeding values
ans <- mixed.solve(y,K=A.mat(M))
accuracy <- cor(g,ans$u)

library(G2P)
library(randomForest)

data(GYSS)
########## predicting breeding value
predRes <- G2P(Markers[1:200,],phenotype[1:200],Markers[201:242,],
               phenotype[201:242],modelMethods = c("rrBLUP", "RFC"),
               outputModel = FALSE)

library(SpATS)
data(wheatdata)
summary(wheatdata)

# Create factor variable for row and columns
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

m0 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20), degree = 3, pord = 2), 
            genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, data = wheatdata, 
            control =  list(tolerance = 1e-03))
# Brief summary
m0
# More information: dimensions
summary(m0) # summary(fit.m2, which = "dimensions") 
# More information: variances
summary(m0, which = "variances") 
# More information: all
summary(m0, which = "all") 

# Plot results
plot(m0)
plot(m0, all.in.one = FALSE)

# Variogram
var.m0 <- variogram(m0)
plot(var.m0)


