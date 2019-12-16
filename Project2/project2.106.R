######################################################
######       STA 106 / Project1 / Liya Li       ######
######                 Appendix                 ######
######################################################
library(ggplot2) # for ggplots
library(car)     # for BF.test
library(EnvStats)# for boxcox()

######               Topic I - Q2               ######
# The goal is to assess if these hawks can be easily distinguished by the length of their wing feathers.
NewHawk <- read.csv("~/Desktop/Winter Quarter 2018/STA 106/106_project2/NewHawk.csv")
View(NewHawk)

######################################################
###### explore the data
total1 <- dim(NewHawk)[1] # nt
total1 
str(NewHawk)              # variables
a1 = length(unique(NewHawk$Species)) 
a1                        # a

# number of observations, mean, sd, var in each group
n.i1 <- aggregate(Wing ~ Species, data = NewHawk, length)[,2]
mu.i1 <- aggregate(Wing ~ Species, data = NewHawk, mean)[,2]
sd.i1 <- aggregate(Wing ~ Species, data = NewHawk, sd)[,2]
var.i1 <- aggregate(Wing ~ Species, data = NewHawk, var)[,2]

dataframe1 <- data.frame(groups=c("Cooper's","Red-tailed","Sharp-Shinned"), n.i1, mu.i1, sd.i1, var.i1)
dataframe1

# plot the original data
plot(NewHawk$Wing, NewHawk$Species, xlab = "Wing length(mm)", ylab = "Species (1.0:CH, 2.0:RT, 3.0:SS)", main = "Relationship between Wing length and Species")

######################################################
###### check ANOVA assumptions on original data
model1 <- lm(Wing~Species, data = NewHawk) # fit the model
model1
par(mfrow = c(2,2))
plot(model1)

# add residuals column to the original dataset and view the dataset 
NewHawk$ei <- model1$residuals
View(NewHawk)

# use tests to check normality and equal variance assumptions on original data
the.shapiro.pval1 <- shapiro.test(model1$residuals)$p.val
the.shapiro.pval1
the.BFtest1 <- leveneTest(ei ~ Species, data=NewHawk, center=median)
the.BFtest1
the.BF.pval1 <- the.BFtest1[[3]][1]
the.BF.pval1

######################################################
###### check ANOVA assumptions on data without outliers
# find the outliers
# outliers via boxplot
par(mfrow = c(1,1))
boxplot(Wing ~ factor(Species), data = NewHawk, xlab = "Species (CH,RT,SS)", ylab = "Wing Length (mm) ", main = "Boxplot: Wing Length vs Species")

# outliers via semi-studentized/standardized residuals
# find SSE, MSE, eij.star
SSE1 <- sum(NewHawk$ei^2)
SSE1
MSE1 <- SSE1 / (total1-a1)
MSE1  
eij.star1 <- model1$residuals / sqrt(MSE1)

# find the outliers
alpha1 <- 0.05
t.cutoff1 <- qt(1-alpha1/(2*total1), total1-a1)
t.cutoff1
CO.eij1 <- which(abs(eij.star1) > t.cutoff1)
CO.eij1

# outliers via studentized/standardized residuals
rij1 = rstandard(model1)
CO.rij1 = which(abs(rij1) > t.cutoff1)
CO.rij1

# remove the outliers, fit new model, check assumptions
outliers1 <- CO.rij1
outliers.percentage1 <- length(outliers1)/total1
outliers.percentage1
NewHawk2 <- NewHawk[-outliers1,]
View(NewHawk2)
model2 <- lm(Wing ~ Species, data = NewHawk2) # fit a model with data without outliers
model2
par(mfrow = c(2,2))
plot(model2)

# add residuals column to the no outliers dataset and view the dataset
NewHawk2$ei <- model2$residuals
View(NewHawk2)

# use tests to check normality and equal variance assumptions on data without outliers
the.shapiro.pval2 = shapiro.test(model2$residuals)$p.val
the.shapiro.pval2
the.BFtest2 = leveneTest(ei ~ Species, data=NewHawk2, center=median)
the.BFtest2
the.BF.pval2 = the.BFtest2[[3]][1]
the.BF.pval2

######################################################
### check ANOVA assumptions on transformed original data
# see if what transformation we can do on original data
# three different ways to find lambda
L1.model1 <- boxcox(model1 ,objective.name = "PPCC",optimize = TRUE)$lambda
L2.model1 <- boxcox(model1 ,objective.name = "Shapiro-Wilk",optimize = TRUE)$lambda
L3.model1 <- boxcox(NewHawk$Wing, objective.name = "Log-Likelihood",optimize = TRUE)$lambda

L.model1 <- c(L1.model1, L2.model1, L3.model1)
L.model1

# L1.model1
NewHawk$Wing.trans3.1 <- NewHawk$Wing^L1.model1
model3.1 <- lm(Wing.trans3.1 ~ Species, data = NewHawk)
NewHawk$ei3.1 <- model3.1$residuals
par(mfrow = c(2,2))
plot(model3.1, main = "Diagnostics Plots for L1.model1")

the.shapiro.pval3.1 = shapiro.test(model3.1$residuals)$p.val
the.shapiro.pval3.1
the.BFtest3.1 = leveneTest(ei3.1 ~ Species, data=NewHawk, center=median)
the.BFtest3.1
the.BF.pval3.1 = the.BFtest3.1[[3]][1]
the.BF.pval3.1

# L2.model1
NewHawk$Wing.trans3.2 <- NewHawk$Wing^L2.model1
model3.2 <- lm(Wing.trans3.2 ~ Species, data = NewHawk)
NewHawk$ei3.2 <- model3.2$residuals
par(mfrow = c(2,2))
plot(model3.2, main = "Diagnostics Plots for L2.model1")

the.shapiro.pval3.2 = shapiro.test(model3.2$residuals)$p.val
the.shapiro.pval3.2
the.BFtest3.2 = leveneTest(ei3.2 ~ Species, data=NewHawk, center=median)
the.BFtest3.2
the.BF.pval3.2 = the.BFtest3.2[[3]][1]
the.BF.pval3.2

# L3.model3
NewHawk$Wing.trans3.3 <- NewHawk$Wing^L3.model1
model3.3 <- lm(Wing.trans3.3 ~ Species, data = NewHawk)
NewHawk$ei3.3 <- model3.3$residuals
par(mfrow = c(2,2))
plot(model3.3, main = "Diagnostics Plots for L3.model1")

the.shapiro.pval3.3 = shapiro.test(model3.3$residuals)$p.val
the.shapiro.pval3.3
the.BFtest3.3 = leveneTest(ei3.3 ~ Species, data=NewHawk, center=median)
the.BFtest3.3
the.BF.pval3.3 = the.BFtest3.3[[3]][1]
the.BF.pval3.3

######################################################
### check ANOVA assumptions on transformed no outliers data
# see if what transformation we can do on no outliers data
# three different ways to find lambda
L1.model2 <- boxcox(model2 ,objective.name = "PPCC",optimize = TRUE)$lambda
L2.model2 <- boxcox(model2 ,objective.name = "Shapiro-Wilk",optimize = TRUE)$lambda
L3.model2 <- boxcox(NewHawk2$Wing, objective.name = "Log-Likelihood",optimize = TRUE)$lambda

L.model2 <- c(L1.model2, L2.model2, L3.model2)
L.model2

# L1.model2
NewHawk2$Wing.trans4.1 <- NewHawk2$Wing^L1.model2
model4.1 <- lm(Wing.trans4.1 ~ Species, data = NewHawk2)
NewHawk2$ei4.1 <- model4.1$residuals
par(mfrow = c(2,2))
plot(model4.1, main = "Diagnostics Plots for L1.model2")

the.shapiro.pval4.1 = shapiro.test(model4.1$residuals)$p.val
the.shapiro.pval4.1
the.BFtest4.1 = leveneTest(ei4.1 ~ Species, data=NewHawk2, center=median)
the.BFtest4.1
the.BF.pval4.1 = the.BFtest4.1[[3]][1]
the.BF.pval4.1

# L2.model2
NewHawk2$Wing.trans4.2 <- NewHawk2$Wing^L2.model2
model4.2 <- lm(Wing.trans4.2 ~ Species, data = NewHawk2)
NewHawk2$ei4.2 <- model4.2$residuals
par(mfrow = c(2,2))
plot(model4.2, main = "Diagnostics Plots for L2.model2")

the.shapiro.pval4.2 = shapiro.test(model4.2$residuals)$p.val
the.shapiro.pval4.2
the.BFtest4.2 = leveneTest(ei4.2 ~ Species, data=NewHawk2, center=median)
the.BFtest4.2
the.BF.pval4.2 = the.BFtest4.2[[3]][1]
the.BF.pval4.2

# L3.model3
NewHawk2$Wing.trans4.3 <- NewHawk2$Wing^L3.model2
model4.3 <- lm(Wing.trans4.3 ~ Species, data = NewHawk2)
NewHawk2$ei4.3 <- model4.3$residuals
par(mfrow = c(2,2))
plot(model4.3, main = "Diagnostics Plots for L3.model2")

the.shapiro.pval4.3 = shapiro.test(model4.3$residuals)$p.val
the.shapiro.pval4.3
the.BFtest4.3 = leveneTest(ei4.3 ~ Species, data=NewHawk2, center=median)
the.BFtest4.3
the.BF.pval4.3 = the.BFtest4.3[[3]][1]
the.BF.pval4.3

### compare the ANOVA assumptions tests results of original data and no outliers data
model1.ANOVA <- c(the.shapiro.pval1,the.BF.pval1,0.01)
model2.ANOVA <- c(the.shapiro.pval2,the.BF.pval2,0.01)
model3.1.ANOVA <- c(the.shapiro.pval3.1,the.BF.pval3.1,0.01)
model3.2.ANOVA <- c(the.shapiro.pval3.2,the.BF.pval3.2,0.01)
model3.3.ANOVA <- c(the.shapiro.pval3.3,the.BF.pval3.3,0.01)
model4.1.ANOVA <- c(the.shapiro.pval4.1,the.BF.pval4.1,0.01)
model4.2.ANOVA <- c(the.shapiro.pval4.2,the.BF.pval4.2,0.01)
model4.3.ANOVA <- c(the.shapiro.pval4.3,the.BF.pval4.3,0.01)
rownames <- c("Shapiro.pval","BF.pval","Alpha")
df.test1 <- data.frame(rownames, model1.ANOVA, model2.ANOVA, model3.1.ANOVA, model3.2.ANOVA, model3.3.ANOVA, model4.1.ANOVA, model4.2.ANOVA, model4.3.ANOVA)
df.test1

######               Topic II - Q1             ######
### read in the data
Salary <- read.csv("~/Desktop/Winter Quarter 2018/STA 106/106_project2/Salary.csv")
head(Salary)

######################################################
### explore the data
total2 <- nrow(Salary)    # nt
total2 
str(Salary)               # variables
a2 <- length(unique(Salary[,2])) 
a2                        # a
b2 <- length(unique(Salary[,3]))
b2                        # b

names(Salary) <- c("Y", "A", "B")
View(Salary)


### find all the means
find.means = function(the.data,fun.name = mean){
  a = length(unique(the.data[,2]))
  b = length(unique(the.data[,3]))
  means.A = by(the.data[,1], the.data[,2], fun.name)
  means.B = by(the.data[,1],the.data[,3],fun.name)
  means.AB = by(the.data[,1],list(the.data[,2],the.data[,3]),fun.name)
  MAB = matrix(means.AB,nrow = b, ncol = a, byrow = TRUE)
  colnames(MAB) = names(means.A)
  rownames(MAB) = names(means.B)
  MA = as.numeric(means.A)
  names(MA) = names(means.A)
  MB = as.numeric(means.B)
  names(MB) = names(means.B)
  MAB = t(MAB)
  results = list(A = MA, B = MB, AB = MAB)
  return(results)
}

means <- find.means(Salary, fun.name = mean)

### find all the sds
find.sds = function(the.data,fun.name = sd){
  a = length(unique(the.data[,2]))
  b = length(unique(the.data[,3]))
  sds.A = by(the.data[,1], the.data[,2], fun.name)
  sds.B = by(the.data[,1],the.data[,3],fun.name)
  sds.AB = by(the.data[,1],list(the.data[,2],the.data[,3]),fun.name)
  MAB = matrix(sds.AB,nrow = b, ncol = a, byrow = TRUE)
  colnames(MAB) = names(sds.A)
  rownames(MAB) = names(sds.B)
  MA = as.numeric(sds.A)
  names(MA) = names(sds.A)
  MB = as.numeric(sds.B)
  names(MB) = names(sds.B)
  MAB = t(MAB)
  results = list(A = MA, B = MB, AB = MAB)
  return(results)
}

sds <- find.sds(Salary, fun.name = sd)

### find all the ns
find.ns = function(the.data,fun.name = length){
  a = length(unique(the.data[,2]))
  b = length(unique(the.data[,3]))
  ns.A = by(the.data[,1], the.data[,2], fun.name)
  ns.B = by(the.data[,1],the.data[,3],fun.name)
  ns.AB = by(the.data[,1],list(the.data[,2],the.data[,3]),fun.name)
  MAB = matrix(ns.AB,nrow = b, ncol = a, byrow = TRUE)
  colnames(MAB) = names(ns.A)
  rownames(MAB) = names(ns.B)
  MA = as.numeric(ns.A)
  names(MA) = names(ns.A)
  MB = as.numeric(ns.B)
  names(MB) = names(ns.B)
  MAB = t(MAB)
  results = list(A = MA, B = MB, AB = MAB)
  return(results)
}

ns <- find.ns(Salary, fun.name = length)

means
sds
ns


### check interaction
interaction.plot(Salary$A, Salary$B, Salary$Y, main = "Interaction Plot")

### boxplot
par(mfrow = c(1,1))
boxplot(Y ~ A*B, data = Salary)

### fit full model
AB <- lm(Y ~ A*B, Salary)

### check outliers using the full model
Salary$ei <- AB$residuals
SSE2 <- sum(Salary$ei^2)
MSE2 <- SSE2 / (total2-a2*b2)
eij.star2 <- AB$residuals / sqrt(MSE2)
alpha2 <- 0.05
t.cutoff2 <- qt(1-alpha2/(2*total2), total2-a2*b2)
t.cutoff2
CO.eij2 <- which(abs(eij.star2) > t.cutoff2)
CO.eij2

### diagnostic
par(mfrow = c(2,2))
plot(AB, main = "Diagnostics Plot for Interaction Model")

### fit other models
A.B <- lm(Y ~ A+B, Salary)
A <- lm(Y ~ A, Salary)
B <- lm(Y ~ B, Salary)
N <- lm(Y ~ 1, Salary)

### check outliers using the A+B model
Salary$eiA.B <- A.B$residuals
SSE3 <- sum(Salary$eiA.B^2)
MSE3 <- SSE3 / (total2-a2*b2)
eij.star3 <- A.B$residuals / sqrt(MSE3)
alpha2 <- 0.05
t.cutoff3 <- qt(1-alpha2/(2*total2), total2-a2*b2)
t.cutoff3
CO.eij3 <- which(abs(eij.star3) > t.cutoff3)
CO.eij3

### find SSEs
all.models <- list(AB,A.B,A,B,N)
SSE <- t(as.matrix(sapply(all.models,function(M) sum(M$residuals^2))))
colnames(SSE) <- c("AB","(A+B)","A","B","Empty/Null")
rownames(SSE) <- "SSE"

### Hypothesis Test
# for interaction
resultAB <- anova(A.B, AB)
resultAB
resultAB[2,6]

# do these two when there is no interaction
# for factor A
resultA <- anova(B, A.B)
resultA
resultA[2,6]

# for factor B
resultB <- anova(A, A.B)
resultB
resultB[2,6]

### find partial R^2
Partial.R2 = function(small.model,big.model){
  SSE1 = sum(small.model$residuals^2)
  SSE2 = sum(big.model$residuals^2)
  PR2 = (SSE1 - SSE2)/SSE1
  return(PR2)
}

R2.interaction <- Partial.R2(A.B, AB)
R2.interaction
R2.A <- Partial.R2(B, AB)
R2.A
R2.B <- Partial.R2(A, AB)
R2.B

# use no interaction model: Y ~ A + B, check diagnostics again
### check diagnostic again
Salary$ei <- A.B$residuals
par(mfrow = c(2,2))
plot(A.B, main = "Diagnostics Plot for No Interaction Model")

the.shapiro.pval.A.B <- shapiro.test(A.B$residuals)$p.val
the.shapiro.pval.A.B
the.BFtest.A.B <- leveneTest(ei ~ A*B, data=Salary, center=median)
the.BFtest.A.B
the.BF.pval.A.B <- the.BFtest.A.B[[3]][1]
the.BF.pval.A.B

overall <- sum(Salary[,1])/total2
overall

# equal weights
######################################################
find.mult = function(alpha,a,b,dfSSE,g,group){
  if(group == "A"){
    Tuk = round(qtukey(1-alpha,a,dfSSE)/sqrt(2),3)
    Bon = round(qt(1-alpha/(2*g), dfSSE ) ,3)
    Sch = round(sqrt((a-1)*qf(1-alpha, a-1, dfSSE)),3) 
  }else if(group == "B"){
    Tuk = round(qtukey(1-alpha,b,dfSSE)/sqrt(2),3)
    Bon = round(qt(1-alpha/(2*g), dfSSE ) ,3)
    Sch = round(sqrt((b-1)*qf(1-alpha, b-1, dfSSE)),3) 
  }else if(group == "AB"){
    Tuk = round(qtukey(1-alpha,a*b,dfSSE)/sqrt(2),3)
    Bon = round(qt(1-alpha/(2*g), dfSSE ) ,3)
    Sch = round(sqrt((a*b-1)*qf(1-alpha, a*b-1, dfSSE)),3) 
  }
  results = c(Bon, Tuk,Sch)
  names(results) = c("Bonferroni","Tukey","Scheffe")
  return(results)
}

scary.CI = function(the.data,MSE,equal.weights = TRUE,multiplier,group,cs){
  if(sum(cs) != 0 & sum(cs !=0 ) != 1){
    return("Error - you did not input a valid contrast")
  }else{
    the.means = find.means(the.data)
    the.ns =find.means(the.data,length)
    nt = nrow(the.data)
    a = length(unique(the.data[,2]))
    b = length(unique(the.data[,3]))
    if(group =="A"){
      if(equal.weights == TRUE){
        a.means = rowMeans(the.means$AB)
        est = sum(a.means*cs)
        mul = rowSums(1/the.ns$AB)
        SE = sqrt(MSE/b^2 * (sum(cs^2*mul)))
        N = names(a.means)[cs!=0]
        CS = paste("(",cs[cs!=0],")",sep = "")
        fancy = paste(paste(CS,N,sep =""),collapse = "+")
        names(est) = fancy
      } else{
        a.means = the.means$A
        est = sum(a.means*cs)
        SE = sqrt(sum(MSE/b^2*cs^2*(1/the.ns$A)))
        N = names(a.means)[cs!=0]
        CS = paste("(",cs[cs!=0],")",sep = "")
        fancy = paste(paste(CS,N,sep =""),collapse = "+")
        names(est) = fancy
      }
    }else if(group == "B"){
      if(equal.weights == TRUE){
        b.means = colMeans(the.means$AB)
        est = sum(b.means*cs)
        mul = colSums(1/the.ns$AB)
        SE = sqrt(MSE/a^2 * (sum(cs^2*mul)))
        N = names(b.means)[cs!=0]
        CS = paste("(",cs[cs!=0],")",sep = "")
        fancy = paste(paste(CS,N,sep =""),collapse = "+")
        names(est) = fancy
      } else{
        b.means = the.means$B
        est = sum(b.means*cs)
        SE = sqrt(MSE/a^2*sum(cs^2*(1/the.ns$B)))
        N = names(b.means)[cs!=0]
        CS = paste("(",cs[cs!=0],")",sep = "")
        fancy = paste(paste(CS,N,sep =""),collapse = "+")
        names(est) = fancy
      }
    } else if(group == "AB"){
      est = sum(cs*the.means$AB)
      SE = sqrt(MSE*sum(cs^2/the.ns$AB))
      names(est) = "someAB"
    }
    the.CI = est + c(-1,1)*multiplier*SE
    results = c(est,the.CI)
    names(results) = c(names(est),"lower bound","upper bound")
    return(results)
  }
}

df.A.B <- total2 - a2 - b2 +1
SSE.A.B <- SSE[1,2]
MSE.A.B <- SSE.A.B / df.A.B
alpha.A.B <- 0.05

### four pairwise 95% CIs
# u11 - u12
multiplier1 <- find.mult(alpha = alpha.A.B, a = a2, b = b2, dfSSE = df.A.B, g = 1, group = "AB")
mul1 <- multiplier1[["Tukey"]]
c1 <- matrix(0, nrow = a2, ncol = b2)
c1[1,1] <- 1
c1[1,2] <- -1
CI1 <- scary.CI(Salary, MSE.A.B, equal.weights = T, mul1, "AB", c1)

# u1. - u2.
multiplier2 <- find.mult(alpha = alpha.A.B, a = a2, b = b2, dfSSE = df.A.B, g = 1, group = "A")
mul2 <- multiplier2[["Tukey"]]
mul3 <- multiplier2[["Scheffe"]]

c2 <- c(1, -1, 0)
CI2<- scary.CI(Salary, MSE.A.B, equal.weights = T, mul2, "A", c2)

# u1. - u3.
c3 <- c(1, 0, -1)
CI3 <- scary.CI(Salary, MSE.A.B, equal.weights = T, mul2, "A", c3)

# u2. - u3.
c4 <- c(0, 1, -1)
CI4 <- scary.CI(Salary, MSE.A.B, equal.weights = T, mul2, "A", c4)

### two contracst CIs
# u1. - (u2.+ u3.)/2
c5 <- c(1, -1/2, -1/2)
CI5 <- scary.CI(Salary, MSE.A.B, equal.weights = T, mul3, "A", c5)

# u21 - (u11 + u31)/2
c6 <- matrix(0, nrow = a2, ncol = b2)
c6[2,1] <- 1
c6[1,1] <- -1/2
c6[3,1] <- -1/2
CI6 <- scary.CI(Salary, MSE.A.B, equal.weights = T, mul3, "AB", c6)

# summary of CI
summary.CI <- data.frame(CI1, CI2, CI3, CI4, CI5, CI6)
summary.CI
rbind(c("u11-u12","u1.-u2.","u1.-u3.","u2.-u3.","u1.-(u2.+u3.)/2","u21-(u11+u31)/2"), summary.CI)







