######################################################
######       STA 106 / Project1 / Liya Li       ######
######                 Appendix                 ######
######################################################

# The goal of this experiment was to see if different nests for sparrows on Kent Island attracted different size sparrows. 
library(ggplot2)

###### read in the data
sparrow <- read.csv("~/Desktop/Winter Quarter 2018/STA 106/106_project1/sparrow.csv", header = TRUE)
View(sparrow)
head(sparrow)
tail(sparrow)

###### explore the data
total <- dim(sparrow)[1] # nt
total 
str(sparrow)             # variables
a = length(unique(sparrow$Treatment)) 
a                        # a

# number of observations, mean, sd, var in each group
n.i <- aggregate(Weight ~ Treatment, data = sparrow, length)[,2]
mu.i <- aggregate(Weight ~ Treatment, data = sparrow, mean)[,2]
sd.i <- aggregate(Weight ~ Treatment, data = sparrow, sd)[,2]
var.i <- aggregate(Weight ~ Treatment, data = sparrow, var)[,2]

dataframe <- data.frame(groups=c("control","enlarged","reduced"), n.i, mu.i, sd.i, var.i)
dataframe

# plot the data
plot(sparrow$Weight, sparrow$Treatment, xlab = "Weight (gram)", ylab = "Treatment (1.0:Control, 2.0:Enlarged, 3.0:Reduced)", main = "Relationship between Weight and Treatment")

### outliers via boxplot
# boxplot of weight and treatment
boxplot(Weight ~ factor(Treatment), data = sparrow, xlab = "Treatment (C,E,R)", ylab = "Weight (gram) ", main = "Boxplot: Weight vs Treatments")

### outliers via semi-studentized/standardized residuals
# fit the model 
the.model <- lm(Weight ~ Treatment, data = sparrow)
the.model


# add residuals column to the dataset and view the dataset 
sparrow$ei <- the.model$residuals
View(sparrow)

# find SSE, MSE, eij.star
SSE <- sum(sparrow$ei^2)
SSE
MSE <- SSE / (total-a)
MSE  
eij.star <- the.model$residuals / sqrt(MSE)
  
# find the outliers
alpha <- 0.05
t.cutoff <- qt(1-alpha/(2*total), total-a)
t.cutoff
CO.eij <- which(abs(eij.star) > t.cutoff)
CO.eij


### outliers via studentized/standardized residuals, create a histogram of studentized residuals
rij = rstandard(the.model)
CO.rij = which(abs(rij) > t.cutoff)
CO.rij




### the.model
# create diagnostic plot for each of the two main assumptions about eij
plot(the.model$fitted.values, the.model$residuals, main = "Old Model: Errors vs. Fitted Values", xlab = "Fitted Values", ylab = "Errors" )
abline(h = 0,col = "purple")

qqnorm(the.model$residuals, pch = 21, font = 1, font.lab =1, cex =1.5, cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
qqline(the.model$residuals)

# use tests to check normality and equal variance assumptions
the.shapiro.pval = shapiro.test(the.model$residuals)$p.val
the.shapiro.pval
library(car)
the.BFtest = leveneTest(ei ~ Treatment, data=sparrow, center=median)
the.BFtest
the.BF.pval = the.BFtest[[3]][1]
the.BF.pval

# remove the outliers, create new model
outliers <- CO.rij
new.sparrow <- sparrow[-outliers,]
View(new.sparrow)

new.model <- lm(Weight ~ Treatment, data = new.sparrow)
new.model

### new.model
# create diagnostic plot for each of the two main assumptions about eij
plot(new.model$fitted.values, new.model$residuals, main = "New Model: Errors vs. Fitted Values", xlab = "Fitted Values", ylab = "Errors" )
abline(h = 0,col = "blue")

qqnorm(new.model$residuals, pch = 21, font = 1, font.lab =1, cex =1.5, cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
qqline(new.model$residuals)

# use tests to check normality and equal variance assumptions
new.shapiro.pval = shapiro.test(new.model$residuals)$p.val
new.shapiro.pval
library(car)
new.BFtest = leveneTest(ei ~ Treatment, data=new.sparrow, center=median)
new.BFtest
new.BF.pval = new.BFtest[[3]][1]
new.BF.pval


old <- c(the.shapiro.pval,the.BF.pval,alpha)
new <- c(new.shapiro.pval,new.BF.pval,alpha)
rownames <- c("Shapiro.pval","BF.pval","Alpha")
df.test <- data.frame(rownames,old,new)
df.test


# report back the model fit and other computations
# SFA
# create ANOVA table (obtain SS, MS, df{})
anova.table <- anova(the.model)
anova.table
SSE <- anova.table[2,2]
SSA <- anova.table[1,2]
MSA <- anova.table[1,3]
MSE <- anova.table[2,3]
# SSTO = var(sparrow$Weight)*(total -1)
SSTO <- colSums(anova.table)[2]
SSTO
Fs <- anova.table[1,4]
Fs
Pval <- anova.table[1,5]
Pval
# alpha = 0.05

give.me.power = function(ybar,ni,MSE,alpha){
  a = length(ybar) 
  nt = sum(ni) 
  overall.mean = sum(ni*ybar)/nt 
  phi = (1/sqrt(MSE))*sqrt( sum(ni*(ybar - overall.mean)^2)/a) #Finds the books value of phi
  phi.star = a *phi^2 #Finds the value of phi we will use for R 
  Fc = qf(1-alpha,a-1,nt-a) #The critical value of F, use in R's function
  power = 1 - pf(Fc, a-1, nt-a, phi.star)# The power, calculated using a non-central F
  return(power)
}

give.me.power(mu.i,n.i,MSE,alpha=0.05)

give.me.CI = function(ybar,ni,ci,MSE,multiplier){
  if(sum(ci) != 0 & sum(ci !=0 ) != 1){
    return("Error - you did not input a valid contrast")
  } else if(length(ci) != length(ni)){
    return("Error - not enough contrasts given")
  }
  else{
    estimate = sum(ybar*ci)
    SE = sqrt(MSE*sum(ci^2/ni))
    CI = estimate + c(-1,1)*multiplier*SE
    result = c(estimate,CI)
    names(result) = c("Estimate","Lower Bound","Upper Bound")
    return(result)
  }
}


# For this project we use Bonferroni multiplier
# CI for group reduced
c1 <- c(0,0,1)
g1=1
Bon1 = qt(1-alpha/(2*g1),total-a)
mu.Reduced <- give.me.CI(mu.i,n.i,c1,MSE,Bon1)


# CI for group control vs enlarged
c2 <- c(1,-1,0)
g2=2
Bon2 = qt(1-alpha/(2*g2),total-a)
mu.Control.Enlarged <- give.me.CI(mu.i,n.i,c2,MSE,Bon2)

# CI for group control vs reduced
c3 <- c(1,0,-1)
g3=2
Bon3 = qt(1-alpha/(2*g3),total-a)
mu.Control.Reduced <- give.me.CI(mu.i,n.i,c3,MSE,Bon3)

df.CI <- data.frame(mu.Reduced,mu.Control.Enlarged,mu.Control.Reduced)
df.CI












  
  






