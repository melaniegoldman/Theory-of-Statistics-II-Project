####################################################################
####################################################################
## Outline
##
## 1) Packages and Setup
## 2) Data Preparation
## 3) BART Evaluation
## 4) Linear Fitting Evaluation
## 5) Propensity Score Evaluation
## 6) Comparitive Metrics
## 7) Plotting
## 8) Other



####################################################################
## 1) Packages and Setup
library("maditr")
library("tidyverse")
"%!in%" <- function(x,y)!('%in%'(x,y))
set.seed(23)



####################################################################
## 2) Data Preparation
twins <- read.csv("TWINS.csv")
#Select rows with babies weighing less than 2 kgs
twins <- twins[twins$dbirwt_0 < 2500 & twins$dbirwt_1 < 2500,]

#Randomly select which twin 
chosen_twin <- rbinom(n=nrow(twins), size=1, prob=0.5)
twins$chosen_twin <- chosen_twin

#define treatment as being born heavier
twins$treatment <- ifelse(twins$chosen_twin == 0, twins$dbirwt_0 > twins$dbirwt_1,twins$dbirwt_1 > twins$dbirwt_0)
twins$treatment <- as.integer(twins$treatment)
#define outcome as whether or not selected twin survived first year of life
twins$outcome <- ifelse(twins$chosen_twin == 0,twins$mort_0, twins$mort_1)

#select birth order of the chosen twin
twins$bord <- ifelse(twins$chosen_twin == 0, twins$bord_0,twins$bord_1)

n <- nrow(twins)
#Percent of Lighter Twins that Died
(sum(twins$mort_0)/n)*100
# ~8.6%
#Percent of Heavier Twins that Died
(sum(twins$mort_1)/n)*100
# ~7.2%

removed <- c("mort_0","mort_1","dbirwt_0","dbirwt_1", "infant_id_0", 
             "infant_id_1","bord_0", "bord_1")
covs=c("pldel","birattnd","brstate","stoccfipb","mager8","ormoth", "mrace","meduc6","dmar", 
       "mplbir","mpre5", "adequacy", "orfath", "frace", "birmon", "gestat10", "csex", "anemia", 
       "cardiac", "lung", "diabetes", "herpes", "hydra", "hemo", "chyper", "phyper", "eclamp", 
       "incervix" ,"pre4000", "preterm", "renal", "rh", "uterine", "othermr",
       "tobacco", "alcohol", "cigar6", "drink5",   "crace","data_year", "nprevistq", 
       "dfageq", "feduc6",  "dlivord_min", "dtotord_min", "bord",
       "brstate_reg","stoccfipb_reg", "mplbir_reg") 
ncovs=length(covs)
#export data
write.csv(twins, file = "twins_data.csv", row.names = F)





####################################################################
## 3) BART Evaluation
# read in filtered data
twins <- read.csv("twins_data.csv")

#Estimate Casual Effects using BART (Melanie 4/27/23)---------------------------

#remove NA values
twins <- na.omit(twins)

#Create Training and Test Sets
train_sample <- sample(1:nrow(twins),0.8*nrow(twins), replace = F)
data_train <- twins[train_sample,]
data_test <- twins[-train_sample,]

#Create vectors of outcome and treatment and matrix of confs
outcome <- data_train$outcome
treatment <- data_train$treatment
confs <- data_train[,covs]

#fit model using bartc in bartCause package
library(bartCause)
bart_fit <- bartc(response = outcome, treatment = treatment, confounders = confs, keepTrees = TRUE )
CATE <- fitted(bart_fit, type = "cate")
PATE <- fitted(bart_fit, type = "pate")
SATE <- fitted(bart_fit, type = "sate")

#END ---------------------------------------

# Run BART with binary response
library(BART)
usek = twins[,c("outcome","treatment",covs)]
xt = as.matrix(sapply(data.frame(usek[,-1]), as.double))
xp = as.matrix(sapply(data.frame(usek[usek$treatment == 1,-1]), as.double))
xp[,1] = 0
y = usek[,1]


bart.tot <- pbart(x.train=xt,   y.train=y,  x.test=xp)
# first just effect of treatment on treated
diffs=bart.tot$yhat.train[,usek$treatment==1]-bart.tot$yhat.test
mndiffs=apply(diffs,1,mean)
mean(mndiffs)
sd(mndiffs)

# get a sense of t.e. heterogeneity
hist(apply(diffs,2,mean))

usek$outcome <- factor(usek$outcome) #convert to factor
lin_mod <- glm( outcome~., family = "binomial" ,data = usek)



# Run BART with ‘BayesTree’, used in the primary papers package
library("BayesTree")

bart(x.train = xt, y.train = y, binaryOffset = 1)


f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0 #y = f(x) + sigma*z , z~N(0,1)
n = 100 #number of observations
set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y=Ey+sigma*rnorm(n)
lmFit = lm(y~.,data.frame(x,y)) #compare lm fit to BART later
##run BART
set.seed(99)
bartFit = bart(x,y,ndpost=200) #default is ndpost=1000, this is to run example fast.
plot(bartFit) # plot bart fit
##compare BART fit to linear matter and truth = Ey
fitmat = cbind(y,Ey,lmFit$fitted,bartFit$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))



# Other stuff from pbart example
geweke <- gewekediag(bart.tot$yhat.train)
plot(geweke$z, pch='.', cex=2, ylab='z', xlab='i',
     sub=paste0('N:', n, ', k:', k),
     xlim=c(1, n), ylim=c(-5, 5))
lines(1:n, rep(-1.96, n), type='l', col=6)
lines(1:n, rep(+1.96, n), type='l', col=6)
lines(1:n, rep(-2.576, n), type='l', col=5)
lines(1:n, rep(+2.576, n), type='l', col=5)
lines(1:n, rep(-3.291, n), type='l', col=4)
lines(1:n, rep(+3.291, n), type='l', col=4)
lines(1:n, rep(-3.891, n), type='l', col=3)
lines(1:n, rep(+3.891, n), type='l', col=3)
lines(1:n, rep(-4.417, n), type='l', col=2)
lines(1:n, rep(+4.417, n), type='l', col=2)
text(c(1, 1), c(-1.96, 1.96), pos=2, cex=0.6, labels='0.95')
text(c(1, 1), c(-2.576, 2.576), pos=2, cex=0.6, labels='0.99')
text(c(1, 1), c(-3.291, 3.291), pos=2, cex=0.6, labels='0.999')
text(c(1, 1), c(-3.891, 3.891), pos=2, cex=0.6, labels='0.9999')
text(c(1, 1), c(-4.417, 4.417), pos=2, cex=0.6, labels='0.99999')




####################################################################
## 4) Linear Fitting Evaluation
twinsSubset   <- data_train[,c(names(twins)[names(twins) %in% covs], "outcome", "treatment")]

linearFitting <- glm(outcome ~ ., data = twinsSubset)

linearFitting$coefficients["treatment"]



####################################################################
## 5) Propensity Score Evaluation
library("Matching")

# Using the linear fitted model above
qx = linearFitting$fitted

## genetic matching should also explicitly control for qx
twinsProp = cbind.data.frame(twinsSubset, "qx" = qx)
## use this just to get the design matrix needed for GenMatch
modqx2 = glm(outcome ~ ., data = twinsProp)

## need formula for balance statistics to explicitly control what quadratic
## terms are included so we don't waste time with squared binary variables
form.quadz <- as.formula("outcome ~ qx + (pldel + birattnd + brstate + stoccfipb + 
                          mager8 + ormoth + mrace + meduc6 + dmar + mplbir + mpre5 + 
                          adequacy + orfath + frace + birmon + gestat10 + csex + 
                          anemia + cardiac + lung + diabetes + herpes + hydra +
                          hemo + chyper + phyper + eclamp + incervix + pre4000 + 
                          preterm + renal + rh + uterine + othermr + tobacco + 
                          alcohol + cigar6 + drink5 + crace + data_year + nprevistq + 
                          dfageq + feduc6 + dlivord_min + dtotord_min + bord + 
                          brstate_reg + stoccfipb_reg + mplbir_reg)^2")

form.quadz <- as.formula("dose400 ~ qx + (bw + momage + nnhealth + birth.o + parity + moreprem + cigs + alcohol + ppvt.imp + bwg + female + mlt.birtF + b.marryF + livwhoF + languageF + whenprenF + drugs + othstudy + momed4F + siteF + momraceF + workdur.imp)^2")
mod.quad=glm(formula=form.quadz,data=usek2,x=TRUE)
ncol(mod.quad$x)


####################################################################
## 6) Comparitive Metrics

####################################################################
## 7) Plotting

####################################################################
## 8) Other
