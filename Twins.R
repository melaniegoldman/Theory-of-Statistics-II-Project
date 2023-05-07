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
library("bartCause")
library("BART")
library("Matching")
library("MatchIt")
library("rbounds")
library("rgenoud")
library("broom")
library("plotly")
library("pROC")
library("remotes")
library("plotBart")
  # remotes::install_github("priism-center/plotBart")

"%!in%" <- function(x,y)!('%in%'(x,y))

set.seed(23)



####################################################################
## 2) Data Preparation
twins <- read.csv("TWINS.csv")
#Select rows with babies weighing less than 2 kgs
twins <- twins[twins$dbirwt_0 < 2500 & twins$dbirwt_1 < 2500, ]

#Randomly select which twin 
chosen_twin <- rbinom(n = nrow(twins), size = 1, prob = 0.5)
twins$chosen_twin <- chosen_twin

#define treatment as being born heavier
twins$treatment <- ifelse(twins$chosen_twin == 0, twins$dbirwt_0 > twins$dbirwt_1, twins$dbirwt_1 > twins$dbirwt_0)
twins$treatment <- as.integer(twins$treatment)
#define outcome as whether or not selected twin survived first year of life
twins$outcome <- ifelse(twins$chosen_twin == 0, twins$mort_0, twins$mort_1)

#select birth order of the chosen twin
twins$bord <- ifelse(twins$chosen_twin == 0, twins$bord_0, twins$bord_1)

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



############ alternate
twins <- read.csv("TWINS.csv")
names(twins)[1] <- "pairID"

# Find the twin pairs with a 20% or more weight difference and define
#     treatment = 1 for heavier twin and = 0 for the lighter twin
# for(i in 1:nrow(twins)){
#   weights <- twins[i, c("dbirwt_0", "dbirwt_1")]
#   twins$weightDif[i] <- ifelse(max(weights)/min(weights) >= 1.2, 1, 0)
# }

weight_pct_diff <- abs(twins$dbirwt_1 - twins$dbirwt_0)/((twins$dbirwt_1 + twins$dbirwt_0)/2)
twins$weightDif <- ifelse(weight_pct_diff > 0.2, 1, 0)


# Subset base on the twins weight difference
twins <- twins[twins$weightDif == 1, names(twins) %!in% "weightDif"]



# Rows = 1 have one twin < 2.5kg, = 2 have both twins < 2.5kg, 0 = none < 2.5kg
for(i in 1:nrow(twins)){
  weights <- twins[i, c("dbirwt_0", "dbirwt_1")]
  twins$wght2500[i] <- ifelse(max(weights) >= 2500 & min(weights) < 2500, 1,
                              ifelse(max(weights) < 2500 & min(weights) < 2500, 2, 0))
}

# Subset base on the twins weight minimum
#twinsMed <- twins[twins$wght2500 == 1, names(twins) %!in% "wght2500"]
twinsMed <- twins[twins$wght2500 == 2, names(twins) %!in% "wght2500"]



# Elongate the dataset so each twin has its own row
twin0 <- names(twinsMed)[names(twinsMed) %!in% covs] %>%
         .[str_detect(., "0")] %>% twinsMed[, .] %>% 
         `colnames<-`(c("infant_id", "bord", "outcome", "dbirwt"))

twin1 <- names(twinsMed)[names(twinsMed) %!in% covs] %>%
         .[str_detect(., "1")] %>% twinsMed[, .] %>% 
         `colnames<-`(c("infant_id", "bord", "outcome", "dbirwt"))

elongated <- rbind(cbind(twinsMed$pairID, twin0, twinsMed[, names(twinsMed) %in% covs]),
                   cbind(twinsMed$pairID, twin1, twinsMed[, names(twinsMed) %in% covs]) )
names(elongated)[1] <- "pairID"

for(i in 1:length(unique(elongated$pairID)) ){
  ID <- unique(elongated$pairID)[i]
  
  one <- rbinom(1, 1, 0.5)
  two <- ifelse(one == 1, 0, 1)
  
  elongated[elongated$pairID == ID,"chosen_twin"]  <- c(one, two)
}


# Define the heavier twin as treatment = 1 and the lighter as treatment = 0
for(i in 1:length(unique(elongated$pairID)) ){ 
  subset <- elongated[elongated$pairID == unique(elongated$pairID)[i], ]
  elongated[rownames(subset), "treatment"] <- ifelse(subset$dbirwt == max(subset$dbirwt), 1, 0) 
}

# for(i in 1:length(unique(elongated$pairID)) ){
#   elongated[elongated$pairID == unique(elongated$pairID)[i], "treatment"] <- 
#     ifelse(subset$dbirwt == max(subset$dbirwt), 1, 0)
# }


#write.csv(elongated, file = "twins_One2500.csv", row.names = F)
write.csv(elongated, file = "twins_Both2500.csv", row.names = F)





####################################################################
## 3) BART Evaluation

# read in filtered data
#twinsAll <- read.csv("twins_One2500.csv")
twinsAll <- read.csv("twins_Both2500.csv")

# select out the randomly chosen twin from the whole dataset to
#     simulate observational study settings
twins <- twinsAll[twinsAll$chosen_twin == 1, ] %>% 
            .[, names(.) %!in% "chosen_twin"] %>% `rownames<-`(NULL)

# remove NA values
twins <- na.omit(twins)



#Estimate Casual Effects using BART (Melanie 4/27/23)---------------------------

#Create Training and Test Sets
train_sample <- sample(1:nrow(twins),0.8*nrow(twins), replace = F)
data_train   <- twins[train_sample,]
data_test    <- twins[-train_sample,]

#Create vectors of outcome and treatment and matrix of confs
outcome   <- data_train$outcome
treatment <- data_train$treatment
confs     <- data_train[, covs]

#fit model using bartc in bartCause package
bart_fit  <- bartc(response = outcome, treatment = treatment, 
                   confounders = confs, keepTrees = TRUE,
                   n.burn = 100, estimand  = "ate")

#######################
# Diagnostic plots

# MCMC
plot_trace(bart_fit, type = c("cate", "sate", "pate", "sigma")) +
  theme_minimal() +
  labs(title = "MCMC Trace Plot - Confirm Convergence")
  # chains should be well mixed to show that there is MCMC convergence

# 
plot_balance(data_train, "treatment", names(confs)) +
  xlim(-0.5, 0.5)

# Plot common support
plot_common_support(bart_fit) +
  geom_point(alpha = 0.4, color = "steelblue",
             position = position_jitter(w = 0, h = 0.15)) +
  scale_shape_manual(values = c(3, 20)) +
  labs(title = "Common Support Checks", y = "Outcome")


# Plots showing the three treatment effect metrics
plot_CATE(.model = bart_fit, type = 'density', ci_80 = TRUE, ci_95 = TRUE, .mean = TRUE) + 
  labs(subtitle = 'My comments on the results') +
  theme_minimal()


# LOESS prediction of treatment heterogeneity
plot_moderator_c_loess(bart_fit, as.numeric(data_train$gestat10), line_color = "blue") +
  theme_minimal() +
  labs(title = "Treatment Effect Heterogeneity", 
       subtitle = "LOESS prediction of ICATEs by a continuous covariate",
       x = "Individual CATE", y = "Gestation 10")

# Branching heterogeneity effects
plot_moderator_search(bart_fit, max_depth = 2) +
  theme_minimal() +
  labs(title = "Exploratory Heterogeneous Effects")

# Waterfall
plot_waterfall(bart_fit, descending = TRUE, .order = NULL, .color = NULL, 
               .alpha = 0.5) + theme_minimal() + 
  geom_point(color = "blue", size = 0.5) + 
  labs(title = "Waterfall Plot", x = "Ordered ICATE's", y = "Individual CATE",
       subtitle = "Point and posterior intervals of each individual’s ICATE",
       caption = "Steeper curve implies that the treatment effect varies across individuals in the sample")



####################################################################
## 4) Linear Fitting Evaluation
twinsSubset   <- data_train[, c(names(twins)[names(twins) %in% covs], 
                               "outcome", "treatment")]

linearFitting <- glm(outcome ~ ., family = "binomial", data = twinsSubset)

linearFitting$coefficients["treatment"]


####################################################################
## 5) Propensity Score Evaluation
## source: https://www.youtube.com/watch?v=rHVGj1F1D_4
## source: https://cran.r-project.org/web/packages/MatchIt/vignettes/estimating-effects.html
attach(twinsSubset)

# Defining variables (Tr is treatment, Y is outcome, X are independent variables)
Tr <- cbind(treatment)
Y  <- cbind(outcome)
X  <- cbind(pldel, birattnd, brstate, stoccfipb, mager8, ormoth, mrace, meduc6, dmar, 
            mplbir, mpre5, adequacy, orfath, frace, birmon, gestat10, csex, anemia, 
            cardiac, lung, diabetes, herpes, hydra, hemo, chyper, phyper, eclamp, 
            incervix, pre4000, preterm, renal, rh, uterine, othermr, tobacco, 
            alcohol, cigar6, drink5, crace, data_year, nprevistq, dfageq, feduc6,
            dlivord_min, dtotord_min, bord, brstate_reg, stoccfipb_reg, mplbir_reg)


# Genetic matching
gen1  <- GenMatch(Tr = Tr, X = X, BalanceMatrix = X, pop.size = 100)
detach(twinsSubset)

# Weight the fitting using the genetic matching weight matrix
f <- as.formula(paste("outcome ~", 
                      paste(names(twinsSubset)[!names(twinsSubset) %in% c("outcome")], 
                            collapse = " + ")))

mgen1 <- matchit(f, data = twinsSubset, method = "full", 
                 link = "probit", estimand = "ate",
                 weights = gen1$Weight.matrix, distance = "mahalanobis")

matched_data1 <- match.data(mgen1)

mgen1_fit <- glm(outcome ~ treatment, data = matched_data1)


# Absolute standardized mean difference between before and after matching
plot(summary(mgen1))


# qqPlots to show before and after fitting alignments for all covariates
plot(mgen1)


# Plot to show the skew of the weights
tempDataset <- twinsSubset

glmMod <- linearFitting %>% augment(type.predict = "response", data = tempDataset) %>%
  mutate(wts = 1/ifelse(treatment == 0, 1 - .fitted, .fitted))

glmMod$method <- "glm"


genmatchProp <- mgen1_fit %>% augment(type.predict = "response", data = tempDataset) %>%
  mutate(wts = 1/ifelse(treatment == 0, 1 - .fitted, .fitted))

genmatchProp$method <- "Propensity"


matching <- rbind(glmMod, genmatchProp)


ggplot(matching, aes(x = wts, group = method, fill = method)) + 
  geom_density(alpha = 0.3) + theme_minimal() + xlim(0, 20) +
  geom_vline(xintercept = median(glmMod$wts), size = 1, color="blue") +
  geom_vline(xintercept = median(genmatchProp$wts), size = 1, color="red") +
  scale_fill_discrete(name ="Modeling Method") +
  labs(title = "Propensity Score Weights Skew Plot", x ="Weights", y = "Density",
       caption = "The vertical lines show the median of both densities. Notice how the 
                 genetic matching did not improved the median placement of the weights 
                 across all covariates compared to the standard glm model")



####################################################################
## 6) Comparative Metrics

# Paper metrics for BART
CATE <- fitted(bart_fit, type = "cate")

#######################
# Predictions

# for BART
# NOTE: type must be in 'mu', 'y', 'mu.0', 'mu.1', 'y.0', 'y.1', 
#      'icate', 'ite', 'p.score'
# NOTE: or standard glm(), when type = "response" then the P(Y = 1|X) is the output
# source: http://www.science.smith.edu/~jcrouser/SDS293/labs/lab4-r.html#:~:text=The%20predict()%20function%20can,information%20such%20as%20the%20logit%20.
testData <- data_test[, c(covs, "treatment")]

names(testData)[names(testData) == "treatment"] <- "z"


data_test$predBARTy0 <- predict(bart_fit, testData, type = "y.0") %>% t() %>%
                        as.data.frame() %>% .[,50]


data_test$predBARTy1 <- predict(bart_fit, testData, type = "y.1") %>% t() %>%
                        as.data.frame() %>% .[,50]


# for glm()
names(testData)[names(testData) == "z"] <- "treatment"

data_test$glmOutcomes <- ifelse(predict.glm(linearFitting, testData, 
                                            type = "response") < 0.5, 0, 1)


# Prediction for standard propensity score matching
data_test$propOutcomes <- ifelse(predict.glm(mgen1_fit, testData, 
                                             type = "response") < 0.5, 0, 1)


####################################################################
## 7) Plotting

## For Melanie:

## Correlation plots with 
##    - treatment ~ outcome
##    - outcome ~ weight

## Histogram overlays - to justify connection with premies etc.?
## counter to the economics paper?
##    - all > 2.5 and < 2.5

#Twins with all weights
twins <- read.csv("TWINS.csv")
chosen_twin <- rbinom(n = nrow(twins), size = 1, prob = 0.5)
twins$treatment <- ifelse(chosen_twin == 0, twins$dbirwt_0 > twins$dbirwt_1, twins$dbirwt_1 > twins$dbirwt_0)
twins$treatment <- as.integer(twins$treatment)
twins$outcome <- ifelse(chosen_twin == 0, twins$mort_0, twins$mort_1)
twins$bord <- ifelse(chosen_twin == 0, twins$bord_0, twins$bord_1)
twins$dbirwt <- ifelse(chosen_twin == 0, twins$dbirwt_0, twins$dbirwt_1)
twins_allWeights <- twins

#Twins with weight both less than 2500 grams
twins <- read.csv("twins_Both2500.csv")
twins_both2500 <- twins[twins$chosen_twin == 1, ] %>% .[, names(.) %!in% "chosen_twin"]
twins_both2500 <- na.omit(twins_both2500)


#plot margin
m <- list(
  l = 50,
  r = 50,
  b = 100,
  t = 100,
  pad = 4
)

#Create histograms for twins of all weights
filtered_twins <- twins_allWeights
hists_allWeights <- lapply(seq_along(filtered_twins), function(x){
  plot_ly(data = filtered_twins, alpha = 0.6,type = "histogram", x = filtered_twins[,x], split =~outcome) %>%
    layout(
      barmode="stack",
      bargap=0.1,
      margin = m,
      legend=list(title=list(text='<b> Outcome </b>')),
      title = paste('Histogram of', names(filtered_twins[x])), xaxis = list(title =names(filtered_twins[x]) ), yaxis = list(title = "Frequency"))
  #hist(filtered_twins[,x], main = paste("Histogram of",names(filtered_twins[x])), xlab = names(filtered_twins[x]))
})
hists_allWeights[[17]] #gestat10

#Create histograms for twins with weights less than 2500 grams
filtered_twins <- twins_both2500
hists_both2500 <- lapply(seq_along(filtered_twins), function(x){
  plot_ly(data = filtered_twins, alpha = 0.6,type = "histogram", x = filtered_twins[,x], split =~outcome) %>%
    layout(
      barmode="stack",
      bargap=0.1,
      margin = m,
      legend=list(title=list(text='<b> Outcome </b>')),
      title = paste('Histogram of', names(filtered_twins[x])), xaxis = list(title =names(filtered_twins[x]) ), yaxis = list(title = "Frequency"))
  #hist(filtered_twins[,x], main = paste("Histogram of",names(filtered_twins[x])), xlab = names(filtered_twins[x]))
})
hists_both2500[[21]]

#Export Histograms - uncomment first to export new plots
# for (i in 1:length(hists_allWeights)){
#   suppressWarnings(
#     export(hists_allWeights[[i]], file = paste0("Histograms (All Twins)/",names(twins_allWeights)[i],"_hist.png"))
#   )
# }
# 
# 
# for (i in 1:length(hists_both2500)){
#   suppressWarnings(
#     export(hists_both2500[[i]], file = paste0("Histograms (Both Less Than 2500g)/",names(twins_both2500)[i],"_hist.png"))
#   )
# }



## For Shelby:

## Difference plot of the prediction on test and expectation

twinsAll

data_test %>% names()


## ROC
# source: https://www.digitalocean.com/community/tutorials/plot-roc-curve-r-programming
roc_bart = data_test %>% roc(outcome, predBARTy1)
roc_glm  = data_test %>% roc(outcome, glmOutcomes)
roc_prop = data_test %>% roc(outcome, propOutcomes)

auc_bart = round(auc(data_test[, 4], data_test$predBARTy1), 4)
auc_glm  = round(auc(data_test[, 4], data_test$glmOutcomes), 4)
auc_prop = round(auc(data_test[, 4], data_test$propOutcomes), 4)


ggroc(list(glm = roc_glm, BART = roc_bart, Propensity = roc_prop), size = 1) +
  ggtitle('ROC Curve') +
  labs(caption = paste0('glm AUC = ', auc_glm, ' and ', 
                       'BART AUC = ', auc_bart, ' and ', 
                       'Propensity AUC = ', auc_prop),
       x = "Specificity", y = "Sensitivity") +
  guides(color = guide_legend(title = "Modeling Method")) +
  theme_minimal()



####################################################################
## 8) Other
# Run BART with binary response
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



# Run BART with ‘BayesTree’, used in the primary papers package
library("BayesTree")

x = as.matrix(data_train[, c(covs, "treatment")]) %>% head()
y = data_train$outcome

BayesTree::bart(x, y, ndpost=200)


f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0 #y = f(x) + sigma*z , z~N(0,1)
n = 100 #number of observations
x = matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y = Ey+sigma*rnorm(n)
sampleData = data.frame(x,y)
lmFit = lm(y ~ ., data = sampleData) #compare lm fit to BART later
##run BART
bartFit = BayesTree::bart(x, y, ndpost=200) #default is ndpost=1000, this is to run example fast.
plot(bartFit) # plot bart fit
##compare BART fit to linear matter and truth = Ey
fitmat = cbind(y,Ey,lmFit$fitted,bartFit$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))


#######################
# Previous matching method:

# Propensity score model 
glm1 <- glm(Tr ~ X, family = binomial(link = "probit"), data = twinsSubset)
summary(glm1)

# Average treatment on the treated effect
rr1 <- Match(Y = Y, Tr = Tr, X = glm1$fitted)
summary(rr1)

# Checking the balancing property
var <- brstate
var <- frace
var <- cigar6
var <- drink5 # MatchBalance did not improve, but qqPlot shows acceptable matching
var <- feduc6
var <- dlivord_min # interesting qqPlot; shows poor matching at higher dlibord_min

mBal_beforeGen <- MatchBalance(Tr ~ X, data = twinsSubset, match.out = rr1, nboots = 0)
qqplot(var[rr1$index.control], var[rr1$index.treated])
abline(coef = c(0, 1), col = 2)

# after gen match
mgen1 <- Match(Y = Y, Tr = Tr, X = X, Weight.matrix = gen1, estimand = "ATE")
summary.Match(mgen1)

mBal <- MatchBalance(Tr ~ X, data = twinsSubset, match.out = mgen1, nboots = 0)

mgen1$weights %>% unique()


# Plot mean treatment and control effect before/after gen matching
a1 <- lapply(mBal_beforeGen$AfterMatching, function(x) x$mean.Tr) %>% unlist()
b1 <- lapply(mBal_beforeGen$BeforeMatching, function(x) x$mean.Tr) %>% unlist()

c1 <- lapply(mBal_beforeGen$AfterMatching, function(x) x$mean.Co) %>% unlist()
d1 <- lapply(mBal_beforeGen$BeforeMatching, function(x) x$mean.Co) %>% unlist()

diff1 <- data.frame("Difference" = b1 - a1, "Eval" = "beforeGen.Tr")
diff2 <- data.frame("Difference" = d1 - c1, "Eval" = "beforeGen.Co")


a2 <- lapply(mBal$AfterMatching, function(x) x$mean.Tr) %>% unlist()
b2 <- lapply(mBal$BeforeMatching, function(x) x$mean.Tr) %>% unlist()

c2 <- lapply(mBal$AfterMatching, function(x) x$mean.Co) %>% unlist()
d2 <- lapply(mBal$BeforeMatching, function(x) x$mean.Co) %>% unlist()

diff3 <- data.frame("Difference" = b2 - a2, "Eval" = "afterGen.Tr")
diff4 <- data.frame("Difference" = d2 - c2, "Eval" = "afterGen.Co")

diff <- rbind(diff1, diff2, diff3, diff4)


# Before/after mean treatment effect
ggplot(diff, aes(Difference)) +
  geom_density(data = subset(diff, Eval == 'beforeGen.Tr'), 
               fill = "red", alpha = 0.2) + 
  geom_density(data = subset(diff, Eval == 'afterGen.Tr'), 
               fill = "blue", alpha = 0.2) + 
  theme_minimal() +
  labs(title = "Xi Treatment Effects", 
       x ="Before - After Matching Mean Treatment Effect", y = "Density",
       caption = "The red curve reflect matching results using standard glm() matching technique and the blue reflects 
       matching results using Genetic Matching with default settings")

ggplot(diff, aes(Difference)) + 
  geom_density(data = subset(diff, Eval == 'beforeGen.Co'), 
               fill = "red", alpha = 0.2) + 
  geom_density(data = subset(diff, Eval == 'afterGen.Co'), 
               fill = "blue", alpha = 0.2) + 
  xlim(-0.25, 1) +
  theme_minimal() +
  labs(title = "Xi Control Effects", 
       x ="Before - After Matching Mean Control Effect", y = "Density",
       caption = "The red curve reflect matching results using standard glm() matching technique and the blue reflects 
       matching results using Genetic Matching with default settings")

