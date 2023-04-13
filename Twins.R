library(maditr)
twins <- read.csv("TWINS.csv")
#Select rows with babies waying less than 2 kgs
twins <- twins[twins$dbirwt_0 < 2500 & twins$dbirwt_1 < 2500,]

#Randomly select which twin 
chosen_twin <- rbinom(n=nrow(twins), size=1, prob=0.05)
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

#read in filtered data
twins <- read.csv("twins_data.csv")


# Run BART with binary response
library(BayesTree)
library(BART)
usek=twins[,c("outcome","treatment",covs)]
xt=as.matrix(sapply(data.frame(usek[,-1]), as.double))
xp=as.matrix(sapply(data.frame(usek[usek$treatment==1,-1]), as.double))
xp[,1]=0
y=usek[,1]


bart.tot <- pbart(x.train=xt,   y.train=y,  x.test=xp)

#pbart function
library(BART)
B <- getOption('mc.cores', 1)
figures = getOption('figures', default='NONE')

data(arq)
str(arq)
arth <- as.matrix(arq)

twins_data <- as.matrix(twins_data)
N <- length(twins_data[ , 'gestat10'])

gestat10 <- 0:9
H <- length(gestat10)

post1 <- mc.pbart(x.train=xt, y.train=y,
                  mc.cores=B, seed=99)

for(i in 1:2)
  for(j in 1:H) {
    x. <- xt
    x.[ , 'treatment'] <- i
    x.[ , 'gestat10'] <- gestat10[j]
    if(i==1 && j==1) x.test <- x.
    else x.test <- rbind(x.test, x.)
  }

pred1 <- predict(post1, newdata=x.test, mc.cores=B)

M <- nrow(pred1$prob.test)
##Friedman's partial dependence function
pd1 <- matrix(nrow=M, ncol=H)
pd2 <- matrix(nrow=M, ncol=H)

par(mfrow=c(1, 2))

for(j in 1:H) {
  h <- (j-1)*N
  pd1[ , j] <- apply(pred1$prob.test[ , h+1:N], 1, mean)
  h <- h+N*H
  pd2[ , j] <- apply(pred1$prob.test[ , h+1:N], 1, mean)
}

pd1.mean <- apply(pd1, 2, mean)
pd2.mean <- apply(pd2, 2, mean)
pd1.025 <- apply(pd1, 2, quantile, probs=0.025)
pd2.025 <- apply(pd2, 2, quantile, probs=0.025)
pd1.975 <- apply(pd1, 2, quantile, probs=0.975)
pd2.975 <- apply(pd2, 2, quantile, probs=0.975)

plot(gestat10, pd1.mean, type='l', col='blue',
     ylim=0:1, xlab='No. of gestation weeks prior to birth.', ylab=expression(p(x)),
     sub='Treatment: 1(blue) vs. 0(red)')
##sub='Low-back/buttock pain: M(blue) vs. F(red)')
lines(gestat10, pd1.025, type='l', col='blue', lty=2)
lines(gestat10, pd1.975, type='l', col='blue', lty=2)
lines(gestat10, pd2.mean, type='l', col='red')
lines(gestat10, pd2.025, type='l', col='red', lty=2)
lines(gestat10, pd2.975, type='l', col='red', lty=2)
lines(gestat10, rep(0, H), type='l')
lines(gestat10, rep(1, H), type='l')

