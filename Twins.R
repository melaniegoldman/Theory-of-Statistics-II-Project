library(maditr)
twins <- read.csv("TWINS.csv")
#Select rows with babies waying less than 2 kgs
twins <- twins[twins$dbirwt_0 < 2000 & twins$dbirwt_1 < 2000,]
#Randomly select which twin 
chosen_twin <- rbinom(n=nrow(twins), size=1, prob=0.05)
twins$chosen_twin <- chosen_twin
#define treatment as being born heavier
twins$treatment <- ifelse(twins$chosen_twin == 0, twins$dbirwt_0 > twins$dbirwt_1,twins$dbirwt_1 > twins$dbirwt_0)
twins$treatment <- as.integer(twins$treatment)
#define outcome as whether or not selected twin survived first year of life
twins$outcome <- ifelse(twins$chosen_twin == 0,twins$mort_0, twins$mort_1)

covs=c("alcohol", "anemia", "birattnd","birmon", "cardiac", "chyper", "cigar6", "crace", "csex", 
       "dfageq", "diabetes", "dlivord_min", "dmar", "drink5", "dtotord_min", "eclamp", "feduc6", 
       "frace", "gestat10", "hemo", "herpes", "hydra", "incervix", "lung", "mager8", "meduc6",
       "mpre5", "mrace", "nprevistq", "orfath", "ormoth", "othermr" ,"phyper", "pre4000", "preterm", 
       "renal", "rh", "tobacco", "uterine")
ncovs=length(covs)
#remove data with NA values
twins <- na.omit(twins)


# Try with binary response
library(BayesTree)
usek=twins[,c("outcome","treatment",covs)]
twins_data <- usek
write.csv(twins_data, file = "twins_data.csv")
xt=as.matrix(sapply(data.frame(usek[,-1]), as.double))
xp=as.matrix(sapply(data.frame(usek[usek$treatment==1,-1]), as.double))
xp[,1]=0
y=as.double(usek[,1])
bart.tot <- bart(x.train=xt,   y.train=y,  x.test=xp)


