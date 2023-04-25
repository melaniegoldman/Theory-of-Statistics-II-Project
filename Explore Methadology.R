
location <- file.choose()

twins <- read.csv(location)


head(twins)



###############
## Paper
library("dbarts")
library("remotes")
library("BayesTree")
remotes::install_github("vdorie/dbarts")


xt=as.matrix(sapply(data.frame(usek[,-1]), as.double))
xp=as.matrix(sapply(data.frame(usek[usek$treatment==1,-1]), as.double))
xp[,1]=0
y=as.double(usek[,1])
bart.tot <- bart(x.train=xt,   y.train=y,  x.test=xp)


usek = twins[, c("outcome","treatment", covs)]
twins_data <- usek


xt = as.matrix(sapply(data.frame(usek[, -1]), as.double))
xp = as.matrix(sapply(data.frame(usek[usek$treatment == 1, -1]), as.double))
xp[,1] = 0
y = as.double(usek[,1])
bart.tot <- bart(x.train = xt, y.train = y, x.test = xp)



bart(
  x.train, y.train, x.test = matrix(0.0, 0, 0),
  sigest = NA, sigdf = 3, sigquant = 0.90,
  k = 2.0, power = 2.0, base = 0.95, 
  splitprobs = 1 / numvars,
  binaryOffset = 0.0, weights = NULL,
  ntree = 200, ndpost = 1000, nskip = 100,
  printevery = 100, keepevery = 1, keeptrainfits = TRUE,
  usequants = FALSE, numcut = 100, printcutoffs = 0,
  verbose = TRUE, nchain = 1, nthread = 1, combinechains = TRUE,
  keeptrees = FALSE, keepcall = TRUE, sampleronly = FALSE,
  seed = NA_integer_, proposalprobs = NULL,
  keepsampler = keeptrees
)


###############
## bartcs()
library("bartcs")

data(ihdp, package = "bartcs")

single_bart(
  Y = ihdp$y_factual,
  trt = ihdp$treatment,
  X = ihdp[, 6:30],
  num_tree = 10,
  num_chain = 2,
  num_post_sample = 20,
  num_burn_in = 10,
  verbose = FALSE
)

