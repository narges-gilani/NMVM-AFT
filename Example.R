
WD.PATH = paste("/Users/gthb3/Desktop/New folder")

source(paste(WD.PATH, '/Gibbs-sampling.r', sep=""))
source(paste(WD.PATH, '/MH_nu.r', sep=""))
source(paste(WD.PATH, '/NMV_CLR_Bayes.r', sep=""))
source(paste(WD.PATH, '/NMVMdist.r', sep=""))

#############################################################

library(survival)
head(gbsg)#
table(is.na(gbsg))
y=log(gbsg$rfstime); cc=gbsg$status
x=cbind(1, gbsg$age, gbsg$meno,(gbsg$size),gbsg$grade, gbsg$nodes, gbsg$hormon)


f1 = NMV.CLR.Bayes(y,
                   x,
                   Family = "GHST",
                   cens = "Right",
                   cc,
                   influence = TRUE,
                   spacing = 1,
                   hyper_set = NA,
                   prior = "Unif",
                   hyper= 1,
                   n.thin = 1,
                   burnin = 2000,
                   n.iter = 22000,
                   n.chains = 1,
                   sample.out = TRUE,
                   table = NULL,
                   level = 95,
                   Prog.Bar = TRUE,
                   Print = TRUE,
                   plot = FALSE)


