#------------------------------------------------------------------------------------------
#Generalized RDA approach

mm<-read.csv("Target_genetic_data_with_climate_and_distance.csv",header=TRUE)

genetic<- mm[, X:X] #columns with genetic data (insert numbers where there is a X)
climate<- mm[, X:X] #columns with climate variables
corrcheck<-cor(climate) #check variables to remove highly correlated ones

install.packages("vegan") #package needed for RDA
library(vegan) #calling package for RDA
RDA <- rda(genetic ~ chosen_climate_variables , climate) #model with just climate, use a + between chosen variables
anova(RDA)
plot(RDA)
#making conditional on the effect of distance of neutral genetic structure
pRDA1<-rda(genetic ~ chosen_climate_variables+Condition(distance_variables+variables_accounting_for_genetic_structure), climate) #accounts for distance will look at outliers from this output
anova(pRDA1)
plot(pRDA1)
#partioning variance, looking at distance constraining other sources of variation
pRDA2<-rda(genetic ~ distance_variables + Condition(chosen_climate_variables+variables_accounting_for_genetic_structure), climate)
anova(pRDA2)
plot(pRDA2)
#partioning variance, looking at neutral genetic structure constraining other sources of variation
pRDA3<-rda(genetic ~ variables_accounting_for_genetic_structure + Condition(chosen_climate_variables+distance_variables), climate)
anova(pRDA3)
plot(pRDA3)
#Full model with all factors as variables
RDAfull<-rda(genetic ~ distance_variables +variables_accounting_for_genetic_structure+chosen_climate_variables+distance_variables, climate)
anova(RDAfull)
plot(RDAfull)
#Look at inertia values from these to understand how variance is partioned and the relative importance of climate, distance, and netural genetic structure
(RDAfull)
(pRDA1)
(pRDA2)
(pRDA3)
anova(pRDA1, by="axis", permutations = 999)#how many signifcant axes, there are two for constrained

#Finding outlier SNPs, in this instance just looking at the model using climate assuming two significant axes
load.rda <- scores(pRDA1, choices=c(1:2), display="species")  # Species scores for the first two constrained axes
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],2.5) #outliers for axis 1
cand2 <- outliers(load.rda[,2],2.5) #outliers for axis 2


ncand <- length(cand1) + length(cand2)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))

colnames(cand1) <- colnames(cand2) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=X)  #assuming X important climate variables 
colnames(foo) <- c("X") #insert the names of chosen climates variables using a comma between entries

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- genetic[,nam]
  foo[i,] <- apply(annual,2,function(x) cor(x,snp.gen))
}
cand <- cbind.data.frame(cand,foo)  
check<-head(cand)

#Check for duplicates
length(cand$snp[duplicated(cand$snp)]) #1 duplicate
foo <- cbind(cand$axis, duplicated(cand$snp)) #axis 1 check
table(foo[foo[,1]==1,2])
table(foo[foo[,1]==2,2]) #axis 2 check
cand <- cand[!duplicated(cand$snp),] #remove the duplicate


#which of the predictors each candidate SNP is most strongly correlated with:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,X] <- names(which.max(abs(bar[X:X]))) # gives the variable, replace instances of X with appropiate column numbers 
  cand[i,X] <- max(abs(bar[X:X]))              # gives the correlation
}

colnames(cand)[X] <- "predictor"
colnames(cand)[X] <- "correlation"

table(cand$predictor) 

