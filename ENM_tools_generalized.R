install.packages("devtools")
library(devtools)

install.packages("ENMTools", type = "source")
library(ENMTools)

install.extras()
library(raster)
library(sp)
climate_var_1 <- raster("climate_file.tif")

env<-stack(X)#all climate variables of interest seperate by a comma
env <- crop(env, extent(X, X, X, X))# crop stack of climate variables to area of interest
plot(env)

target_species <- enmtools.species(species.name = "target_species", 
                              presence.points = read.csv("target_presence_points.csv"))
target_species <- check.species(target_species)
interactive.plot.enmtools.species(target_species)#check presence data
check<-raster.cor.matrix(env)#check co-lin of predictors
raster.cor.plot(env)#now plotted


target_species.glm <- enmtools.glm(species = target_species, env = env,nback=X, test.prop = X)#running a basic GLM approach, change Xs to desired parameter value

# to test to see if the niche is identical or not for two different target groups

target_one <- enmtools.species(species.name = "target_one", 
                                presence.points = read.csv("pres_points_one.csv"))
target_one <- check.species(target_one)

target_two <- enmtools.species(species.name = "target_two", 
                               presence.points = read.csv("8pres_points_two.csv"))
target_two <- check.species(target_two)

id.glm <- identity.test(species.1 = target_one, species.2 = target_two, env = env, type = "glm", nreps = X,nback=X) #idenity test, change Xs to desired parameter value