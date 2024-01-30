# Set our directory
setwd("")

#Packages we need
library(geiger)
library(phytools)
library(plyr)
library(ggplot2)
library(Hmisc)
library(tidyverse)

library(mvMORPH)
library(RPANDA)
library(OUwie)
library(emdbook)
library(phangorn)
library(ggtext)
library(glue)
library(cowplot)

## For Bio_1
#tree<-read.nexus("AnolisFullyResolvedCloseToConsensus.tre")
tree<-read.nexus("mcc_thinned_allruns.trees") # We're just reading the tree we need here
# # # fixing tree 
#tree$edge.length <- tree$edge.length * 1000
#tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE)
#is.ultrametric(tree)

# dataset
datos <- read.csv("pgls_data_bio1.csv") # Here we need, the species names, and the bio 1 values for each species
regions <- read.csv("pgls_data_bio1_regions.csv") # Here we need, the species names and if they are island or mainland species (for some reason I have to upload to separate data frames)
datos.regions <- merge(datos, regions, "species")
datos <- datos.regions[,c(1,5,11)]
rownames(datos) <- datos$species


# match tree with data
TreeOnly <- setdiff(tree$tip.label,rownames(datos))
TreeOnly # Enter the name of the object we just created to see what's in it.
DataOnly <- setdiff(rownames(datos), tree$tip.label)
DataOnly # Enter to see what species are in the data set but not the tree

# tree pruned
pruned <- drop.tip(tree, tip=TreeOnly)

# select areas
reg <- datos$sregions
names(reg)<-row.names(datos)

# fit data to other models
bio_01<-datos$Bio_01
names(bio_01) <- datos$species
sd<-sd(bio_01)

# fit data to OUWie
datos.ouwie <- as.data.frame(datos)
rownames(datos.ouwie) <- c()
datos.ouwie <- datos.ouwie[,c(1,3,2)]
datos.ouwie["mserr"]<-sd # Here we are calculated the standard deviation for all the species in our dataset


# SIMMAP
simmap.reg <- make.simmap(pruned, reg, model="ER", nsim=500)
summary.simmap <-summary(simmap.reg)


# Model fitting
# Geographical using OUwie
# bm <- OUwie(simmap.reg[[1]], datos.ouwie, model="BM1", simmap.tree=TRUE, mserr = "known")
# bms <- OUwie(simmap.reg[[1]], datos.ouwie, model="BMS", simmap.tree=TRUE, mserr = "known")
# ou.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OU1", simmap.tree=TRUE, mserr = "known")
# oum.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUM", simmap.tree=TRUE, mserr = "known")
# oumv.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMV", simmap.tree=TRUE, mserr = "known")
# ouma.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMA", simmap.tree=TRUE, mserr = "known")
# oumva.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMVA", simmap.tree=TRUE, mserr = "known")
# 
# aic.all <- c(bm$AIC, bms$AIC, ou.ouwie$AIC, oum.ouwie$AIC, oumv.ouwie$AIC, ouma.ouwie$AIC, oumva.ouwie$AIC)
# aic.all


## We chose the model with the lowest AIC value, in this case it was OUMVA

OUMVA.fitted.pruned <- lapply(1:500, function(x) OUwie(simmap.reg[[x]], datos.ouwie, model="OUMVA", simmap.tree=TRUE, mserr = "known")) # We're going to run this analysis 500 times
