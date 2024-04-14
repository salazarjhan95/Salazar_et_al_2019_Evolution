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

## For Bio_6
#tree<-read.nexus("AnolisFullyResolvedCloseToConsensus.tre")
tree<-read.nexus("mcc_thinned_allruns.trees")
# # # fixing tree 
#tree$edge.length <- tree$edge.length * 1000
#tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE)
#is.ultrametric(tree)


# dataset
datos <- read.csv("pgls_data_bio6.csv")
regions <- read.csv("pgls_data_bio6_regions.csv")
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
reg <- datos$second_devision1.x
names(reg)<-row.names(datos)

# fit data to other models
bio_06<-datos$Bio_06
names(bio_06) <- datos$species
sd<-sd(bio_06)

# fit data to OUWie
datos.ouwie <- as.data.frame(datos)
rownames(datos.ouwie) <- c()
datos.ouwie <- datos.ouwie[,c(1,3,2)]
datos.ouwie["mserr"]<-sd


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
oumva.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMVA", simmap.tree=TRUE, mserr = "known")
oumva.ouwie
# 
# aic.all <- c(bm$AIC, bms$AIC, ou.ouwie$AIC, oum.ouwie$AIC, oumv.ouwie$AIC, ouma.ouwie$AIC, oumva.ouwie$AIC)
# aic.all

OUMVA.fitted.pruned <- lapply(1:500, function(x) OUwie(simmap.reg[[x]], datos.ouwie, model="OUMVA", simmap.tree=TRUE, mserr = "known"))

#For highland
OUMVA.fitted.pruned[[1]]$solution[,1]
sigma2.oumva.mainland <- lapply(1:500, function(x) OUMVA.fitted.pruned[[x]]$solution[,1])
df_main <- data.frame(t(sapply(sigma2.oumva.mainland,c)))
df_main

landmass <- rep(c("highland"), times = 500)
landmass <- data.frame(landmass)
landmass

df_main2 <- cbind(df_main, landmass)
df_main2
#write.csv(df_main, file = "data_mainland_jul20-2023.csv")


#For OUM.fitted.pruned
OUMVA.fitted.pruned[[1]]$solution[,2]
sigma2.oumva.island <- lapply(1:500, function(x) OUMVA.fitted.pruned[[x]]$solution[,2])
df_island <- data.frame(t(sapply(sigma2.oumva.island,c)))
df_island

landmass <- rep(c("lowland"), times = 500)
landmass <- data.frame(landmass)
landmass

df_island2 <- cbind(df_island, landmass)
df_island2

df_bio6_alphasigma <- rbind(df_main2, df_island2)
colnames(df_bio6_alphasigma) <- c("sigma.sq", "alpha", "elevation")
df_bio6_alphasigma

write.csv(df_bio6_alphasigma, file = "df_bio6_alphasigma_bio6_Cal_2.csv")

#df_bio6_alphasigma <- read.csv(file = "df_bio6_alphasigma_700m_Cal.csv")



#For highland
OUMVA.fitted.pruned[[1]]$theta[1,]
theta.oumva.mainland <- lapply(1:500, function(x) OUMVA.fitted.pruned[[x]]$theta[1,])
df_main <- data.frame(t(sapply(theta.oumva.mainland,c)))
df_main

landmass <- rep(c("lowland"), times =500)
landmass <- data.frame(landmass)
landmass

df_main2 <- cbind(df_main, landmass)
df_main2

colnames(df_main2) <- c("theta", "se", "elevation")
df_main2

#write.csv(df_main, file = "data_mainland.csv")


#For lowland
OUMVA.fitted.pruned[[1]]$theta[2,]
theta.oumva.island <- lapply(1:100, function(x) OUMVA.fitted.pruned[[x]]$theta[2,])
df_island <- data.frame(t(sapply(theta.oumva.island,c)))
df_island

landmass <- rep(c("highland"), times = 500)
landmass <- data.frame(landmass)
landmass

df_island2 <- cbind(df_island, landmass)
df_island2

colnames(df_island2) <- c("theta", "se", "elevation")
df_island2

df_bio6_theta <- rbind(df_main2, df_island2)
colnames(df_bio6_theta) <- c("theta", "se", "elevation")
df_bio6_theta

write.csv(df_bio6_theta, file = "df_bio6_theta_NewJan_2.csv")




## Box plot

df_bio6_theta <- read.csv("df_bio6_theta_NewJan_2.csv")

attach(df_bio6_theta)
head(df_bio6_theta)

theta_bio6 <- ggplot(df_bio6_theta, aes(x=elevation, y=theta, color = elevation, alpha = 0.8)) +
  geom_boxplot(aes(color = elevation), size = .75, width=0.8) +
  #geom_boxplot(aes(color = elevation))+
  geom_jitter(alpha = 0.5,  size = 5, position=position_jitter(0.1)) +
  scale_y_continuous(breaks=seq(0, 44, 4), limit = c(0, 44)) +
  scale_color_manual(values=c("#1E52F6","#D2D2D2"),
                     name  = "elevation",
                     
                     breaks = c("highland","lowland"),
                     
                     labels = c("Highland","Lowland")) +
  #4E84C4
  scale_fill_manual(values=c("#1E52F6","#D2D2D2"),
                    name  = "elevation",
                    
                    breaks = c("highland","lowland"),
                    
                    labels = c("Highland","Lowland")) +
  
  
  
  # theme_bw(base_size = 12, base_family = "")+
  theme_classic()+
  
  theme(legend.position = "none")+
  
  
  #theme(panel.border = element_rect(colour = "black", size = 0.5))+
  
  theme(panel.grid.major = element_line(colour =NA),
        panel.grid.minor =element_line(colour = NA))+
  
  
  theme(axis.title.y= element_text( color = "black")) +
  
  labs(x="", y = "Temperature (°C)",
       title = expression(Optimal~Trait~Value~(theta))) +
  
  scale_x_discrete(limits=c("highland","lowland"),  
                   
                   labels= c("Highland","Lowland"))+
  
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        plot.title    = element_text(family = "serif", size = 15)) +
  
  theme(axis.text.x = element_text(family="serif", size = 13, lineheight = 0.9, vjust = 1))+
  theme(axis.text.y = element_text(family="serif",size = 13, lineheight = 0.9, vjust = 1))+
  theme(axis.title.y = element_text(family="serif",vjust=1.2, size=15, color = "#696969"))+
  theme(axis.title.x = element_text(family="serif",vjust=-0.5, size=15)) 

detach(df_bio6_theta)

print(theta_bio6)








sigma_allgroups <- read.csv("df_bio6_alphasigma_bio6_Cal_2.csv")

## Box plot
attach(sigma_allgroups)
head(sigma_allgroups)

sigma_allgroups <- ggplot(sigma_allgroups, aes(x=elevation, y=sigma.sq, color = elevation, alpha = 0.8)) +
  geom_boxplot(aes(color = elevation), size = .75, width=0.8) +
  #geom_boxplot(aes(color = elevation))+
  geom_jitter(alpha = 0.5,  size = 5, position=position_jitter(0.1)) +
  scale_y_continuous(breaks=seq(-0.05, 0.20, 0.05), limit = c(-0.05, 0.20)) +
  scale_color_manual(values=c("#1E52F6","#D2D2D2"),
                     name  = "elevation",
                     
                     breaks = c("highland","lowland"),
                     
                     labels = c("Highland","Lowland")) +
  #4E84C4
  scale_fill_manual(values=c("#1E52F6","#D2D2D2"),
                    name  = "elevation",
                    
                    breaks = c("highland","lowland"),
                    
                    labels = c("Highland","Lowland")) +
  
  
  
  # theme_bw(base_size = 12, base_family = "")+
  theme_classic()+
  
  theme(legend.position = "none")+
  
  
  #theme(panel.border = element_rect(colour = "black", size = 0.5))+
  
  theme(panel.grid.major = element_line(colour =NA),
        panel.grid.minor =element_line(colour = NA))+
  
  
  theme(axis.title.y= element_text( color = "black")) +
  
  labs(x="", y = "MAT (Bio 06; °C)<br><span style = 'font-size:13pt'><span style ='color:#696969;'> 
         Change per million years (°C)</span></span>",
       title = expression(Evolutionary~Rate~(sigma^2)),
       subtitle = "All groups")  +
  
  theme(axis.title.y = element_textbox_simple(
    width = NULL, family="serif", size=18, hjust = 0.5,
    orientation = "left-rotated")) +
  
  scale_x_discrete(limits=c("highland","lowland"),  
                   
                   labels= c("Highland","Lowland"))+
  
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        plot.title    = element_text(family = "serif", size = 15),
        plot.subtitle = element_text(family = "serif", size = 12, color = "#696969")) +
  
  theme(axis.text.x = element_text(family="serif", size = 13, lineheight = 0.9, vjust =0.5, angle = 10))+
  theme(axis.text.y = element_text(family="serif",size = 13, lineheight = 0.9, vjust = 1))+
  #theme(axis.title.y = element_text(family="serif",vjust=1.2, size=15, color = "#696969"))+
  theme(axis.title.x = element_text(family="serif",vjust=-0.5, size=15)) 



#pdf("for_evolution_sigma_theta_bio6.pdf", height = 6, width = 10)

plot_grid(sigma_allgroups, theta_bio6, labels = c("", ""), 
          ncol = 2, nrow = 1, label_fontfamily = "serif")

#dev.off()  
