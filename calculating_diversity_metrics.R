setwd("C:/TFE/Week 10 Biodiversity")

############ MEASURING COMPOSITIONAL SIMILARITY #########

install.packages("safedata")
install.packages("vegan")
install.packages("vegetarian")
library(safedata)
set_safe_dir("C:/TFE/safedata")
library(vegan)
library(vegetarian)


###### getting data ########

culicid<- search_taxa("Culicidae") #search taxa
brant <- search_authors("Brant") #search authour
combined <- culicid & brant #combine search to reveal single record
combined

mozzies <- load_safe_data(record_id = combined, worksheet ='DailyHLC2013-2014')
attach(mozzies)

###### calculating Jaccard similarity ######

ground <- mozzies[Disturbance=='Primary'&Height=='Ground', -(1:12) ]
canopy <- mozzies[Disturbance=='Primary'&Height=='Canopy', -(1:12) ]

ground.mega <- colSums(ground)
canopy.mega <- colSums(canopy)
full.mega <- ground.mega+canopy.mega

full.names <- names(full.mega)[full.mega != 0]

a <- length(names(ground.mega)[ground.mega > 0 & canopy.mega > 0])
b <- length(names(ground.mega)[ground.mega > 0 & canopy.mega == 0])
c <- length(names(canopy.mega)[ground.mega == 0 & canopy.mega > 0])

J <- a/(a+b+c)
J

####### shortway to get jaccard ##########

comm <- rbind(ground.mega, canopy.mega)

1 - vegdist(comm, method = "jaccard", binary = TRUE)


###### 26 options in vegdist function #######
indices <- c("canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "chao")

similarity <- index <- type <- NULL
for(i in 1:length(indices)){
  similarity <- c(similarity, 1 - vegdist(comm, method = indices[i], binary = TRUE))
  index <- c(index, indices[i])
  type <- c(type, "incidence")
  similarity <- c(similarity, 1 - vegdist(comm, method = indices[i], binary = FALSE))
  index <- c(index, indices[i])
  type <- c(type, "abundance")
}

summary(similarity)
hist(similarity)

######## calculating similarity profiles ########

ground.mat <- t(as.matrix(ground.mega))
canopy.mat <- t(as.matrix(canopy.mega))
sim.groups(ground.mat, canopy.mat, q=0)

sorenson <- 1 - vegdist(comm, method = "bray", binary = TRUE)
sorenson
##### sim.groups and vegdist for sorenson give same answer 

C1N <- sim.groups(ground.mat, canopy.mat, q=1)
C2N <- sim.groups(ground.mat, canopy.mat, q=2)
C1N
C2N

####### partitioning diversity ###########

?d
alpha <- d(comm, lev="alpha", q=0)

alpha.ground <- d(ground.mega, lev="alpha", q=0)
aplha.canopy <- d(canopy.mega, lev="alpha", q=0)

gamma <- d(comm, lev="gamma", q=0)
beta <- d(comm, lev="beta", q=0)

beta.1 <- d(comm, lev="beta", q=1)
beta.2 <- d(comm, lev="beta", q=2)
