setwd("C:/TFE/Week 10 Biodiversity")
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

#### CREATING DIVERSITY PROFILES ####


##### long way to calculate shannon diversity #########

twice.comm <-mozzies[Disturbance=='Twice', -(1:12) ]
twice.mega <- colSums(twice.comm)

twice.prop <- twice.mega/(sum(twice.mega))
sum(twice.prop)

val <- twice.prop*log(twice.prop)
total <- sum(val, na.rm=TRUE)

H <- -(total)
H


########## easy way to calculate shannon diversity ########
diversity(twice.mega, index='shannon')

#### calculating shannon evenness ####

S.obs <- colSums(twice.comm) > 0
S <- sum(S.obs)

J <- H/log(S)
J

#################### Calculate Berger-Parker relative dominance #######

twice.mega

Nmax <- max(twice.mega) 
Ntot <- sum(twice.mega)

Crel <- Nmax/Ntot 
Crel

###### Calculate the % Singletons ###########

Ssingleton <- sum(twice.mega ==1)

Rsingleton <- Ssingleton/S

######### calculate fisher's alpha ########

alpha <- fisher.alpha(twice.mega)


####### true diversity #######

qD <- matrix(NA, nrow=nrow(twice.comm), ncol=3)

for(i in 1:nrow(twice.comm)){
  for(j in 1:3){
    qD[i,j] <- d(twice.comm[i, ], q=j-1)
  }
}

qD

qD <- data.frame(qD)
names(qD)<-c("D0", "D1", "D2")
attach(qD)
head(qD)

?specnumber
Sestimates <- specnumber(twice.comm)
Sestimates
plot(Sestimates, qD$D0)


shannon <- diversity(twice.comm, index='shannon')
plot(shannon, qD$D1)


simpson <- diversity(twice.comm, index='simpson')
plot(simpson, qD$D2)


######### diversity profiles #######

D0mega <- sum(qD$D0)
D1mega <- sum(qD$D1)
D2mega <- sum(qD$D2)

plot(c(D0mega, D1mega, D2mega), type = "b", xlab = "alpha", ylab = "Diversity")


####### test for heterogeneity among samples ##########

ground <- mozzies[Disturbance=="Twice"&Height=="Ground", -(1:12)]

canopy <- mozzies[Disturbance=="Twice"&Height=="Canopy", -(1:12)]

ground.mega <- colSums(ground)
canopy.mega <- colSums(canopy)

Lc <- H*(sum(ground.mega)+sum(canopy.mega))
Lc

H11 <- diversity(ground.mega, index='shannon')
H12 <- diversity(canopy.mega, index='shannon')

L2 <- sum(ground.mega) * H11 + sum(canopy.mega) * H12
L2

AICc <- 2*(Lc+S)
AIC2 <- 2*(L2+(2*S))

AICc
AIC2
AICc < AIC2

## so AIC2 is greater >2 bigger than AICc, so treating as two separate communities is best model ##

########### Comparing diversity profiles #############

qDground <- c(d(ground.mega, q = 0), d(ground.mega, q = 1), d(ground.mega, q = 2))
qDcanopy <- c(d(canopy.mega, q = 0), d(canopy.mega, q = 1), d(canopy.mega, q = 2))

plot(rep(0:2, 2), c(qDground, qDcanopy),col= rep(c(1,2), 1, each = 3), xlab = "alpha", ylab = "Diversity")


############ MEASURING COMPOSITIONAL SIMILARITY #########


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
