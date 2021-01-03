setwd("C:/Sorrell Cowen")
Ad_Lib <- read.csv("C:/Sorrell Cowen/Data/AAd_Libs_JC_2019.csv")
Census <- read.csv("C:/Sorrell Cowen/Data/CCensus_JC_2019.csv")
HOCensus_CS <- Census[Census$Status==c("IN CS"), c("AnimalID","Age","Sex","Group","Status")]

firstmerge <- merge(Ad_Lib,HOCensus_CS,by="AnimalID")
View(firstmerge)
firstmergee <- subset(firstmerge, select = -c(Obs,Notes,IOR,Response,Actor2,Actor3,Recip2,Recip3))
summary(firstmergee)

no_FUN <- firstmergee[firstmergee$Recipient!="FUN", ]  ####firstmergee is subset data for ad libs and census###
summary(no_FUN)

####graph for female mount proximity####
summary(no_FUN)
no_FUN_Mount <- no_FUN[no_FUN$Action=="Mount", ]
no_FUN_Mount <- no_FUN_Mount[!(no_FUN_Mount$Prox_10M_AF=="No Value"), ]
no_FUN_Mount <- no_FUN_Mount[!(no_FUN_Mount$Prox_10M_AF=="UN"), ]
subno_FUN_Mount_ProxF <- subset(no_FUN_Mount, select = c(AnimalID,Action,Context,Prox_10M_AF,Prox_10M_All))
summary(subno_FUN_Mount_ProxF)
unique(subno_FUN_Mount_ProxF$Prox_10M_AF)
#####Remove Levels####
is.na(subno_FUN_Mount_ProxF$Prox_10M_AF) <- subno_FUN_Mount_ProxF$Prox_10M_AF == "No Value"
subno_FUN_Mount_ProxF$Prox_10M_AF <- factor(subno_FUN_Mount_ProxF$Prox_10M_AF)
is.na(subno_FUN_Mount_ProxF$Prox_10M_AF) <- subno_FUN_Mount_ProxF$Prox_10M_AF == "UN"
subno_FUN_Mount_ProxF$Prox_10M_AF <- factor(subno_FUN_Mount_ProxF$Prox_10M_AF)
###reorder levels####
subno_FUN_Mount_ProxF$Prox_10M_AF <- factor(subno_FUN_Mount_ProxF$Prox_10M_AF,levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15+"))

plot(subno_FUN_Mount_ProxF$Prox_10M_AF,xlim=c(0,20), ylim=c(0,300), ylab= "Frequency of \nMounts", xlab= "No. of Females within 10m of Homosexual Mount", col="#999999", cex.names ="1")

#####graph for unknown female mounting female proximity####
FUN <- firstmergee[firstmergee$Recipient=="FUN", ]
FUN_Mount <- FUN[FUN$Action=="Mount", ]
FUN_Mount <- FUN_Mount[!(FUN_Mount$Prox_10M_AF=="No Value"), ]
FUN_Mount <- FUN_Mount[!(FUN_Mount$Prox_10M_AF=="UN"), ]
sub_FUN_Mount_ProxF <- subset(FUN_Mount, select = c(AnimalID,Action,Context,Prox_10M_AF,Prox_10M_All))
summary(sub_FUN_Mount_ProxF)
#####Remove Levels####
is.na(sub_FUN_Mount_ProxF$Prox_10M_AF) <- sub_FUN_Mount_ProxF$Prox_10M_AF == "No Value"
sub_FUN_Mount_ProxF$Prox_10M_AF <- factor(sub_FUN_Mount_ProxF$Prox_10M_AF)
is.na(sub_FUN_Mount_ProxF$Prox_10M_AF) <- sub_FUN_Mount_ProxF$Prox_10M_AF == "UN"
sub_FUN_Mount_ProxF$Prox_10M_AF <- factor(sub_FUN_Mount_ProxF$Prox_10M_AF)
###reorder levels####
sub_FUN_Mount_ProxF$Prox_10M_AF <- factor(sub_FUN_Mount_ProxF$Prox_10M_AF,levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15+"))

plot(sub_FUN_Mount_ProxF$Prox_10M_AF, ylim=c(0,300), ylab= "Frequency of \nMounts", xlab= "No. of Females within 10m of Heterosexual Mount", col="#999999", cex= "1")

summary(subno_FUN_Mount_ProxF)

ex_subnoFUN <- subno_FUN_Mount_ProxF[ , c("AnimalID","Prox_10M_AF")]
ex_subFUN <- sub_FUN_Mount_ProxF[ , c("AnimalID","Prox_10M_AF")]

ex_subFUN <- na.omit(ex_subFUN)
ex_subnoFUN <- na.omit(ex_subnoFUN)

#### need to change levels so 15+ is 15 in order to gain mean ####
library(plyr)

levels(ex_subFUN$Prox_10M_AF)
levels(ex_subFUN$Prox_10M_AF)[levels(ex_subFUN$Prox_10M_AF)=="15+"] <- 15
ex_subFUN[,2] <- as.numeric(as.character(ex_subFUN[,2]))
View(ex_subFUN)

levels(ex_subnoFUN$Prox_10M_AF)
levels(ex_subnoFUN$Prox_10M_AF)[levels(ex_subnoFUN$Prox_10M_AF)=="15+"] <- 15
ex_subnoFUN[,2] <- as.numeric(as.character(ex_subnoFUN[,2]))
View(ex_subnoFUN)

t.test(ex_subnoFUN$Prox_10M_AF,ex_subFUN$Prox_10M_AF)

mean(ex_subFUN$Prox_10M_AF)
sd(ex_subFUN$Prox_10M_AF)
mean(ex_subnoFUN$Prox_10M_AF)
sd(ex_subnoFUN$Prox_10M_AF)



#### producing graphs of linear regression to compare mount frquency against frequency of females observed around each individual male ####

byID <- group_by(zcanID, ID)
byIDF <- summarise(byID, mean_F = mean(F, na.rm = TRUE))
byIDF

meanf <- tapply(zcanID$F, zcanID$ID, mean, na.rm=TRUE)
meanf
class(meanf)
dfmeanf <- as.data.frame.table(meanf)
View(dfmeanf)
colnames(dfmeanf) <- c("AnimalID","Mean_F")

summary(ex_subnoFUN)

ybaby <- count(ex_subnoFUN, AnimalID)
View(ybaby)

reghome <- merge(dfmeanf, ybaby, by="AnimalID", all=TRUE)
reghome[is.na(reghome)] <- 0

colnames(reghome) <- c("AnimalID","Mean_F_Prox","Freq_of_Homo_Mount")
reggieboy <- lm(reghome$Mean_F_Prox ~ reghome$Freq_of_Homo_Mount)

plot(reghome$Freq_of_Homo_Mount, reghome$Mean_F_Prox)
abline(lm(reghome$Mean_F_Prox ~ reghome$Freq_of_Homo_Mount), col="red")



plot(reghome$Mean_F_Prox, reghome$Freq_of_Homo_Mount)
abline(lm(reghome$Freq_of_Homo_Mount ~ reghome$Mean_F_Prox), col="red")
opp <- lm(reghome$Freq_of_Homo_Mount ~ reghome$Mean_F_Prox)
summary(opp)


new_homo <- reghome[reghome$Mean_F_Prox + reghome$Freq_of_Homo_Mount > 0, ]


par(mar = c(5, 5, 5, 1))
plot(new_homo$Freq_of_Homo_Mount, new_homo$Mean_F_Prox, frame.plot= FALSE, pch=21, bg="#FFFFFF",xlim=c(0,60),ylim=c(0,0.8),ylab="Average Number of Females within 2m \nof Male Across All Behaviour" , xlab="Frequency the Male was Observed as the Actor in a Homosexual Mount")



abline(lm(new_homo$Mean_F_Prox ~ new_homo$Freq_of_Homo_Mount), col="black")
newhomoooo <- lm(new_homo$Mean_F_Prox ~ new_homo$Freq_of_Homo_Mount)
summary(newhomoooo)

library(dplyr)

##### including frequency they were recipient also ####

reciphomo <- no_FUN_Mount[ ,c("Recipient","Prox_10M_AF")]
counthomo <- count(reciphomo, Recipient)
counthomo <- na.omit(counthomo)
counthomo <- counthomo[-1, ]  ##### was a blank #####

counthomo <- counthomo %>%
  rename( AnimalID = Recipient
  )
recki <- merge(new_homo, counthomo, by="AnimalID", all=TRUE)
recki[is.na(recki)] <- 0
recki$Freq_Mount <- recki$Freq_of_Homo_Mount + recki$n
reckit <- recki[ ,c("AnimalID","Mean_F_Prox","Freq_Mount")]

reckit <- reckit[reckit$Mean_F_Prox + reckit$Freq_Mount > 0, ]

plot(reckit$Freq_Mount, reckit$Mean_F_Prox, frame.plot= FALSE, pch=21, bg="#FFFFFF",xlim=c(0,120),ylim=c(0,1),ylab="Average Number of Females within 2m \nof Male Across All Behaviour" , xlab="Frequency the Male was Observed in a Homosexual Mount")
reckitlm <- lm(reckit$Mean_F_Prox ~ reckit$Freq_Mount)
abline(reckitlm, col="black")
summary(reckitlm)

summary(reckit)
summary(heterolm)
summary(ratio1)

##### transforming to log #####


reckit$logfreqmount <- log10(reckit$Freq_Mount + 1)
plot(reckit$logfreqmount, reckit$Mean_F_Prox, frame.plot= FALSE, pch=21, bg="#FFFF00",ylab="Average Number of Females within 2m \nof Male Across All Behaviour" , xlab="log(x+1) \nFrequency the Male was Observed in a Homosexual Mount", main="Proximity of Females to Males Involved in Homosexual Mounts")
logreckitlm <- lm(reckit$Mean_F_Prox ~ reckit$Freq_Mount) 
summary(logreckitlm)
abline(logreckitlm)


###### just for recipient ####

reciprecki <- recki[ ,c("AnimalID","Mean_F_Prox","n")]
reciprecki <- reciprecki[reciprecki$Mean_F_Prox + reciprecki$n > 0, ]

reciprecki$Mean_F_Prox = jitter(reciprecki$Mean_F_Prox)
reciprecki$n = jitter(reciprecki$n)

plot(reciprecki$n, reciprecki$Mean_F_Prox, frame.plot= FALSE, pch=21, bg="#FFCC00",xlim=c(0,60),ylim=c(0,1),ylab="Average Number of Females within 2m \nof Male Across All Behaviour" , xlab="Frequency the Male was Observed as a Recipient in a Homosexual Mount", main="Proximity of Females to Homosexually Mounted Males")
reciponlylm <- lm(reciprecki$Mean_F_Prox ~ reciprecki$n)
abline(reciponlylm, col="black")
summary(reciponlylm)


#### other attempts at merging actor and recipient ####

recipCensus <- HOCensus_CS %>% 
  rename( Recipient = AnimalID
  )
recipmerge <- merge(no_FUN_Mount,recipCensus, by="Recipient")
summary(recipmerge)
recipmerge <- recipmerge[ , c("Recipient","Action")]
yeahbaby <- count(recipmerge, Recipient)
yeahbaby <- yeahbaby %>% 
  rename( AnimalID = Recipient
  )

actandrecip <- merge(new_homo, yeahbaby, by="AnimalID", all=TRUE)
actandrecip[is.na(actandrecip)] <- 0

actandrecip$Freq_in_a_Mount <- actandrecip$Freq_of_Homo_Mount + actandrecip$n

thefinale <- actandrecip[ ,c("AnimalID","Mean_F_Prox","Freq_in_a_Mount")]

thefinale <- thefinale[-195, ]

par(mar = c(5, 5, 5, 1))
plot(thefinale$Freq_in_a_Mount, thefinale$Mean_F_Prox, frame.plot= FALSE, pch=21, bg="#00FF99",xlim=c(0,100),ylim=c(0,1),ylab="Average Number of Females within 2m \nof Male Across All Behaviour" , xlab="Frequency the Male was Observed in a Homosexual Mount", main="Proximity of Females to Actively Homosexual Males")
thefinalelm <- lm(thefinale$Mean_F_Prox ~ thefinale$Freq_in_a_Mount)

abline(thefinalelm, col="black")
summary(thefinalelm)


#### including data right from beginning gives same results ####

specialmerge <- merge(reghome, yeahbaby, by="AnimalID", all=TRUE)
specialmerge[is.na(specialmerge)] <- 0
specialmerge$join <- specialmerge$Freq_of_Homo_Mount + specialmerge$n
spackymerge <- specialmerge[ , c("AnimalID","Mean_F_Prox","join")]
spackymerge[is.na(spackymerge)] <- 0
spackymerge <- spackymerge[spackymerge$Mean_F_Prox + spackymerge$join > 0, ]

sm <- lm(spackymerge$Mean_F_Prox ~ spackymerge$join)
plot(spackymerge$join, spackymerge$Mean_F_Prox, frame.plot= FALSE, pch=21, bg="#00FF99",xlim=c(0,100),ylim=c(0,1),ylab="Average Number of Females within 2m \nof Male Across All Behaviour" , xlab="Frequency the Male was Observed in a Homosexual Mount", main="Proximity of Females to Actively Homosexual Males")
abline(sm, col="black")
summary(sm)


##### ratio plot ########

summary(reckit)
count_fun_mount <- count(FUN_Mount, AnimalID)
count_fun_mount <- na.omit(count_fun_mount)

summary(new_homo)
ratio1 <- merge(count_fun_mount, new_homo, by="AnimalID", all=TRUE)
ratio1[is.na(ratio1)] <- 0
colnames(ratio1) <- c("AnimalID","Freq_Hetero_Mount","Mean_F_Prox","Freq_Homo_Mount")

ratio1$ratio <- ratio1$Freq_Homo_Mount/(ratio1$Freq_Hetero_Mount + ratio1$Freq_Homo_Mount)
ratio1 <- na.omit(ratio1)


plot(ratio1$ratio, ratio1$Mean_F_Prox, frame.plot= FALSE, pch=21, bg="#FFFFFF",xlim=c(0,1),ylim=c(0,0.8),ylab="Average Number of Females within 2m \nof Male Across All Behaviour" , xlab="Homosexual Activity Ratio of Male \n(Homosexual Mounts : Total Mounts)")
ratio1lm <- lm(ratio1$Mean_F_Prox ~ ratio1$ratio)
abline(ratio1lm)
summary(ratio1lm)
summary(ratio1)

###### heterosexual lm plot to compare to count homo lm plot #######

heterolm <- ratio1[ ,c("AnimalID","Mean_F_Prox","Freq_Hetero_Mount")]
plot(heterolm$Freq_Hetero_Mount, heterolm$Mean_F_Prox, frame.plot= FALSE, pch=21, bg="#FFFFFF",xlim=c(0,120),ylim=c(0,0.8),ylab="Average Number of Females within 2m \nof Male Across All Behaviour" , xlab="Frequency Male was Observed in a Heterosexual Mount")
lmheterolm <- lm(heterolm$Mean_F_Prox ~ heterolm$Freq_Hetero_Mount)
abline(lmheterolm)
summary(lmheterolm)
summary(heterolm)
summary(reckit)


#### log graph #####

heterolm$logfreqmount <- log(heterolm$Freq_Hetero_Mount + 1)
plot(heterolm$logfreqmount, heterolm$Mean_F_Prox, frame.plot= FALSE, pch=21, bg="#009966",xlim=c(0,7),ylim=c(0,0.8),ylab="Average Number of Females within 2m \nof Male Across All Behaviour" , xlab="log(x+1) \nFrequency Male was Observed in a Heterosexual Mount", main="Proximity of Females to Male Heterosexual Mounting")
logheterolm <- lm(heterolm$Mean_F_Prox ~ heterolm$logfreqmount)
abline(lmheterolm, col="black")
summary(logheterolm)
