#Start-up
setwd("C:\\Users\\geekb\\Documents\\School\\Internships\\F022\\Report\\Treponema_Serology\\Validation_Treponema_Pallidum\\Analysis")

library(ModelGood)
library(openxlsx)
library(dplyr)
library(plyr)
library(ggplot2)
library(ROCit)
library(ggpol)
library(scales)
library(tidyr)
library(ggalt)
library(magrittr)
library(data.table)
library(grid)

options(scipen=5)


rawmfiDF <- read.xlsx("SAS2.xlsx", sheet = "Sera", colNames = TRUE)


rawmfiDF$`BK VP1` <- rawmfiDF$`_62`
rawmfiDF$`JC VP1` <- rawmfiDF$`_63`
rawmfiDF$`WU VP1` <- rawmfiDF$`_64`
rawmfiDF$`pGP-3` <- rawmfiDF$`_65`
rawmfiDF$`mgG unique` <- rawmfiDF$`_66`
rawmfiDF$`MgPa N-Term` <- rawmfiDF$`_67`
rawmfiDF$`MgPa` <- rawmfiDF$`_70`
rawmfiDF$`GST-Tag` <- rawmfiDF$`_90`
rawmfiDF$Tp15 <- rawmfiDF$`_86`
rawmfiDF$Tp17 <- rawmfiDF$`_87`
rawmfiDF$Tp44 <- rawmfiDF$`_88`
rawmfiDF$Tp47 <- rawmfiDF$`_89`
rawmfiDF$Dilution <- factor(rawmfiDF$Dilution)
rawmfiDF["status"] <- lapply(rawmfiDF["Status"], as.numeric)
rawmfiDF$Status <- factor(rawmfiDF$status>0,levels=c(FALSE,TRUE),labels=c("Negative","Positive"))

mfiCols <- names(select(rawmfiDF,contains("_")))

# Set dilution for mfiDF
mfiDF <- rawmfiDF
mfiDF <- rawmfiDF[which(rawmfiDF$Dilution == "100"),]
cutoff100 <- c(280,2400,330,1500)

mfiDF <- rawmfiDF[which(rawmfiDF$Dilution == "1000"),]
cutoff1000 <- c(50,250,40,190)
sum(mfiDF$Status == 1)


## Scatterplot with colors for positive and negative
antigen <- "Tp15"
antigen <- "Tp17"
antigen <- "Tp44"
antigen <- "Tp47"

ggplot(mfiDF, aes_string(x="FortNr", y=antigen, color="Status")) + geom_point() + scale_y_log10()
ggplot(mfiDF, aes_string(x = antigen)) + geom_histogram() + scale_x_log10()


## Sensitivity and specificity stats

cut <- 250

antigen <- mfiDF$Tp17

mfiDF$mficut <- factor(antigen>cut,levels=c(TRUE,FALSE),labels=c("Positive","Negative"))

with(mfiDF,table(Cutoff=mfiDF$mficut,Syphilis=mfiDF$status))
with(mfiDF,Sensitivity(mficut,status))
with(mfiDF,Specificity(mficut,status))
with(mfiDF,NPV(mficut,status))
with(mfiDF,PPV(mficut,status))



## Optimal


mfiDF100 <- mfiDF[which(rawmfiDF$Dilution == "100"),]
mfiDF1000 <- mfiDF[which(rawmfiDF$Dilution == "1000"),]


mfiDF100$Tp15pos100 <- as.numeric(mfiDF100$Tp15>cutoff100[1])
mfiDF100$Tp17pos100 <- as.numeric(mfiDF100$Tp17>cutoff100[2])
mfiDF100$Tp44pos100 <- as.numeric(mfiDF100$Tp44>cutoff100[3])
mfiDF100$Tp47pos100 <- as.numeric(mfiDF100$Tp47>cutoff100[4])
mfiDF100 %<>% select(ID, Status,Tp15pos100,
                     Tp17pos100,
                     Tp44pos100,
                     Tp47pos100
                     )

mfiDF1000$Tp15pos1000 <- as.numeric(mfiDF1000$Tp15>cutoff1000[1])
mfiDF1000$Tp17pos1000 <- as.numeric(mfiDF1000$Tp17>cutoff1000[2])
mfiDF1000$Tp44pos1000 <- as.numeric(mfiDF1000$Tp44>cutoff1000[3])
mfiDF1000$Tp47pos1000 <- as.numeric(mfiDF1000$Tp47>cutoff1000[4])
mfiDF1000 %<>% select(ID, Tp15pos1000,
                      Tp17pos1000,
                      Tp44pos1000,
                      Tp47pos1000
                      )




whichAntigens <- c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE)

TpOpt <- inner_join(mfiDF100, mfiDF1000, by = "ID")

TpOpt <- TpOpt[,whichAntigens,drop=FALSE]

TpOpt$mficut <- factor(rowSums(select(TpOpt,contains("Tp")))>=2
                       ,levels=c(TRUE,FALSE),labels=c("+","-"))

TpOpt$FalseResult <- revalue(TpOpt$mficut, c("+"="Positive","-"="Negative"))

TpOpt$FalseResult <- TpOpt$FalseResult == TpOpt$Status

with(TpOpt,table(Cutoff=mficut,Syphilis=Status))
with(TpOpt,Sensitivity(mficut,Status))
with(TpOpt,Specificity(mficut,Status))
with(TpOpt,NPV(mficut,Status))
with(TpOpt,PPV(mficut,Status))

TpOpt %>% filter(FalseResult == FALSE)



mfiDF$mficut <- factor(antigen>cut,levels=c(TRUE,FALSE),labels=c("Positive","Negative"))




mfiDF$FortNr[which(mfiDF$mficut != mfiDF$status)]

## ROC curves


curve <- Roc(list(Tp15=Status~Tp15,Tp17=Status~Tp17,Tp44=Status~Tp44,Tp47=Status~Tp47),data=mfiDF)
curve <- Roc(Status~Tp15,data=mfiDF)
curve <- Roc(Status~Tp17,data=mfiDF)
curve <- Roc(Status~Tp44,data=mfiDF)
curve <- Roc(Status~Tp47,data=mfiDF)
plot(curve,auc=TRUE)

curve


## Sensitivity and specificity

##ROCit

antigen <- mfiDF$Tp15
sensMin <- 0.0

mfiDF$status <- factor(mfiDF$Status>0,levels=c(FALSE,TRUE),labels=c("Negative","Positive"))

roc <- rocit(score = antigen, class = mfiDF$status,negref = "Negative") 
summary(roc)

sumCut <- as.data.frame(cbind(Cutoff=roc$Cutoff, TPR=roc$TPR, TNR=1-roc$FPR, Youdan=roc$TPR+(1-roc$FPR)-1))
sumCut

which.max(sumCut$Youdan)

tail(subset(sumCut, sumCut$TNR>sensMin))

plot(roc, values = F, col = c(2,4))

# Box plots

boxplot(Tp15~Status,data=mfiDF, 
        xlab="Disease Status", ylab="MFI") 

p <- ggplot(longDF, aes(x=Beadset, y=MFIvalue))
p <- ggplot(mfiDF, aes(x=Status, y=Tp17))
p <- ggplot(mfiDF, aes(x=Status, y=Tp44))
p <- ggplot(mfiDF, aes(x=Status, y=Tp47))

p + geom_boxplot(fill=NA,width=0.2,outlier.color = NA,) + scale_y_log10() +labs(y= "Antigen MFI", x = "Reference Assay") +
    geom_jitter(shape=16, position=position_jitter(0.15))


#Grouped box-scatter plots
longDF <- gather(mfiDF, key="Beadset",value="MFIvalue",contains("Tp"), -FortNr)


cutoff_line <- data.frame(Beadset = c("Tp15", "Tp17", "Tp44","Tp47"), val = cutoff[-1])


ggplot(longDF, aes(x = Status, 
                   y = MFIvalue)) +
  geom_hline(aes(yintercept = val), cutoff_line,color='red', size=1.5)+
  theme_minimal()+
  geom_jitter( width = .15,
                      shape = 21, 
                      #jitter.color = NA,
                     #position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8) ,
                 alpha=.7, size=2, stroke = 1, color="#000000", fill="#444444"
                 ) +
  
  
  facet_grid(.~ Beadset, scales="free",space="free")+
  
  scale_y_log10() +
  annotation_logticks(sides = "l", outside = TRUE) + 

  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust=0.5, vjust=0.5),
        axis.text.x = element_text(size = 14),
        plot.caption = element_text(vjust = 5),
        plot.title = element_text(hjust = 0.5,vjust=0, size=12)
        )+
labs(y= "MFI", x = "Serostatus")

 
# Sample ages







# Grouped bar chart for sensitivity and specificity

sensspecDF <- read.xlsx("Cutoff Stats.xlsx", sheet = "Raw", colNames = TRUE)

ssDF <- sensspecDF[which(sensspecDF$Dilution == "100"),]
ssDF <- sensspecDF[which(sensspecDF$Dilution == "1000"),]

ggplot(ssDF, aes(x=Antigen, y=Value, fill=Stat)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2,position=position_dodge(.9)) +
  coord_cartesian(ylim = c(85, 100))+
  scale_fill_manual(values = c("#abcbff", "#0047b9")) +labs(y= "Percentage (%)", x = "Antigen (1:1000)")


# Seroprevalance of all antigens


cutoff1000 <- c(250,250,250,95,180,150,150,250)
antigens <- c("BK VP1","JC VP1","WU VP1","pGP-3","mgG unique","MgPa N-Term","MgPa","Tp17")

seromfiDF <- rawmfiDF[which(rawmfiDF$Dilution == "1000"),]
seromfiDF <- seromfiDF[antigens]

seromfiDF[antigens] <- as.numeric(seromfiDF[antigens] > cutoff1000)

seropos <- 100*colSums(seromfiDF)/nrow(seromfiDF)
posPercent <- data.frame(antigens,seropos)
posPercent %<>% rename(Antigen = antigens)
posPercent %<>% rename(Seropositivity = seropos)
posPercent <- posPercent[-c(7), ] 
posPercent$Antigen[5] <- "mgG\nunique"
posPercent$Antigen[6] <- "MgPa N-Term;\nrMgPa"
antigens <- antigens[-c(7)]
posPercent$Population <- "Patient"

posPercent$Antigen  <- factor(posPercent$Antigen , levels = posPercent$Antigen)

pop <- data.frame(posPercent$Antigen,c(90,34,89,4,12.8,6.5,0.007),rep("Public",7))
pop %<>% setNames(colnames(posPercent)) 

posPercent %<>% rbind(pop)

posPercent$Population  <- factor(posPercent$Population , levels = c("Public", "Patient"))

ggplot(data = posPercent, aes(x=Antigen,y=Seropositivity,fill=Population)) +
geom_bar(stat="identity", position="dodge", color="black",width=0.7) +
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust=0.5, vjust=0.5),
        plot.caption = element_text(vjust = 5),
        plot.title = element_text(hjust = 0.5,vjust=0, size=12),
        legend.position = c(.9,.8),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
  )+
  #annotation_logticks(sides = "l", outside = TRUE) + 
  annotation_ticks(sides = 'l', ticks_per_base = c(20))+
  scale_fill_manual(values=c("#ffffff","#AAAAAA")) +
  labs(y= "Seroprevalence (%)", x = "Antigen")


publicSera <- posPercent


###  Seroprevalance between patients






cutoff1000 <- c(250,250,250,95,180,150,150,250)
antigens <- c("BK VP1","JC VP1","WU VP1","pGP-3","mgG unique","MgPa N-Term","MgPa","Tp17")

seromfiDF <- rawmfiDF[which(rawmfiDF$Dilution == "1000"),]
seromfiDF <- seromfiDF[antigens]

seromfiDF[antigens] <- as.numeric(seromfiDF[antigens] > cutoff1000)

seropos <- filter(seromfiDF, seromfiDF$Tp17 == 1)
seroneg <- filter(seromfiDF, seromfiDF$Tp17 == 0)

seropos <- 100*colSums(seropos)/nrow(seropos)
seropos <- data.frame(antigens,seropos)
seropos$Status <- "Positive"
seropos %<>% dplyr::rename(Seropositivity = seropos)
seropos <- seropos[-c(7,8), ] 
seropos$antigens[5] <- "mgG\nunique"
seropos$antigens[6] <- "MgPa N-Term;\nrMgPa"

seropos$antigens  <- factor(seropos$antigens , levels = seropos$antigens)

seroneg <- 100*colSums(seroneg)/nrow(seroneg)
seroneg <- data.frame(antigens,seroneg)
seroneg$Status <- "Negative"
seroneg %<>% dplyr::rename(Seropositivity = seroneg)
seroneg <- seroneg[-c(7,8), ] 
seroneg$antigens[5] <- "mgG\nunique"
seroneg$antigens[6] <- "MgPa N-Term;\nrMgPa"

posPercent <- rbind(seropos, seroneg)

posPercent %<>% dplyr::rename(Antigen = antigens)

posPercent %<>% dplyr::rename(Serostatus = Status)

ggplot(data = posPercent, aes(x=Antigen,y=Seropositivity,fill=Serostatus)) +
  geom_bar(stat="identity", position="dodge", color="black",width=0.7) +
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black",fill = NA, size=1),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust=0.5, vjust=0.5),
        plot.caption = element_text(vjust = 5),
        plot.title = element_text(hjust = 0.5,vjust=0, size=12),
        legend.position = c(.87,.8),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
  )+
  #annotation_logticks(sides = "l", outside = TRUE) + 
  annotation_ticks(sides = 'l', ticks_per_base = c(20))+
  scale_fill_manual(values=c("#ffffff","#AAAAAA")) +
  labs(y= "Seroprevalence (%)", x = "Antigen")






patientSera <- posPercent




###### Demographics

demoDF <- read.xlsx("T.pallidum_ReferenceSera.xlsx", sheet = "Filtered_Selection", colNames = TRUE)
demoDF <- demoDF[complete.cases(demoDF), ]

agebreaks <- c(0,20,30,40,50,60,70,80,100)
agelabels <- c("0-20","21-30","31-40","41-50","51-60","61-70","71-80","80+")

setDT(demoDF)[ , agegroups := cut(ALTER, 
                                breaks = agebreaks, 
                                right = FALSE, 
                                labels = agelabels)]




ggplot(data = demoDF, aes(x=agegroups,fill=ERGEBNIST1)) +
  geom_histogram(stat="count", position="dodge", color = "Black", width=0.7) +
  
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black",fill = NA, size=1),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust=0.5, vjust=0.5),
        plot.caption = element_text(vjust = 5),
        plot.title = element_text(hjust = 0.5,vjust=0, size=12),
        legend.position = c(.945,.8),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
  )+
  #stat_density(aes(group = agegroups),position="identity",geom="line") +
  facet_wrap(~ SEX) + 
  #annotation_logticks(sides = "l", outside = TRUE) + 
  annotation_ticks(sides = 'l', ticks_per_base = c(16))+
  scale_fill_manual(values=c("#ffffff","#AAAAAA"),name = "Serostatus") +
  labs(y= "Count", x = "Age Group")



mean(filter(demoDF, demoDF$SEX=="Female")$ALTER)


count(filter(demoDF, demoDF$SEX=="Male", demoDF$ERGEBNIST1=="Positive"))
count(filter(demoDF, demoDF$SEX=="Male", demoDF$ERGEBNIST1=="Negative"))
count(filter(demoDF, demoDF$SEX=="Female", demoDF$ERGEBNIST1=="Positive"))
count(filter(demoDF, demoDF$SEX=="Female", demoDF$ERGEBNIST1=="Negative"))

################### ELISA


measDF <- read.xlsx("200627 - ELISA results.xlsx", sheet = "Data", colNames = TRUE,check.names = FALSE)

measDF <- gather(measDF, Protein, Measurement, `400`:`0.01`, factor_key=TRUE )

measDF$Protein <- as.numeric(as.character(measDF$Protein))

errorDF <- read.xlsx("200627 - ELISA results.xlsx", sheet = "Error", colNames = TRUE,check.names = FALSE)
errorDF <- gather(errorDF, Protein, Measurement, `400`:`0.01`, factor_key=TRUE )
errorDF$Protein <- as.numeric(as.character(errorDF$Protein))
errorDF %<>% dplyr::rename(Error = Measurement)


elisaDF <- merge(measDF, errorDF, by=c("Antigen","Protein")) # NA's match
elisaDF %<>% arrange(Protein)


dodge <- 0.15

ggplot(elisaDF, aes(x=Protein,y=Measurement,color=factor(Antigen),shape=factor(Antigen))) + 
  
  geom_point(size = 3,position = position_dodge(width = dodge)) + 
  scale_x_log10(breaks=c(0.01, 0.1, 1, 10, 100),labels=c(0.01, 0.1, 1, 10, 100)) +
  annotation_logticks(sides = "b") + 
  geom_path(size = 1, position = position_dodge(width = dodge)) +
  #geom_smooth(method = "loess", se = FALSE) + 
  annotation_ticks(sides = 'l', ticks_per_base = c(20))+
  labs(y= "A480", x = "Protein (µg/well)") + 
   geom_errorbar(aes(ymin=Measurement-Error, ymax=Measurement+Error), width=.2,
                 position = position_dodge(width = dodge)) + 
  scale_fill_manual(values=c("#ffffff","#AAAAAA"))+
theme_minimal()+
  theme(panel.border = element_rect(colour = "black",fill = NA, size=1),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust=0.5, vjust=0.5),
        plot.caption = element_text(vjust = 5),
        plot.title = element_text(hjust = 0.5,vjust=0, size=12),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        legend.title = element_blank()
  )


labels=c(0.01, 0.1, 1, 10, 100)
