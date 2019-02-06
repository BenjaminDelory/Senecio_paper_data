#######################################################################################################################################################
#R code used in Delory et al (2019) The exotic species Senecio inaequidens pays the price for arriving late in temperate European grassland communities
#######################################################################################################################################################

#########################################################################
#Set the path to the Senecio_paper_data directory
#Example: "C:/Users/Benjamin Delory/Documents/GitHub/Senecio_paper_data"
#########################################################################

path<-"SET_WORKING_DIRECTORY_HERE"

##################
#Run the code :-)
##################

#Set working directory

setwd(path)

#Load R packages

library(RColorBrewer)
library(readr)
library(multcomp)
library(beanplot)
#library(jpeg)
library(effsize)
library(car)

#Create bootstrap function to calculate confidence intervals using the percentile method

CIbootstrap<-function(x, n=1000, CI=0.95){
  
  #x is a numeric vector containing sample values
  #n is the number of bootsrapping
  #CI is the confidence interval type
  
  if (TRUE %in% is.na(x)){x<-x[-which(is.na(x)==TRUE)]}
  
  bstrap<-c()
  
  for (i in 1:n){
    
    bsample<-sample(x, length(x), replace=TRUE)
    bestimate<-mean(bsample)
    bstrap<-c(bstrap, bestimate)}
  
  return(c(quantile(bstrap, (1-CI)/2), quantile(bstrap, CI+(1-CI)/2)))}

#Load  data
data <- read_delim("Data/Delory_et_al_2019_Senecio.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
data<-as.data.frame(data)
data$Arrival<-as.factor(data$Arrival)
data$Community<-as.factor(data$Community)
SA<-0.28*0.28 #Surface area of a pot
data$BGrasses<-data$BGrasses/SA #Express data in g/m²
data$BSenecio<-data$BSenecio/SA #Express data in g/m²
data$BLegumes<-data$BLegumes/SA #Express data in g/m²
View(data)

#Correlation between the number of Senecio individuals and the productivity of the invasive species
plot(data$NumSenecio, data$BSenecio, pch=16, las=1)
cor.test(data$NumSenecio, data$BSenecio)

#Relative abundance of Senecio
data$percSenecio<-100*data$BSenecio/(data$BSenecio+data$BLegumes+data$BGrasses)

#Relative abundance of Natives
data$percNatives<-100*(data$BLegumes+data$BGrasses)/(data$BSenecio+data$BLegumes+data$BGrasses)

#################################################
#################################################
#Analysis SENECIO: Early arrival - Synchronous 1
#################################################
#################################################

data1<-data[data$Arrival=="Early"|data$Arrival=="Synchronous1",]
data1$Arrival<-factor(data1$Arrival)
data1$Community<-factor(data1$Community)

#######################
#Stats biomass Senecio
#######################

dotchart(data1$BSenecio) #Possible outlier detected. Check for measurement error. Confirm that this is not a measurement error.

#Fit linear model with all observations
mod1<-lm(BSenecio~Community*Arrival, data1)
summary(mod1)
par(mfrow=c(3,2))
plot(mod1)
plot(cooks.distance(mod1), type="h")
anova(mod1)
outlierTest(mod1)

#Fit GLM
mod2<-glm(BSenecio~Community*Arrival, data=data1, family=Gamma(link="log"))
mod3<-glm(BSenecio~Community*Arrival, data=data1, family=Gamma(link="identity"))

mod4<-glm(BSenecio~Community*Arrival+NumSenecio, data=data1, family=Gamma(link="log"))
mod5<-glm(BSenecio~Community*Arrival+NumSenecio, data=data1, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3, mod4, mod5)
drop1(mod4, test="F")
anova(mod2, mod4, test="F")

#Mod4 does not seem to be a significant improvement of mod2 based on AIC value

E2<-resid(mod2, type="pearson") #Pearson residuals
F2<-fitted(mod2) #Fitted values (same scale as the response variable)
eta<-predict(mod2, type="link") #Fitted values (same scale as the predictor function: log scale)
par(mfrow=c(3,2))
plot(F2, E2, xlab="Fitted values", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(eta, E2, xlab="Eta", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(data1$NumSenecio, E2, xlab="Number of Senecio individual", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
boxplot(E2~Arrival, data=data1, ylab="Pearson residuals", xlab="Arrival")
abline(h=0, lty=2)
boxplot(E2~Community, data=data1, ylab="Pearson residuals", xlab="Community")
abline(h=0, lty=2)
plot(cooks.distance(mod2), type="h") #All values lower than 1, so no influential observation

E4<-resid(mod4, type="pearson") #Pearson residuals
F4<-fitted(mod4) #Fitted values (same scale as the response variable)
eta<-predict(mod4, type="link") #Fitted values (same scale as the predictor function: log scale)
par(mfrow=c(3,2))
plot(F4, E4, xlab="Fitted values", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(eta, E4, xlab="Eta", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(data1$NumSenecio, E4, xlab="Number of Senecio individual", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
boxplot(E4~Arrival, data=data1, ylab="Pearson residuals", xlab="Arrival")
abline(h=0, lty=2)
boxplot(E4~Community, data=data1, ylab="Pearson residuals", xlab="Community")
abline(h=0, lty=2)
plot(cooks.distance(mod4), type="h") #All values lower than 1, so no influential observation

mod2.nointer<-glm(BSenecio~Arrival+Community, data=data1, family=Gamma(link="log"))

summary(mod2)
drop1(mod2, test="F") #Significant Arrival*Community interaction
drop1(mod2.nointer, test="F")
DevianceExplained2<-100*(mod2$null.deviance-mod2$deviance)/mod2$null.deviance #% of the variation explained by the model

layout(1)
interaction.plot(x.factor=data1$Community, trace.factor=data1$Arrival, response=data1$BSenecio)

tmp<-expand.grid(Arrival=levels(data1$Arrival), Community=levels(data1$Community)) #We first compute the means of the response for all 8 combinations of the levels of Arrival and Community
X<-model.matrix(~Community*Arrival, data=tmp)
glht(mod2, linfct=X)
predict(mod2, newdata=tmp, type="link")
#We now construct a contrast matrix based on Tukey-contrasts for Arrival in a block-diagonal way (for each level of Community)
Tukey<-contrMat(table(data1$Arrival), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(data1$Community)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(data1$Community)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))
K<-rbind(K, matrix(0, nrow=2, ncol=4))
rownames(K)<-c(rownames(K)[1:2], "G:Early - GL:Early", "G:Synchronous1 - GL:Synchronous1")
K[3,]<-c(1,0,-1,0)
K[4,]<-c(0,1,0,-1)

summary(glht(mod2, linfct=K%*%X))

########################
#Stats N content senecio
########################

mod1<-lm(NSenecio~Community*Arrival, data1)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
anova(mod1) #No effect

########################
#Stats C content senecio
########################

mod1<-lm(CSenecio~Community*Arrival, data1)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
anova(mod1)

tmp<-expand.grid(Arrival=levels(data1$Arrival), Community=levels(data1$Community)) #We first compute the means of the response for all 8 combinations of the levels of Arrival and Community
X<-model.matrix(~Community*Arrival, data=tmp)
glht(mod1, linfct=X)
predict(mod1, newdata=tmp, type="response")
#We now construct a contrast matrix based on Tukey-contrasts for Arrival in a block-diagonal way (for each level of Community)
Tukey<-contrMat(table(data1$Arrival), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(data1$Community)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(data1$Community)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))
K<-rbind(K, matrix(0, nrow=2, ncol=4))
rownames(K)<-c(rownames(K)[1:2], "G:Early - GL:Early", "G:Synchronous1 - GL:Synchronous1")
K[3,]<-c(1,0,-1,0)
K[4,]<-c(0,1,0,-1)

summary(glht(mod1, linfct=K%*%X))

##############
#Effect sizes
##############

data1$treat<-paste(data1$Arrival, data1$Community, sep="-")

GEarly.GS1<-data1$BSenecio[data1$treat=="Early-G"|data1$treat=="Synchronous1-G"]
f<-as.factor(as.character(data1$treat[data1$treat=="Early-G"|data1$treat=="Synchronous1-G"]))
es.GEarly.GS1<-cohen.d(GEarly.GS1, f, hedges.correction=TRUE)

GLEarly.GLS1<-data1$BSenecio[data1$treat=="Early-GL"|data1$treat=="Synchronous1-GL"]
f<-as.factor(as.character(data1$treat[data1$treat=="Early-GL"|data1$treat=="Synchronous1-GL"]))
es.GLEarly.GLS1<-cohen.d(GLEarly.GLS1, f, hedges.correction=TRUE)

GEarly.GLEarly<-data1$BSenecio[data1$treat=="Early-G"|data1$treat=="Early-GL"]
f<-as.factor(as.character(data1$treat[data1$treat=="Early-G"|data1$treat=="Early-GL"]))
es.GEarly.GLEarly<-cohen.d(GEarly.GLEarly, f, hedges.correction=TRUE)

GS1.GLS1<-data1$BSenecio[data1$treat=="Synchronous1-G"|data1$treat=="Synchronous1-GL"]
f<-as.factor(as.character(data1$treat[data1$treat=="Synchronous1-G"|data1$treat=="Synchronous1-GL"]))
es.GS1.GLS1<-cohen.d(GS1.GLS1, f, hedges.correction=TRUE)

es.senecio1<-data.frame(Hedges=c(es.GEarly.GS1$estimate, es.GLEarly.GLS1$estimate, es.GEarly.GLEarly$estimate, es.GS1.GLS1$estimate),
                      CIinf=c(es.GEarly.GS1$conf.int[1], es.GLEarly.GLS1$conf.int[1], es.GEarly.GLEarly$conf.int[1], es.GS1.GLS1$conf.int[1]),
                      CIsup=c(es.GEarly.GS1$conf.int[2], es.GLEarly.GLS1$conf.int[2], es.GEarly.GLEarly$conf.int[2], es.GS1.GLS1$conf.int[2]))

rownames(es.senecio1)<-c("G:Early-G:Synch1", "GL:Early-GL:Synch1", "G:Early-GL:Early", "G:Synch1-GL:Synch1")

#################################################
#################################################
#Analysis SENECIO: Late arrival - Synchronous 2
#################################################
#################################################

data2<-data[data$Arrival=="Late"|data$Arrival=="Synchronous2",]
data2$Arrival<-factor(data2$Arrival)
data2$Community<-factor(data2$Community)

#######################
#Stats biomass Senecio
#######################

dotchart(data2$BSenecio)

#Fit linear model with all observations
mod1<-lm(BSenecio~Community*Arrival, data2)
summary(mod1)
par(mfrow=c(3,2))
plot(mod1) #We clearly have non-homogeneity of variances between groups
plot(cooks.distance(mod1), type="h")

mod2<-glm(BSenecio~Community*Arrival, data=data2, family=Gamma(link="log"))
mod3<-glm(BSenecio~Community*Arrival, data=data2, family=Gamma(link="identity"))

mod4<-glm(BSenecio~Community*Arrival+NumSenecio, data=data2, family=Gamma(link="log"))
mod5<-glm(BSenecio~Community*Arrival+NumSenecio, data=data2, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3, mod4, mod5)
drop1(mod4, test="F")
anova(mod2, mod4, test="F")

#Mod4 does not seem to be a significant improvement of mod2 based on AIC value

E2<-resid(mod2, type="pearson") #Pearson residuals
F2<-fitted(mod2) #Fitted values (same scale as the response variable)
eta<-predict(mod2, type="link") #Fitted values (same scale as the predictor function: log scale)
par(mfrow=c(3,2))
plot(F2, E2, xlab="Fitted values", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(eta, E2, xlab="Eta", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(data2$NumSenecio, E2, xlab="Number of Senecio individual", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
boxplot(E2~Arrival, data=data2, ylab="Pearson residuals", xlab="Arrival")
abline(h=0, lty=2)
boxplot(E2~Community, data=data2, ylab="Pearson residuals", xlab="Community")
abline(h=0, lty=2)
plot(cooks.distance(mod2), type="h") #All values lower than 1, so no influential observation

E4<-resid(mod4, type="pearson") #Pearson residuals
F4<-fitted(mod4) #Fitted values (same scale as the response variable)
eta<-predict(mod4, type="link") #Fitted values (same scale as the predictor function: log scale)
par(mfrow=c(3,2))
plot(F4, E4, xlab="Fitted values", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(eta, E4, xlab="Eta", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(data2$NumSenecio, E4, xlab="Number of Senecio individual", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
boxplot(E4~Arrival, data=data2, ylab="Pearson residuals", xlab="Arrival")
abline(h=0, lty=2)
boxplot(E4~Community, data=data2, ylab="Pearson residuals", xlab="Community")
abline(h=0, lty=2)
plot(cooks.distance(mod4), type="h") #All values lower than 1, so no influential observation

summary(mod2)
drop1(mod2, test="F") #No significant interaction between Arrival and Community
drop1(glm(BSenecio~Community+Arrival, data=data2, family=Gamma(link="log")), test="F") #Only Arrival is significant
par(mfrow=c(2,2))
plot(mod2) #By default, plot(glm object) uses the deviance residuals (and not the Pearson residuals)
DevianceExplained2<-100*(mod2$null.deviance-mod2$deviance)/mod2$null.deviance #% of the variation explained by the model

tmp<-expand.grid(Arrival=levels(data2$Arrival), Community=levels(data2$Community)) #We first compute the means of the response for all 8 combinations of the levels of Arrival and Community
X<-model.matrix(~Community*Arrival, data=tmp)
glht(mod2, linfct=X)
predict(mod2, newdata=tmp, type="link")
#We now construct a contrast matrix based on Tukey-contrasts for Arrival in a block-diagonal way (for each level of Community)
Tukey<-contrMat(table(data2$Arrival), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(data2$Community)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(data2$Community)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))
K<-rbind(K, matrix(0, nrow=2, ncol=4))
rownames(K)<-c(rownames(K)[1:2], "G:Late - GL:Late", "G:Synchronous2 - GL:Synchronous2")
K[3,]<-c(1,0,-1,0)
K[4,]<-c(0,1,0,-1)

summary(glht(mod2, linfct=K%*%X))

########################
#Stats N content senecio
########################

mod1<-lm(NSenecio~Community*Arrival, data2)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
anova(mod1) #Significant effect of Arrival
Anova(mod1, type="II")
Anova(mod1, type="III")

tmp<-expand.grid(Arrival=levels(data2$Arrival), Community=levels(data2$Community)) #We first compute the means of the response for all 8 combinations of the levels of Arrival and Community
X<-model.matrix(~Community*Arrival, data=tmp)
glht(mod1, linfct=X)
predict(mod1, newdata=tmp, type="response")
#We now construct a contrast matrix based on Tukey-contrasts for Arrival in a block-diagonal way (for each level of Community)
Tukey<-contrMat(table(data2$Arrival), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(data2$Community)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(data2$Community)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))
K<-rbind(K, matrix(0, nrow=2, ncol=4))
rownames(K)<-c(rownames(K)[1:2], "G:Late - GL:Late", "G:Synchronous2 - GL:Synchronous2")
K[3,]<-c(1,0,-1,0)
K[4,]<-c(0,1,0,-1)

summary(glht(mod1, linfct=K%*%X))

########################
#Stats C content senecio
########################

mod1<-lm(CSenecio~Community*Arrival, data2)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
anova(mod1)

#######################
#Nitrogen facilitation
#######################

newdata<-data[-which(data$Arrival=="Late"),]
model<-lm(deltaNSenecio~Community*Arrival, newdata)
plot(fitted(model), resid(model)); abline(h=0)
anova(model) #delta 15N not different between treatments
summary(glht(model, linfct=mcp(Community="Tukey")))
summary(glht(model, linfct=mcp(Arrival="Tukey")))

##############
#Effect sizes
##############

data2$treat<-paste(data2$Arrival, data2$Community, sep="-")

GLate.GS2<-data2$BSenecio[data2$treat=="Late-G"|data2$treat=="Synchronous2-G"]
f<-as.factor(as.character(data2$treat[data2$treat=="Late-G"|data2$treat=="Synchronous2-G"]))
es.GLate.GS2<-cohen.d(GLate.GS2, f, hedges.correction=TRUE)

GLLate.GLS2<-data2$BSenecio[data2$treat=="Late-GL"|data2$treat=="Synchronous2-GL"]
f<-as.factor(as.character(data2$treat[data2$treat=="Late-GL"|data2$treat=="Synchronous2-GL"]))
es.GLLate.GLS2<-cohen.d(GLLate.GLS2, f, hedges.correction=TRUE)

GLate.GLLate<-data2$BSenecio[data2$treat=="Late-G"|data2$treat=="Late-GL"]
f<-as.factor(as.character(data2$treat[data2$treat=="Late-G"|data2$treat=="Late-GL"]))
es.GLate.GLLate<-cohen.d(GLate.GLLate, f, hedges.correction=TRUE)

GS2.GLS2<-data2$BSenecio[data2$treat=="Synchronous2-G"|data2$treat=="Synchronous2-GL"]
f<-as.factor(as.character(data2$treat[data2$treat=="Synchronous2-G"|data2$treat=="Synchronous2-GL"]))
es.GS2.GLS2<-cohen.d(GS1.GLS1, f, hedges.correction=TRUE)

es.senecio2<-data.frame(Hedges=c(es.GLate.GS2$estimate, es.GLLate.GLS2$estimate, es.GLate.GLLate$estimate, es.GS2.GLS2$estimate),
                       CIinf=c(es.GLate.GS2$conf.int[1], es.GLLate.GLS2$conf.int[1], es.GLate.GLLate$conf.int[1], es.GS2.GLS2$conf.int[1]),
                       CIsup=c(es.GLate.GS2$conf.int[2], es.GLLate.GLS2$conf.int[2], es.GLate.GLLate$conf.int[2], es.GS2.GLS2$conf.int[2]))

rownames(es.senecio2)<-c("G:Late-G:Synch2", "GL:Late-GL:Synch2", "G:Late-GL:Late", "G:Synch2-GL:Synch2")

#################
#################
#Export Figure 2
#################
#################

colors<-c("white", "grey50", "black", "black")

#img<-readJPEG("senecio3crop.jpg")
#imgheight<-3648
#imgwidth<-1872
#imgratio<-imgheight/imgwidth

tiff(filename="Figure2.tif", res=1000, width=15.5, height=12, units="cm", compression="lzw", pointsize=12)
layout(matrix(c(1,1,2,3,4,5), ncol=3, nrow=2), widths=c(3,6.5,6.5))

par(mar=c(5,0,4,0)+0.1)
plot(1:10, ty="n", ylab="", xlab="", bty="n", xaxt="n", yaxt="n", asp=1)
#rasterImage(img,1.7,-2.5,9.7,-2.5+8*imgratio)

par(bty="l")
par(mar=c(2,5,2,0.5)+0.1)
gf<-beanplot(BSenecio~Arrival*Community, data1, log="", what=c(1,0,1,1), border = NA, col = colors, 
             ylab=expression(paste(italic("Senecio")," shoot dry weight (g ", m^-2, ")", sep="")), xlab="", 
             las=1, at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(0,700), names=c("","","",""), yaxt="n")
axis(2, at=seq(0,700,by=100), las=1)
text(1.5, 520, "***", cex=2)
text(4.25, 700, "***", cex=2)
lines(c(2.3, 2.3), c(gf$stats[2], 100), lty=2)
lines(c(4.45, 4.45), c(gf$stats[4], 100), lty=2)
lines(c(2.3, 4.45), c(100, 100), lty=2)
text((2.3+4.45)/2, 100, "*", pos=1, cex=2)
text(0.25, 700, "A", font=2, cex=1.4, pos=4, xpd=TRUE)

par(bty="l")
par(mar=c(4,5,0,0.5)+0.1)
gf<-beanplot(NSenecio~Arrival*Community, data1, log="", what=c(1,0,1,1), border = NA, col = colors, 
             ylab=expression(paste(italic("Senecio")," shoot N content (%)", sep="")), 
             xlab="", las=1, at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(2.5,7), 
             names=rep(c("", expression(Sync[1]), "", expression(Sync[1])), 1), yaxt="n", cex.axis=0.9)
axis(2, at=seq(2.5,7,by=0.5), las=1)
text(1.5, 4.5, "ns", cex=1.1)
text(4.25, 4.5, "ns", cex=1.1)
text(0.25, 7, "C", font=2, cex=1.4, pos=4, xpd=TRUE)
text(1, par("usr")[3]-0.55, "Exotic\nearly", xpd=T, adj=0.5, cex=0.9)
text(3.75, par("usr")[3]-0.55, "Exotic\nearly", xpd=T, adj=0.5, cex=0.9)
text(1.5, par("usr")[3]-1.2, "Grasses", font=2, xpd=T)
text(4.25, par("usr")[3]-1.2, "Grasses-Legumes", font=2, xpd=T)

par(bty="l")
par(mar=c(2,4.5,2,1)+0.1)
gf<-beanplot(BSenecio~Arrival*Community, data2, log="", what=c(1,0,1,1), border = NA, col = colors, ylab="", xlab="", las=1, at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(0,70), names=c("","","",""), yaxt="n")
axis(2, at=seq(0,70,by=10), las=1)
text(1.5, 70, "***", cex=2)
text(4.25, 70, "***", cex=2)
text(0.25, 70, "B", font=2, cex=1.4, pos=4, xpd=TRUE)

par(bty="l")
par(mar=c(4,4.5,0,1)+0.1)
gf<-beanplot(NSenecio~Arrival*Community, data2, log="", what=c(1,0,1,1), border = NA, col = colors, ylab="", xlab="", 
             las=1, at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(2.5,7), 
             names=rep(c("", expression(Sync[2]), "", expression(Sync[2])), 1), 
             yaxt="n", cex.axis=0.9)
axis(2, at=seq(2.5,7,by=0.5), las=1)
text(1.5, 7, "**", cex=2)
text(4.25, 7, "*", cex=2)
text(0.25, 7, "D", font=2, cex=1.4, pos=4, xpd=TRUE)
text(1, par("usr")[3]-0.55, "Exotic\nlate", xpd=T, adj=0.5, cex=0.9)
text(3.75, par("usr")[3]-0.55, "Exotic\nlate", xpd=T, adj=0.5, cex=0.9)
text(1.5, par("usr")[3]-1.2, "Grasses", font=2, xpd=T)
text(4.25, par("usr")[3]-1.2, "Grasses-Legumes", font=2, xpd=T)

dev.off()

#Plot effect sizes
es.senecio<-rbind(es.senecio1, es.senecio2)

layout(1)
par(mar=c(5,10,4,2)+0.1)
gf<-barplot(es.senecio$Hedges, horiz=TRUE, xlim=c(-7,5), xlab="Hedges' g effect size ± 95% CI", border=NA, 
            col="white", xaxt="n", names.arg = rownames(es.senecio), las=1, ylab="")
axis(1, at=seq(-7,5, by=1), las=1)
axis(2, at=gf, labels=FALSE)
abline(v=0, lty=2)
box(bty="l")
arrows(es.senecio$CIinf, gf, es.senecio$CIsup, gf, code=3, angle=90, length=0.05)
points(es.senecio$Hedges, gf, cex=1.4, pch=21, col="black", bg=c(rep("black", 4), rep("white", 4)))

#################################################
#################################################
#Analysis GRASSES: Late arrival - Synchronous 1
#################################################
#################################################

data1<-data[data$Arrival=="Late"|data$Arrival=="Synchronous1",]
data1$Arrival<-factor(data1$Arrival)
data1$Community<-factor(data1$Community)

#######################
#Stats biomass grasses
#######################

dotchart(data1$BGrasses) #No outlier detected

#Fit linear model with all observations
mod1<-lm(BGrasses~Community*Arrival, data1)
summary(mod1)
par(mfrow=c(3,2))
plot(mod1)
plot(cooks.distance(mod1), type="h")

#Fit GLM
mod2<-glm(BGrasses~Community*Arrival, data=data1, family=Gamma(link="log"))
mod3<-glm(BGrasses~Community*Arrival, data=data1, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3)

E2<-resid(mod2, type="pearson") #Pearson residuals
F2<-fitted(mod2) #Fitted values (same scale as the response variable)
eta<-predict(mod2, type="link") #Fitted values (same scale as the predictor function: log scale)
par(mfrow=c(3,2))
plot(F2, E2, xlab="Fitted values", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(eta, E2, xlab="Eta", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
boxplot(E2~Arrival, data=data1, ylab="Pearson residuals", xlab="Arrival")
abline(h=0, lty=2)
boxplot(E2~Community, data=data1, ylab="Pearson residuals", xlab="Community")
abline(h=0, lty=2)
plot(cooks.distance(mod2), type="h") #All values lower than 1, so no influential observation

summary(mod2)
drop1(mod2, test="F") #No significant interaction between Arrival and Community
drop1(glm(BGrasses~Community+Arrival, data=data1, family=Gamma(link="log")), test="F")
par(mfrow=c(2,2))
plot(mod2) #By default, plot(glm object) uses the deviance residuals (and not the Pearson residuals)
DevianceExplained2<-100*(mod2$null.deviance-mod2$deviance)/mod2$null.deviance #% of the variation explained by the model

layout(1)
interaction.plot(x.factor=data1$Community, trace.factor=data1$Arrival, response=data1$BGrasses)

tmp<-expand.grid(Arrival=levels(data1$Arrival), Community=levels(data1$Community)) #We first compute the means of the response for all 8 combinations of the levels of Arrival and Community
X<-model.matrix(~Community*Arrival, data=tmp)
glht(mod2, linfct=X)
predict(mod2, newdata=tmp, type="link")
#We now construct a contrast matrix based on Tukey-contrasts for Arrival in a block-diagonal way (for each level of Community)
Tukey<-contrMat(table(data1$Arrival), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(data1$Community)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(data1$Community)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))
K<-rbind(K, matrix(0, nrow=2, ncol=4))
rownames(K)<-c(rownames(K)[1:2], "G:Late - GL:Late", "G:Synchronous1 - GL:Synchronous1")
K[3,]<-c(1,0,-1,0)
K[4,]<-c(0,1,0,-1)

summary(glht(mod2, linfct=K%*%X))
summary(glht(mod3, linfct=K%*%X))

########################
#Stats N content grasses
########################

dotchart(data1$NGrasses) #No outlier detected

mod1<-lm(NGrasses~Community*Arrival, data1)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
anova(mod1) #Significant effect of Arrival

mod2<-glm(NGrasses~Community*Arrival, data=data1, family=Gamma(link="log"))
mod3<-glm(NGrasses~Community*Arrival, data=data1, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3) #GLM does not seem to be a better model

tmp<-expand.grid(Arrival=levels(data1$Arrival), Community=levels(data1$Community)) #We first compute the means of the response for all 8 combinations of the levels of Arrival and Community
X<-model.matrix(~Community*Arrival, data=tmp)
glht(mod1, linfct=X)
predict(mod1, newdata=tmp, type="response")
#We now construct a contrast matrix based on Tukey-contrasts for Arrival in a block-diagonal way (for each level of Community)
Tukey<-contrMat(table(data1$Arrival), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(data1$Community)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(data1$Community)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))
K<-rbind(K, matrix(0, nrow=2, ncol=4))
rownames(K)<-c(rownames(K)[1:2], "G:Late - GL:Late", "G:Synchronous1 - GL:Synchronous1")
K[3,]<-c(1,0,-1,0)
K[4,]<-c(0,1,0,-1)

summary(glht(mod1, linfct=K%*%X))

########################
#Stats C content grasses
########################

dotchart(data1$CGrasses) #No outlier detected

mod1<-lm(CGrasses~Community*Arrival, data1)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
anova(mod1)

mod2<-glm(CGrasses~Community*Arrival, data=data1, family=Gamma(link="log"))
mod3<-glm(CGrasses~Community*Arrival, data=data1, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3) #GLM does not seem to be a better model

#################################################
#################################################
#Analysis GRASSES: Early arrival - Synchronous 2
#################################################
#################################################

data2<-data[data$Arrival=="Early"|data$Arrival=="Synchronous2",]
data2$Arrival<-factor(data2$Arrival)
data2$Community<-factor(data2$Community)

#######################
#Stats biomass grasses
#######################

dotchart(data2$BGrasses) #No visible outlier

#Fit linear model with all observations
mod1a<-lm(BGrasses~Community*Arrival, data2)
summary(mod1a)
par(mfrow=c(3,2))
plot(mod1a) #We clearly have non-homogeneity of variances between groups
plot(cooks.distance(mod1a), type="h")

mod2<-glm(BGrasses~Community*Arrival, data=data2, family=Gamma(link="log"))
mod3<-glm(BGrasses~Community*Arrival, data=data2, family=Gamma(link="identity"))

AIC(mod1a, mod2, mod3)

E2<-resid(mod2, type="pearson") #Pearson residuals
F2<-fitted(mod2) #Fitted values (same scale as the response variable)
eta<-predict(mod2, type="link") #Fitted values (same scale as the predictor function: log scale)
par(mfrow=c(3,2))
plot(F2, E2, xlab="Fitted values", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(eta, E2, xlab="Eta", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
boxplot(E2~Arrival, data=data2, ylab="Pearson residuals", xlab="Arrival")
abline(h=0, lty=2)
boxplot(E2~Community, data=data2, ylab="Pearson residuals", xlab="Community")
abline(h=0, lty=2)
plot(cooks.distance(mod2), type="h") #All values lower than 1, so no influential observation

summary(mod2)
drop1(mod2, test="F") #There is a significant interaction between Arrival and Community
par(mfrow=c(2,2))
plot(mod2) #By default, plot(glm object) uses the deviance residuals (and not the Pearson residuals)
DevianceExplained2<-100*(mod2$null.deviance-mod2$deviance)/mod2$null.deviance #% of the variation explained by the model

tmp<-expand.grid(Arrival=levels(data2$Arrival), Community=levels(data2$Community)) #We first compute the means of the response for all 8 combinations of the levels of Arrival and Community
X<-model.matrix(~Community*Arrival, data=tmp)
glht(mod2, linfct=X)
predict(mod2, newdata=tmp, type="link")
#We now construct a contrast matrix based on Tukey-contrasts for Arrival in a block-diagonal way (for each level of Community)
Tukey<-contrMat(table(data2$Arrival), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(data2$Community)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(data2$Community)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))
K<-rbind(K, matrix(0, nrow=2, ncol=4))
rownames(K)<-c(rownames(K)[1:2], "G:Early - GL:Early", "G:Synchronous2 - GL:Synchronous2")
K[3,]<-c(1,0,-1,0)
K[4,]<-c(0,1,0,-1)

summary(glht(mod2, linfct=K%*%X))

########################
#Stats N content grasses
########################

dotchart(data2$NGrasses) #No outlier detected

mod1<-lm(NGrasses~Community*Arrival, data2)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
anova(mod1)

mod2<-glm(NGrasses~Community*Arrival, data=data2, family=Gamma(link="log"))
mod3<-glm(NGrasses~Community*Arrival, data=data2, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3) #GLM does not seem to be a better model

tmp<-expand.grid(Arrival=levels(data2$Arrival), Community=levels(data2$Community)) #We first compute the means of the response for all 8 combinations of the levels of Arrival and Community
X<-model.matrix(~Community*Arrival, data=tmp)
glht(mod1, linfct=X)
predict(mod1, newdata=tmp, type="response")
#We now construct a contrast matrix based on Tukey-contrasts for Arrival in a block-diagonal way (for each level of Community)
Tukey<-contrMat(table(data2$Arrival), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(data2$Community)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(data2$Community)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))
K<-rbind(K, matrix(0, nrow=2, ncol=4))
rownames(K)<-c(rownames(K)[1:2], "G:Early - GL:Early", "G:Synchronous2 - GL:Synchronous2")
K[3,]<-c(1,0,-1,0)
K[4,]<-c(0,1,0,-1)

summary(glht(mod1, linfct=K%*%X))

########################
#Stats C content grasses
########################

dotchart(data2$CGrasses)

mod1<-lm(CGrasses~Community*Arrival, data2)
summary(mod1)
par(mfrow=c(3,2))
plot(mod1)
anova(mod1)
plot(cooks.distance(mod1), type="h")

mod2<-glm(CGrasses~Community*Arrival, data=data2, family=Gamma(link="log"))
mod3<-glm(CGrasses~Community*Arrival, data=data2, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3) #GLM does not seem to be a better model

#######################
#Nitrogen facilitation
#######################

model<-lm(deltaNGrasses~Community*Arrival, data)
plot(fitted(model), resid(model)); abline(h=0)
anova(model) #delta 15N not different between treatments

#################
#################
#Export Figure 3
#################
#################

colors<-c("white", "grey50", "black", "black")

#img<-readJPEG("grassescrop.jpg")
#imgheight<-968
#imgwidth<-1382
#imgratio<-imgheight/imgwidth

tiff(filename="Figure3.tif", res=1000, width=13.5, height=14, units="cm", compression="lzw", 
     pointsize=12)
layout(matrix(c(1,2,3,1,4,5), ncol=2, nrow=3), heights=c(1.7,4,4))

par(mar=c(0,0,0,0))
plot(1:10, ty="n", ylab="", xlab="", bty="n", xaxt="n", yaxt="n", asp=1)
#rasterImage(img,1,1,9/imgratio,10)

par(bty="l")
par(mar=c(3,4.5,1,1)+0.1)
gf<-beanplot(BGrasses~Arrival*Community, data1, log="", what=c(1,0,1,1), border = NA, col = colors, 
             ylab=expression(paste("Grass shoot dry weight (g ", m^-2, ")", sep="")), xlab="", las=1, 
             at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(0,700), names=c("","","",""), yaxt="n", 
             cex.lab=1)
axis(2, at=seq(0,700,by=100), las=1, cex.axis=1)
text(1.5, 700, "***", cex=2)
text(4.25, 300, "***", cex=2)
text(0.25, 700, "A", font=2, cex=1.4, pos=4, xpd=TRUE)

par(bty="l")
par(mar=c(4,4.5,0,1)+0.1)
gf<-beanplot(NGrasses~Arrival*Community, data1, log="", what=c(1,0,1,1), border = NA, col = colors, 
             ylab="Grass shoot N content (%)", xlab="", las=1, at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(2.5,6), 
             names=rep(c("", expression(Sync[1]), "", expression(Sync[1])), 1), yaxt="n", cex.axis=0.9)
axis(2, at=seq(2.5,6,by=0.5), las=1, cex.axis=1)
text(1.5, 4.7, "ns", cex=1.2)
text(4.25, 4.7, "ns", cex=1.2)
text(0.25, 6, "C", font=2, cex=1.4, pos=4, xpd=TRUE)
text((3.25+2.50)/2, 6, expression(paste(P[Arrival], "=0.038", sep="")), cex=1, xpd=TRUE)
text(1, par("usr")[3]-0.45, "Exotic\nlate", xpd=T, cex=0.9)
text(3.75, par("usr")[3]-0.45, "Exotic\nlate", xpd=T, cex=0.9)
text(1.5, par("usr")[3]-0.95, "Grasses", font=2, xpd=T, cex=1)
text(4.25, par("usr")[3]-0.95, "Grasses-Legumes", font=2, xpd=T, cex=1)

par(bty="l")
par(mar=c(3,4.5,1,1)+0.1)
gf<-beanplot(BGrasses~Arrival*Community, data2, log="", what=c(1,0,1,1), border = NA, col = colors, ylab="", xlab="", las=1, at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(0,170), names=c("","","",""), yaxt="n", cex.axis=1.1, cex.lab=1.1)
axis(2, at=seq(0,160,by=20), las=1, cex.axis=1)
text(1.5, 170, "***", cex=2)
text(4.25, 90, "***", cex=2)
text(0.25, 170, "B", font=2, cex=1.4, pos=4, xpd=TRUE)

par(bty="l")
par(mar=c(4,4.5,0,1)+0.1)
gf<-beanplot(NGrasses~Arrival*Community, data2, log="", what=c(1,0,1,1), border = NA, col = colors, ylab="", xlab="", 
             las=1, at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(2.5,6), 
             names=rep(c("", expression(Sync[2]), "", expression(Sync[2])), 1), yaxt="n", cex.axis=0.9)
axis(2, at=seq(2.5,6,by=0.5), las=1, cex.axis=1)
text(1.5, 6, "***", cex=2)
text(4.25, 6, "ns", cex=1.2)
lines(c(1.3, 1.3), c(gf$stats[1], 2.75), lty=2)
lines(c(3.45, 3.45), c(gf$stats[3], 2.75), lty=2)
lines(c(1.3, 3.45), c(2.75, 2.75), lty=2)
text((1.3+3.45)/2, 2.75, "*", pos=3, cex=2)
text(0.25, 6, "D", font=2, cex=1.4, pos=4, xpd=TRUE)
text(1, par("usr")[3]-0.45, "Exotic\nearly", xpd=T, cex=0.9)
text(3.75, par("usr")[3]-0.45, "Exotic\nearly", xpd=T, cex=0.9)
text(1.5, par("usr")[3]-0.95, "Grasses", font=2, xpd=T, cex=1)
text(4.25, par("usr")[3]-0.95, "Grasses-Legumes", font=2, xpd=T, cex=1)

dev.off()

#################################################
#################################################
#Analysis LEGUMES: Late arrival - Synchronous 1
#################################################
#################################################

data1<-data[data$Arrival=="Late"|data$Arrival=="Synchronous1",]
data1<-data1[data1$Community=="GL",]
data1$Arrival<-factor(data1$Arrival)
data1$Community<-factor(data1$Community)

#######################
#Stats biomass legumes
#######################

dotchart(data1$BLegumes) #No outlier detected

#Fit linear model with all observations
mod1<-lm(BLegumes~Arrival, data1)
summary(mod1)
par(mfrow=c(3,2))
plot(mod1)
plot(cooks.distance(mod1), type="h")
anova(mod1)

#Fit GLM
mod2<-glm(BLegumes~Arrival, data=data1, family=Gamma(link="log"))
mod3<-glm(BLegumes~Arrival, data=data1, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3)

E2<-resid(mod2, type="pearson") #Pearson residuals
F2<-fitted(mod2) #Fitted values (same scale as the response variable)
eta<-predict(mod2, type="link") #Fitted values (same scale as the predictor function: log scale)
par(mfrow=c(3,2))
plot(F2, E2, xlab="Fitted values", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(eta, E2, xlab="Eta", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
boxplot(E2~Arrival, data=data1, ylab="Pearson residuals", xlab="Arrival")
abline(h=0, lty=2)
boxplot(E2~Community, data=data1, ylab="Pearson residuals", xlab="Community")
abline(h=0, lty=2)
plot(cooks.distance(mod2), type="h") #All values lower than 1, so no influential observation

summary(mod2)
drop1(mod2, test="F") #Significant effect of Arrival
par(mfrow=c(2,2))
plot(mod2) #By default, plot(glm object) uses the deviance residuals (and not the Pearson residuals)
DevianceExplained2<-100*(mod2$null.deviance-mod2$deviance)/mod2$null.deviance #% of the variation explained by the model

########################
#Stats N content legumes
########################

dotchart(data1$NLegumes) #No outlier detected

mod1<-lm(NLegumes~Arrival, data1)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
anova(mod1)

mod2<-glm(NLegumes~Arrival, data=data1, family=Gamma(link="log"))
mod3<-glm(NLegumes~Arrival, data=data1, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3) #GLM does not seem to be a better model

drop1(mod2, test="F")

########################
#Stats C content legumes
########################

dotchart(data1$CLegumes) #No outlier detected

mod1<-lm(CLegumes~Arrival, data1)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
anova(mod1)

mod2<-glm(CLegumes~Arrival, data=data1, family=Gamma(link="log"))
mod3<-glm(CLegumes~Arrival, data=data1, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3) #GLM does not seem to be a better model

#################################################
#################################################
#Analysis LEGUMES: Early arrival - Synchronous 2
#################################################
#################################################

data2<-data[data$Arrival=="Early"|data$Arrival=="Synchronous2",]
data2<-data2[data2$Community=="GL",]
data2$Arrival<-factor(data2$Arrival)
data2$Community<-factor(data2$Community)

#######################
#Stats biomass legumes
#######################

dotchart(data2$BLegumes) #No outlier detected

#Fit linear model with all observations
mod1<-lm(BLegumes~Arrival, data2)
summary(mod1)
par(mfrow=c(3,2))
plot(mod1)
plot(cooks.distance(mod1), type="h")
anova(mod1)

#Fit GLM
mod2<-glm(BLegumes~Arrival, data=data2, family=Gamma(link="log"))
mod3<-glm(BLegumes~Arrival, data=data2, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3)

E2<-resid(mod2, type="pearson") #Pearson residuals
F2<-fitted(mod2) #Fitted values (same scale as the response variable)
eta<-predict(mod2, type="link") #Fitted values (same scale as the predictor function: log scale)
par(mfrow=c(3,2))
plot(F2, E2, xlab="Fitted values", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
plot(eta, E2, xlab="Eta", ylab="Pearson residuals")
abline(v=0, h=0, lty=2)
boxplot(E2~Arrival, data=data2, ylab="Pearson residuals", xlab="Arrival")
abline(h=0, lty=2)
boxplot(E2~Community, data=data2, ylab="Pearson residuals", xlab="Community")
abline(h=0, lty=2)
plot(cooks.distance(mod2), type="h") #All values lower than 1, so no influential observation

summary(mod2)
drop1(mod2, test="F") #Significant effect of Arrival
par(mfrow=c(2,2))
plot(mod2) #By default, plot(glm object) uses the deviance residuals (and not the Pearson residuals)
DevianceExplained2<-100*(mod2$null.deviance-mod2$deviance)/mod2$null.deviance #% of the variation explained by the model

########################
#Stats N content legumes
########################

dotchart(data2$NLegumes) #No outlier detected

mod1<-lm(NLegumes~Arrival, data2)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
anova(mod1)

mod2<-glm(NLegumes~Arrival, data=data2, family=Gamma(link="log"))
mod3<-glm(NLegumes~Arrival, data=data2, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3) #GLM does not seem to be a better model

########################
#Stats C content legumes
########################

dotchart(data2$CLegumes) #No outlier detected

mod1<-lm(CLegumes~Arrival, data2)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
anova(mod1)

mod2<-glm(CLegumes~Arrival, data=data2, family=Gamma(link="log"))
mod3<-glm(CLegumes~Arrival, data=data2, family=Gamma(link="identity"))

AIC(mod1, mod2, mod3) #GLM does not seem to be a better model

#######################
#Nitrogen facilitation
#######################

newdata<-data[which(data$Community=="GL"),]
model<-lm(deltaNLegumes~Arrival, newdata)
plot(fitted(model), resid(model)); abline(h=0)
anova(model) #delta 15N not different between treatments

#################
#################
#Export Figure 4
#################
#################

colors<-c("white", "grey50", "black", "black")

#img<-readJPEG("legumescrop.jpg")
#imgheight<-1022
#imgwidth<-1276
#imgratio<-imgheight/imgwidth

tiff(filename="Figure4.tif", res=1000, width=9, height=11, units="cm", compression="lzw", pointsize=10)
layout(matrix(c(1,2,3,1,4,5), ncol=2, nrow=3), heights=c(1.7,4,4))

par(mar=c(0,0,0,0))
plot(1:10, ty="n", ylab="", xlab="", bty="n", xaxt="n", yaxt="n", asp=1)
#rasterImage(img,1,1,9/imgratio,10)

par(bty="l")
par(mar=c(2,4.5,2,0.5)+0.1)
gf<-beanplot(BLegumes~Arrival, data1, log="", what=c(1,0,1,1), border = NA, col = colors, 
             ylab=expression(paste("Legume shoot dry weight (g ", m^-2, ")", sep="")), xlab="", las=1, ylim=c(0,400), 
             names=c("",""), yaxt="n", cex.axis=1, cex.lab=1, at=c(1,2.5), xlim=c(0.25,3.25), ll=0.11)
axis(2, at=seq(0,400,by=50), las=1, cex.axis=1)
text(1.75, 400, "***", cex=2)
text(0.1, 400, "A", font=2, cex=1.3, pos=4, xpd=TRUE)

par(bty="l")
par(mar=c(3,4.5,1,0.5)+0.1)
gf<-beanplot(NLegumes~Arrival, data1, log="", what=c(1,0,1,1), border = NA, col = colors, ylab="Legume shoot N content (%)", 
             xlab="", las=1, ylim=c(2,5), names=c("Exotic late", expression(Sync[1])), yaxt="n", cex.axis=1, cex.lab=1,
             at=c(1,2.5), xlim=c(0.25,3.25), ll=0.11)
axis(2, at=seq(2,5,by=0.5), las=1, cex.axis=1)
text(1.75, 4, "ns", cex=1.1)
text(0.1, 5, "C", font=2, cex=1.3, pos=4, xpd=TRUE)

par(bty="l")
par(mar=c(2,4,2,1)+0.1)
gf<-beanplot(BLegumes~Arrival, data2, log="", what=c(1,0,1,1), border = NA, col = colors, ylab="", xlab="", 
             las=1, ylim=c(0,80), names=c("",""), yaxt="n", cex.axis=1, cex.lab=1,
             at=c(1,2.5), xlim=c(0.25,3.25), ll=0.11)
axis(2, at=seq(0,80,by=10), las=1, cex.axis=1)
text(1.75, 80, "***", cex=2)
text(0.1, 80, "B", font=2, cex=1.3, pos=4, xpd=TRUE)

par(bty="l")
par(mar=c(3,4,1,1)+0.1)
gf<-beanplot(NLegumes~Arrival, data2, log="", what=c(1,0,1,1), border = NA, col = colors, ylab="", xlab="", 
             las=1, ylim=c(2,5), names=c("Exotic early", expression(Sync[2])), yaxt="n", cex.axis=1, cex.lab=1,
             at=c(1,2.5), xlim=c(0.25,3.25), ll=0.11)
axis(2, at=seq(2,5,by=0.5), las=1, cex.axis=1)
text(1.75, 5, "***", cex=2)
text(0.1, 5, "D", font=2, cex=1.3, pos=4, xpd=TRUE)

dev.off()

#################
#################
#Export Figure 5
#################
#################

data$BNatives<-data$BGrasses+data$BLegumes

###########################
#Benefit of arriving early
###########################

#For Senecio

data1<-data[data$Arrival=="Early"|data$Arrival=="Synchronous1",]
data1$Arrival<-factor(data1$Arrival)
data1$Community<-factor(data1$Community)

B1<-data.frame(Community=c(rep("G", 5), rep("GL", 5)), Replicate=c(c(1:5), c(1:5)), B=rep(NA, 10))

for (i in 1:10){
  Bij<-data1$BSenecio[data1$Arrival=="Early" & data1$Community==B1$Community[i] & data1$Replicate==B1$Replicate[i]]
  Bsync1<-data1$BSenecio[data1$Arrival=="Synchronous1" & data1$Community==B1$Community[i] & data1$Replicate==B1$Replicate[i]]
  B1$B[i]<-(Bij-Bsync1)/(Bij+Bsync1)}

confB1G<-CIbootstrap(x=B1$B[1:5])
confB1GL<-CIbootstrap(x=B1$B[6:10])

#For all natives

data1<-data[data$Arrival=="Late"|data$Arrival=="Synchronous1",]
data1$Arrival<-factor(data1$Arrival)
data1$Community<-factor(data1$Community)

B3<-data.frame(Community=c(rep("G", 5), rep("GL", 5)), Replicate=c(c(1:5), c(1:5)), B=rep(NA, 10))

for (i in 1:10){
  Bji<-data1$BNatives[data1$Arrival=="Late" & data1$Community==B3$Community[i] & data1$Replicate==B3$Replicate[i]]
  Bsync1<-data1$BNatives[data1$Arrival=="Synchronous1" & data1$Community==B3$Community[i] & data1$Replicate==B3$Replicate[i]]
  B3$B[i]<-(Bji-Bsync1)/(Bji+Bsync1)}

confB3G<-CIbootstrap(x=B3$B[1:5])
confB3GL<-CIbootstrap(x=B3$B[6:10])

#Make one dataset

B<-rbind(B1,B3)
B$Species<-factor(c(rep("aSenecio", 10), rep("bNatives", 10)))

modelB<-lm(B~Species*Community, B)
#plot(fitted(modelB), resid(modelB), pch=1)
#abline(h=0)
anova(modelB)

tmp<-expand.grid(Community=levels(B$Community), Species=levels(B$Species))
X<-model.matrix(~Species*Community, data=tmp)
glht(modelB, linfct=X)
predict(modelB, newdata=tmp, type="response")
Tukey<-contrMat(table(B$Community), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(B$Species)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(B$Species)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))

summary(glht(modelB, linfct=K%*%X))

#Plot results for shoot dry weight

tiff(filename="Figure5.tif", res=1000, compression="lzw", width=7, height=12, units="cm", 
     pointsize=7)
layout(matrix(1:2, ncol=1, nrow=2))

par(bty="l", mar=c(0.5,4,4,2)+0.1)
gf<-beanplot(B~Community*Species, B, side="both", log="", what=c(0,0,1,1), border = NA, 
             col = list(c("white", "black", "black", "black"),
                        c("white", "gray50", "gray50", "gray50")), 
             ylab=expression(paste("Benefit of arriving early (", italic(B[E]), " or ", italic(B[N]), ")", sep="")), las=1, ylim=c(-1,1), xlab="", 
             cex.main=1, names=c("", ""), at=c(1,2.5), xlim=c(0.5,3),
             yaxt="n", ll=0.11, cex.main=1.3)
abline(h=0, lty=3, lwd=0.8)
axis(2, at=seq(-1,1, by=0.5), las=1)
arrows(0.75, confB1G[1], 0.75, confB1G[2], code=3, angle=90, length=0.02, col="black")
arrows(1.25, confB1GL[1], 1.25, confB1GL[2], code=3, angle=90, length=0.02, col="gray50")
arrows(2.25, confB3G[1], 2.25, confB3G[2], code=3, angle=90, length=0.02, col="black")
arrows(2.75, confB3GL[1], 2.75, confB3GL[2], code=3, angle=90, length=0.02, col="gray50")
text(0.55, 0.95, "A", font=2, cex=1.4)
text(1, 0.7, "*", cex=2)
text(2.5, 0.6, "ns", cex=1.1)

#######################
#Cost of arriving late
#######################

#For Senecio

data2<-data[data$Arrival=="Late"|data$Arrival=="Synchronous2",]
data2$Arrival<-factor(data2$Arrival)
data2$Community<-factor(data2$Community)

P1<-data.frame(Community=c(rep("G", 5), rep("GL", 5)), Replicate=c(c(1:5), c(1:5)), P=rep(NA, 10))

for (i in 1:10){
  Bji<-data2$BSenecio[data2$Arrival=="Late" & data2$Community==P1$Community[i] & data2$Replicate==P1$Replicate[i]]
  Bsync2<-data2$BSenecio[data2$Arrival=="Synchronous2" & data2$Community==P1$Community[i] & data2$Replicate==P1$Replicate[i]]
  P1$P[i]<-(Bji-Bsync2)/(Bji+Bsync2)}

confP1G<-CIbootstrap(x=P1$P[1:5])
confP1GL<-CIbootstrap(x=P1$P[6:10])

#For natives

data2<-data[data$Arrival=="Early"|data$Arrival=="Synchronous2",]
data2$Arrival<-factor(data2$Arrival)
data2$Community<-factor(data2$Community)

P3<-data.frame(Community=c(rep("G", 5), rep("GL", 5)), Replicate=c(c(1:5), c(1:5)), P=rep(NA, 10))

for (i in 1:10){
  Bij<-data2$BNatives[data2$Arrival=="Early" & data2$Community==P3$Community[i] & data2$Replicate==P3$Replicate[i]]
  Bsync2<-data2$BNatives[data2$Arrival=="Synchronous2" & data2$Community==P3$Community[i] & data2$Replicate==P3$Replicate[i]]
  P3$P[i]<-(Bij-Bsync2)/(Bij+Bsync2)}

confP3G<-CIbootstrap(x=P3$P[1:5])
confP3GL<-CIbootstrap(x=P3$P[6:10])

#Make one dataset

B<-rbind(P1,P3)
B$Species<-factor(c(rep("aSenecio", 10), rep("bNatives", 10)))

modelB<-lm(P~Species*Community, B)
#plot(fitted(modelB), resid(modelB), pch=1)
#abline(h=0)
anova(modelB)

tmp<-expand.grid(Community=levels(B$Community), Species=levels(B$Species))
X<-model.matrix(~Species*Community, data=tmp)
glht(modelB, linfct=X)
predict(modelB, newdata=tmp, type="response")
Tukey<-contrMat(table(B$Community), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(B$Species)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(B$Species)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))

summary(glht(modelB, linfct=K%*%X))

#Plot results for shoot dry weight

par(bty="l", mar=c(3,4,1.5,2)+0.1)
gf<-beanplot(P~Community*Species, B, side="both", log="", what=c(0,0,1,1), border = NA, 
             col = list(c("white", "black", "black", "black"),
                        c("white", "gray50", "gray50", "gray50")), 
             ylab=expression(paste("Cost of arriving late (", italic(P[E]), " or ", italic(P[N]), ")", sep="")), las=1, ylim=c(-1,1), xlab="", 
             cex.main=1, names=c(expression(italic("S. inaequidens")), "Natives"), at=c(1,2.5), xlim=c(0.5,3),
             yaxt="n", main="", ll=0.11, cex.main=1.3)
abline(h=0, lty=3, lwd=0.8)
axis(2, at=seq(-1,1, by=0.5), las=1)
arrows(0.75, confP1G[1], 0.75, confP1G[2], code=3, angle=90, length=0.02, col="black")
arrows(1.25, confP1GL[1], 1.25, confP1GL[2], code=3, angle=90, length=0.02, col="gray50")
arrows(2.25, confP3G[1], 2.25, confP3G[2], code=3, angle=90, length=0.02, col="black")
arrows(2.75, confP3GL[1], 2.75, confP3GL[2], code=3, angle=90, length=0.02, col="gray50")
text(0.55, 0.95, "B", font=2, cex=1.4)
text(1, -0.6, "ns", cex=1.1)
text(2.5, -0.05, "***", cex=2)
legend("topright", legend=c("Grasses", "Grasses-Legumes"), lty=1, lwd=2, col=c("black", "gray50"), bty="n", cex=0.9)

dev.off()

#################
#################
#Export Figure 6
#################
#################

tiff(filename="Figure6.tif", res=1000, compression="lzw", width=14, height=14, units="cm", pointsize=9)
par(mfrow=c(2,2))

#Plot all functional groups
data1<-data
data1<-data1[-which(is.na(data1$deltaNSenecio)==TRUE),]

par(mar=c(2,4.5,4,0)+0.1)
plot(data1$NSenecio, data1$deltaNSenecio, xlab="", ylab=expression(paste(delta^15, "N (","\211",")", sep="")), bty="l", type="n", xlim=c(2,7), ylim=c(0,8), las=1)
points(data1$NSenecio, data1$deltaNSenecio, pch=1, col="black", cex=1)
points(data$NGrasses, data$deltaNGrasses, pch=4, col="black", cex=1)
points(data$NLegumes, data$deltaNLegumes, pch=16, col="black", cex=1)
segments(x0=min(data1$NSenecio, na.rm=TRUE), y0=mean(data1$deltaNSenecio, na.rm=TRUE), x1=max(data1$NSenecio, na.rm=TRUE), 
         y1=mean(data1$deltaNSenecio, na.rm=TRUE), lty=3, col="black")
segments(x0=min(data$NGrasses, na.rm=TRUE), y0=mean(data$deltaNGrasses, na.rm=TRUE), x1=max(data$NGrasses, na.rm=TRUE), 
         y1=mean(data$deltaNGrasses, na.rm=TRUE), lty=2, col="black")
segments(x0=min(data$NLegumes, na.rm=TRUE), y0=mean(data$deltaNLegumes, na.rm=TRUE), x1=max(data$NLegumes, na.rm=TRUE), 
         y1=mean(data$deltaNLegumes, na.rm=TRUE), lty=1, col="black")
axis(2, at=c(1,3,5,7), las=1)
legend("bottomright", legend=c(expression(italic("S. inaequidens")), "Grasses", "Legumes"), pch=c(1,4,16), lty=c(3,2,1),
       col=c("black", "black", "black"), bty="n", pt.cex=1.2, cex=0.8)
text(2.1,8, "A", cex=1.5, font=2)

#Plot per functional group and per treatment
data1$pch<-NA
data1$col<-NA
data$pch<-NA
data$col<-NA

pal<-brewer.pal(name="Paired", 8)

data1$col[which(data1$Arrival=="Synchronous1")]<-adjustcolor(pal[2], alpha.f=0.8)
data1$col[which(data1$Arrival=="Synchronous2")]<-adjustcolor(pal[4], alpha.f=0.8)
data1$col[which(data1$Arrival=="Early")]<-adjustcolor(pal[6], alpha.f=0.8)
data1$col[which(data1$Arrival=="Late")]<-adjustcolor(pal[8], alpha.f=0.8)
data1$pch[which(data1$Community=="G" & data1$Arrival=="Synchronous1")]<-1
data1$pch[which(data1$Community=="G" & data1$Arrival=="Synchronous2")]<-2
data1$pch[which(data1$Community=="G" & data1$Arrival=="Early")]<-5
data1$pch[which(data1$Community=="G" & data1$Arrival=="Late")]<-0
data1$pch[which(data1$Community=="GL" & data1$Arrival=="Synchronous1")]<-16
data1$pch[which(data1$Community=="GL" & data1$Arrival=="Synchronous2")]<-17
data1$pch[which(data1$Community=="GL" & data1$Arrival=="Early")]<-18
data1$pch[which(data1$Community=="GL" & data1$Arrival=="Late")]<-15
data$col[which(data$Arrival=="Synchronous1")]<-adjustcolor(pal[2], alpha.f=0.8)
data$col[which(data$Arrival=="Synchronous2")]<-adjustcolor(pal[4], alpha.f=0.8)
data$col[which(data$Arrival=="Early")]<-adjustcolor(pal[6], alpha.f=0.8)
data$col[which(data$Arrival=="Late")]<-adjustcolor(pal[8], alpha.f=0.8)
data$pch[which(data$Community=="G" & data$Arrival=="Synchronous1")]<-1
data$pch[which(data$Community=="G" & data$Arrival=="Synchronous2")]<-2
data$pch[which(data$Community=="G" & data$Arrival=="Early")]<-5
data$pch[which(data$Community=="G" & data$Arrival=="Late")]<-0
data$pch[which(data$Community=="GL" & data$Arrival=="Synchronous1")]<-16
data$pch[which(data$Community=="GL" & data$Arrival=="Synchronous2")]<-17
data$pch[which(data$Community=="GL" & data$Arrival=="Early")]<-18
data$pch[which(data$Community=="GL" & data$Arrival=="Late")]<-15

#Plot Senecio
#img<-readJPEG("senecio3crop.jpg")
#imgheight<-3648
#imgwidth<-1872
#imgratio<-imgheight/imgwidth

par(mar=c(2,2.5,4,2)+0.1)
plot(data1$NSenecio, data1$deltaNSenecio, xlab="", ylab="", type="n", xlim=c(2,7), ylim=c(0,8), las=1, bty="l")
points(data1$NSenecio, data1$deltaNSenecio, pch=data1$pch, col=data1$col, cex=1.2)
axis(2, at=c(1,3,5,7), las=1)
#xleft<-2
#ytop<-2.5
#ybottom<-0
#xright<-(ytop-ybottom+imgratio*xleft)/imgratio
#deltax<-xright-xleft
#xright<-xleft+deltax*(5/8)
#rasterImage(img, xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop)
text(2.1,8, "B", cex=1.5, font=2)
legend("bottomright", bty="n", legend=c("Synchronous1", "Exotic late", "Synchronous2", "Exotic early"), 
       pch=c(16,15,17,18), col=adjustcolor(pal[c(2,8,4,6)], alpha.f=0.8), pt.cex=1.2, cex=0.8)

#Plot grasses
#img<-readJPEG("grassescrop.jpg")
#imgheight<-968
#imgwidth<-1382
#imgratio<-imgheight/imgwidth

par(mar=c(5,4.5,1,0)+0.1)
plot(data$NGrasses, data$deltaNGrasses, xlab="N content (%)", ylab=expression(paste(delta^15, "N (","\211",")", sep="")), bty="l", type="n", xlim=c(2,7), ylim=c(0,8), las=1)
points(data$NGrasses, data$deltaNGrasses, pch=data$pch, col=data$col, cex=1.2)
axis(2, at=c(1,3,5,7), las=1)
#xright<-7
#ytop<-2
#ybottom<-0
#xleft<-(imgratio*xright-ytop+ybottom)/imgratio
#deltax<-xright-xleft
#xleft<-xright-deltax*(5/8)
#rasterImage(img, xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop)
text(2.1,8, "C", cex=1.5, font=2)

#Plot legumes
#img<-readJPEG("legumescrop.jpg")
#imgheight<-1022
#imgwidth<-1276
#imgratio<-imgheight/imgwidth

par(mar=c(5,2.5,1,2)+0.1)
plot(data$NLegumes, data$deltaNLegumes, xlab="N content (%)", ylab="", type="n", xlim=c(2,7), ylim=c(0,8), las=1, bty="l")
points(data$NLegumes, data$deltaNLegumes, pch=data$pch, col=data$col, cex=1.2)
axis(2, at=c(1,3,5,7), las=1)
#xright<-7
#ytop<-2
#ybottom<-0
#xleft<-(imgratio*xright-ytop+ybottom)/imgratio
#deltax<-xright-xleft
#xleft<-xright-deltax*(5/8)
#rasterImage(img, xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop)
text(2.1,8, "D", cex=1.5, font=2)

dev.off()

###########
###########
#Figure S1
###########
###########

#Senecio
data1<-data[data$Arrival=="Early"|data$Arrival=="Synchronous1",]
data1$Arrival<-factor(data1$Arrival)
data1$Community<-factor(data1$Community)

data2<-data[data$Arrival=="Late"|data$Arrival=="Synchronous2",]
data2$Arrival<-factor(data2$Arrival)
data2$Community<-factor(data2$Community)

colors<-c("white", "grey50", "black", "black")

#img<-readJPEG("senecio3crop.jpg")
#imgheight<-3648
#imgwidth<-1872
#imgratio<-imgheight/imgwidth

tiff(filename="FigureS1.tif", res=1000, width=15, height=7, units="cm", compression="lzw", pointsize=11)
layout(matrix(c(1,2,3), ncol=3, nrow=1), widths=c(3,6.5,6.5))

par(mar=c(5,0,4,0)+0.1)
plot(1:10, ty="n", ylab="", xlab="", bty="n", xaxt="n", yaxt="n", asp=1)
#rasterImage(img,1.7,-2.5,9.7,-2.5+8*imgratio)

par(bty="l")
par(mar=c(4,4,3,1)+0.1)
gf<-beanplot(CSenecio~Arrival*Community, data1, log="", what=c(1,0,1,1), border = NA, col = colors, 
             ylab=expression(paste(italic("Senecio"), " shoot C content (%)", sep="")), 
             xlab="", las=1, at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(25,45), 
             names=c("",expression(Sync[1]),"",expression(Sync[1])), yaxt="n")
axis(2, at=seq(25,45,by=2.5), las=1)
text(1.5, 42.5, "ns", cex=1.1)
text(4.25, 42.5, "*", cex=2)
text(0.25, 45, "A", font=2, cex=1.4, pos=4, xpd=TRUE)
text(1, par("usr")[3]-2.2, "Exotic\nearly", xpd=T, cex=1)
text(3.75, par("usr")[3]-2.2, "Exotic\nearly", xpd=T, cex=1)
text(1.5, par("usr")[3]-4.7, "Grasses", font=2, xpd=T, cex=1)
text(4.25, par("usr")[3]-4.7, "Grasses-Legumes", font=2, xpd=T, cex=1)

par(bty="l")
par(mar=c(4,3,3,2)+0.1)
gf<-beanplot(CSenecio~Arrival*Community, data2, log="", what=c(1,0,1,1), border = NA, col = colors, ylab="", xlab="", 
             las=1, at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(25,45), 
             names=c("",expression(Sync[2]),"",expression(Sync[2])), yaxt="n")
axis(2, at=seq(25,45,by=2.5), las=1)
text(1.5, 40, "ns", cex=1.1)
text(4.25, 38, "ns", cex=1.1)
text(0.25, 45, "B", font=2, cex=1.4, pos=4, xpd=TRUE)
text(1, par("usr")[3]-2.2, "Exotic\nlate", xpd=T, cex=1)
text(3.75, par("usr")[3]-2.2, "Exotic\nlate", xpd=T, cex=1)
text(1.5, par("usr")[3]-4.7, "Grasses", font=2, xpd=T, cex=1)
text(4.25, par("usr")[3]-4.7, "Grasses-Legumes", font=2, xpd=T, cex=1)

dev.off()

###########
###########
#Figure S2
###########
###########

#Grasses
data1<-data[data$Arrival=="Late"|data$Arrival=="Synchronous1",]
data1$Arrival<-factor(data1$Arrival)
data1$Community<-factor(data1$Community)

data2<-data[data$Arrival=="Early"|data$Arrival=="Synchronous2",]
data2$Arrival<-factor(data2$Arrival)
data2$Community<-factor(data2$Community)

colors<-c("white", "grey50", "black", "black")

#img<-readJPEG("grassescrop.jpg")
#imgheight<-968
#imgwidth<-1382
#imgratio<-imgheight/imgwidth

tiff(filename="FigureS2.tif", res=1000, width=14, height=10, units="cm", 
     compression="lzw", pointsize=10)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2), heights=c(1,3))

par(mar=c(0,0,0,0))
plot(1:10, ty="n", ylab="", xlab="", bty="n", xaxt="n", yaxt="n", asp=1)
#rasterImage(img,1,1,9/imgratio,10)

par(bty="l")
par(mar=c(4,4,3,1)+0.1)
gf<-beanplot(CGrasses~Arrival*Community, data1, log="", what=c(1,0,1,1), border = NA, col = colors, 
             ylab="Grass shoot C content (%)", 
             xlab="", las=1, at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(30,45), 
             names=c("",expression(Sync[1]),"",expression(Sync[1])), yaxt="n", cex.axis=1, cex.lab=1)
axis(2, at=seq(30,45,by=2.5), las=1, cex.axis=1)
text(1.5, 43.5, "ns", cex=1.1)
text(4.25, 43.5, "ns", cex=1.1)
text(0.25, 45, "A", font=2, cex=1.4, pos=4, xpd=TRUE)
text(1, par("usr")[3]-1.7, "Exotic\nlate", xpd=T, cex=1)
text(3.75, par("usr")[3]-1.7, "Exotic\nlate", xpd=T, cex=1)
text(1.5, par("usr")[3]-3.8, "Grasses", font=2, xpd=T, cex=1)
text(4.25, par("usr")[3]-3.8, "Grasses-Legumes", font=2, xpd=T, cex=1)

par(bty="l")
par(mar=c(4,3,3,2)+0.1)
gf<-beanplot(CGrasses~Arrival*Community, data2, log="", what=c(1,0,1,1), border = NA, col = colors, ylab="", 
             xlab="", las=1, at=c(1,2,3.75,4.75), xlim=c(0.5,5.25), ylim=c(30,45), 
             names=c("",expression(Sync[2]),"",expression(Sync[2])), yaxt="n", cex.axis=1, cex.lab=1)
axis(2, at=seq(30,45,by=2.5), las=1, cex.axis=1)
text(1.5, 43.5, "ns", cex=1.1)
text(4.25, 41.5, "ns", cex=1.1)
text(0.25, 45, "B", font=2, cex=1.4, pos=4, xpd=TRUE)
text(1, par("usr")[3]-1.7, "Exotic\nearly", xpd=T, cex=1)
text(3.75, par("usr")[3]-1.7, "Exotic\nearly", xpd=T, cex=1)
text(1.5, par("usr")[3]-3.8, "Grasses", font=2, xpd=T, cex=1)
text(4.25, par("usr")[3]-3.8, "Grasses-Legumes", font=2, xpd=T, cex=1)

dev.off()

###########
###########
#Figure S3
###########
###########

#Legumes
data1<-data[data$Arrival=="Late"|data$Arrival=="Synchronous1",]
data1<-data1[data1$Community=="GL",]
data1$Arrival<-factor(data1$Arrival)
data1$Community<-factor(data1$Community)

data2<-data[data$Arrival=="Early"|data$Arrival=="Synchronous2",]
data2<-data2[data2$Community=="GL",]
data2$Arrival<-factor(data2$Arrival)
data2$Community<-factor(data2$Community)

colors<-c("white", "grey50", "black", "black")

#img<-readJPEG("legumescrop.jpg")
#imgheight<-968
#imgwidth<-1382
#imgratio<-imgheight/imgwidth

tiff(filename="FigureS3.tif", res=1000, width=14, height=10, units="cm", 
     compression="lzw", pointsize=10)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2), heights=c(1,3))

par(mar=c(0,0,0,0))
plot(1:10, ty="n", ylab="", xlab="", bty="n", xaxt="n", yaxt="n", asp=1)
#rasterImage(img,1,1,9/imgratio,10)

par(bty="l")
par(mar=c(4,4,3,1)+0.1)
gf<-beanplot(CLegumes~Arrival*Community, data1, log="", what=c(1,0,1,1), border = NA, col = colors, 
             ylab="Legume shoot C content (%)", 
             xlab="", las=1, at=c(1,2.5), ylim=c(35,44), xlim=c(0.25,3.25), ll=0.11,
             names=c("Exotic late",expression(Sync[1])), yaxt="n", cex.axis=1, cex.lab=1)
axis(2, at=seq(35,44,by=1), las=1, cex.axis=1)
text(1.75, 43, "ns", cex=1.1)
text(0.15, 44, "A", font=2, cex=1.4, pos=4, xpd=TRUE)

par(bty="l")
par(mar=c(4,3,3,2)+0.1)
gf<-beanplot(CLegumes~Arrival*Community, data2, log="", what=c(1,0,1,1), border = NA, col = colors, ylab="", 
             xlab="", las=1, at=c(1,2.5), ylim=c(35,44), xlim=c(0.25,3.25), ll=0.11,
             names=c("Exotic early", expression(Sync[2])), yaxt="n", cex.axis=1, cex.lab=1)
axis(2, at=seq(35,44,by=1), las=1, cex.axis=1)
text(1.75, 42, "ns", cex=1.1)
text(0.15, 44, "B", font=2, cex=1.4, pos=4, xpd=TRUE)

dev.off()

###########
###########
#Figure S4
###########
###########

###########################
#Benefit of arriving early
###########################

data$NconcSenecio<-(data$BSenecio*(data$NSenecio/100))
data$NconcGrasses<-(data$BGrasses*(data$NGrasses/100))
data$NconcLegumes<-(data$BLegumes*(data$NLegumes/100))

data$NconcNatives<-NA
for (i in 1:nrow(data)){
  if(is.na(data$NLegumes[i])==TRUE){data$NconcNatives[i]<-(data$BGrasses[i]*(data$NGrasses[i]/100))}
  else {data$NconcNatives[i]<-(data$BGrasses[i]*(data$NGrasses[i]/100)+data$BLegumes[i]*(data$NLegumes[i]/100))}}

#For Senecio

data1<-data[data$Arrival=="Early"|data$Arrival=="Synchronous1",]
data1$Arrival<-factor(data1$Arrival)
data1$Community<-factor(data1$Community)

B1<-data.frame(Community=c(rep("G", 5), rep("GL", 5)), Replicate=c(c(1:5), c(1:5)), B=rep(NA, 10))

for (i in 1:10){
  Bij<-data1$NconcSenecio[data1$Arrival=="Early" & data1$Community==B1$Community[i] & data1$Replicate==B1$Replicate[i]]
  Bsync1<-data1$NconcSenecio[data1$Arrival=="Synchronous1" & data1$Community==B1$Community[i] & data1$Replicate==B1$Replicate[i]]
  B1$B[i]<-(Bij-Bsync1)/(Bij+Bsync1)}

confB1G<-CIbootstrap(x=B1$B[1:5])
confB1GL<-CIbootstrap(x=B1$B[6:10])

#For all natives

data1<-data[data$Arrival=="Late"|data$Arrival=="Synchronous1",]
data1$Arrival<-factor(data1$Arrival)
data1$Community<-factor(data1$Community)

B3<-data.frame(Community=c(rep("G", 5), rep("GL", 5)), Replicate=c(c(1:5), c(1:5)), B=rep(NA, 10))

for (i in 1:10){
  Bji<-data1$NconcNatives[data1$Arrival=="Late" & data1$Community==B3$Community[i] & data1$Replicate==B3$Replicate[i]]
  Bsync1<-data1$NconcNatives[data1$Arrival=="Synchronous1" & data1$Community==B3$Community[i] & data1$Replicate==B3$Replicate[i]]
  B3$B[i]<-(Bji-Bsync1)/(Bji+Bsync1)}

confB3G<-CIbootstrap(x=B3$B[1:5])
confB3GL<-CIbootstrap(x=B3$B[6:10])

#Make one dataset

B<-rbind(B1,B3)
B$Species<-factor(c(rep("aSenecio", 10), rep("bNatives", 10)))

modelB<-lm(B~Species*Community, B)
#plot(fitted(modelB), resid(modelB), pch=1)
#abline(h=0)
anova(modelB)

tmp<-expand.grid(Community=levels(B$Community), Species=levels(B$Species))
X<-model.matrix(~Species*Community, data=tmp)
glht(modelB, linfct=X)
predict(modelB, newdata=tmp, type="response")
Tukey<-contrMat(table(B$Community), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(B$Species)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(B$Species)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))

summary(glht(modelB, linfct=K%*%X))

#Plot results for shoot N content

tiff(filename="FigureS4.tif", res=1000, compression="lzw", width=7, height=12, units="cm", 
     pointsize=7)
layout(matrix(1:2, ncol=1, nrow=2))

par(bty="l", mar=c(0.5,4,4,2)+0.1)
gf<-beanplot(B~Community*Species, B, side="both", log="", what=c(0,0,1,1), border = NA, 
             col = list(c("white", "black", "black", "black"),
                        c("white", "gray50", "gray50", "gray50")), 
             ylab=expression(paste("Benefit of arriving early (", italic(B[E]), " or ", italic(B[N]), ")", sep="")), las=1, ylim=c(-1,1), xlab="", 
             cex.main=1, names=c("", ""), at=c(1,2.5), xlim=c(0.5,3),
             yaxt="n", ll=0.11, cex.main=1.3)
abline(h=0, lty=3, lwd=0.8)
axis(2, at=seq(-1,1, by=0.5), las=1)
arrows(0.75, confB1G[1], 0.75, confB1G[2], code=3, angle=90, length=0.02, col="black")
arrows(1.25, confB1GL[1], 1.25, confB1GL[2], code=3, angle=90, length=0.02, col="gray50")
arrows(2.25, confB3G[1], 2.25, confB3G[2], code=3, angle=90, length=0.02, col="black")
arrows(2.75, confB3GL[1], 2.75, confB3GL[2], code=3, angle=90, length=0.02, col="gray50")
text(0.55, 0.95, "A", font=2, cex=1.4)
text(1, 0.7, "*", cex=2)
text(2.5, 0.6, "ns", cex=1.1)

#######################
#Cost of arriving late
#######################

#For Senecio

data2<-data[data$Arrival=="Late"|data$Arrival=="Synchronous2",]
data2$Arrival<-factor(data2$Arrival)
data2$Community<-factor(data2$Community)

P1<-data.frame(Community=c(rep("G", 5), rep("GL", 5)), Replicate=c(c(1:5), c(1:5)), P=rep(NA, 10))

for (i in 1:10){
  Bji<-data2$NconcSenecio[data2$Arrival=="Late" & data2$Community==P1$Community[i] & data2$Replicate==P1$Replicate[i]]
  Bsync2<-data2$NconcSenecio[data2$Arrival=="Synchronous2" & data2$Community==P1$Community[i] & data2$Replicate==P1$Replicate[i]]
  P1$P[i]<-(Bji-Bsync2)/(Bji+Bsync2)}

confP1G<-CIbootstrap(x=P1$P[1:5])
confP1GL<-CIbootstrap(x=P1$P[6:10])

#For natives

data2<-data[data$Arrival=="Early"|data$Arrival=="Synchronous2",]
data2$Arrival<-factor(data2$Arrival)
data2$Community<-factor(data2$Community)

P3<-data.frame(Community=c(rep("G", 5), rep("GL", 5)), Replicate=c(c(1:5), c(1:5)), P=rep(NA, 10))

for (i in 1:10){
  Bij<-data2$NconcNatives[data2$Arrival=="Early" & data2$Community==P3$Community[i] & data2$Replicate==P3$Replicate[i]]
  Bsync2<-data2$NconcNatives[data2$Arrival=="Synchronous2" & data2$Community==P3$Community[i] & data2$Replicate==P3$Replicate[i]]
  P3$P[i]<-(Bij-Bsync2)/(Bij+Bsync2)}

confP3G<-CIbootstrap(x=P3$P[1:5])
confP3GL<-CIbootstrap(x=P3$P[6:10])

#Make one dataset

B<-rbind(P1,P3)
B$Species<-factor(c(rep("aSenecio", 10), rep("bNatives", 10)))

modelB<-lm(P~Species*Community, B)
#plot(fitted(modelB), resid(modelB), pch=1)
#abline(h=0)
anova(modelB)

tmp<-expand.grid(Community=levels(B$Community), Species=levels(B$Species))
X<-model.matrix(~Species*Community, data=tmp)
glht(modelB, linfct=X)
predict(modelB, newdata=tmp, type="response")
Tukey<-contrMat(table(B$Community), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(B$Species)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(B$Species)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))

summary(glht(modelB, linfct=K%*%X))

#Plot results for shoot N content

par(bty="l", mar=c(3,4,1.5,2)+0.1)
gf<-beanplot(P~Community*Species, B, side="both", log="", what=c(0,0,1,1), border = NA, 
             col = list(c("white", "black", "black", "black"),
                        c("white", "gray50", "gray50", "gray50")), 
             ylab=expression(paste("Cost of arriving late (", italic(P[E]), " or ", italic(P[N]), ")", sep="")), las=1, ylim=c(-1,1), xlab="", 
             cex.main=1, names=c(expression(italic("S. inaequidens")), "Natives"), at=c(1,2.5), xlim=c(0.5,3),
             yaxt="n", main="", ll=0.11, cex.main=1.3)
abline(h=0, lty=3, lwd=0.8)
axis(2, at=seq(-1,1, by=0.5), las=1)
arrows(0.75, confP1G[1], 0.75, confP1G[2], code=3, angle=90, length=0.02, col="black")
arrows(1.25, confP1GL[1], 1.25, confP1GL[2], code=3, angle=90, length=0.02, col="gray50")
arrows(2.25, confP3G[1], 2.25, confP3G[2], code=3, angle=90, length=0.02, col="black")
arrows(2.75, confP3GL[1], 2.75, confP3GL[2], code=3, angle=90, length=0.02, col="gray50")
text(0.55, 0.95, "B", font=2, cex=1.4)
text(1, -0.6, "ns", cex=1.1)
text(2.5, -0.1, "***", cex=2)
legend("topright", legend=c("Grasses", "Grasses-Legumes"), lty=1, lwd=2, col=c("black", "gray50"), bty="n", cex=0.9)

dev.off()

###########
###########
#Figure S5
###########
###########

SA<-0.28*0.28
data$NconcSenecio<-(SA*data$BSenecio*(data$NSenecio/100))/(SA*data$BSenecio)
data$NconcNatives<-NA
for (i in 1:nrow(data)){
  if(is.na(data$NLegumes[i])==TRUE){data$NconcNatives[i]<-(SA*data$BGrasses[i]*(data$NGrasses[i]/100))/(SA*data$BGrasses[i])}
  else {data$NconcNatives[i]<-(SA*data$BGrasses[i]*(data$NGrasses[i]/100)+SA*data$BLegumes[i]*(data$NLegumes[i]/100))/(SA*data$BGrasses[i]+SA*data$BLegumes[i])}}

#Synchronous 1

newdata<-data[data$Arrival=="Synchronous1",]
newdata<-data.frame(Species=c(rep("aSenecio", 10), rep("bNatives", 10)),
                    Community=rep(newdata$Community, 2),
                    N=c(newdata$NconcSenecio, newdata$NconcNatives))

model<-lm(N*1000~Community*Species, newdata)
#plot(fitted(model), resid(model)); abline(h=0)
anova(model)

tmp<-expand.grid(Species=levels(newdata$Species), Community=levels(newdata$Community))
X<-model.matrix(~Community*Species, data=tmp)
glht(model, linfct=X)
predict(model, newdata=tmp, type="response")
Tukey<-contrMat(table(newdata$Species), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(newdata$Community)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(newdata$Community)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))

summary(glht(model, linfct=K%*%X))

tiff(filename="FigureS5.tif", res=1000, compression="lzw", width=12, height=8, units="cm", pointsize=7)
layout(matrix(1:2, nrow=1, ncol=2))
par(bty="l", mar=c(5,4.5,4,1)+0.1)
beanplot(N*1000~Species*Community, newdata, side="both", what=c(1,0,1,1), border = NA, 
         col = list(c("white", "black", "black", "black"),
                    c("white", "gray50", "gray50", "gray50")), las=1, 
         ylab=expression(paste("N content (mg N ", g^{-1}, " shoot dry weight)")), names=c(expression("Grasses", "Grasses-Legumes")), 
         ylim=c(20,80), at=c(1,2.5), xlim=c(0.5,3), main="Synchronous 1", xlab="Plant age: 60 days")
text(1, 43, "ns", cex=1.1)
text(2.5, 43, "***", cex=2)
legend("topright", legend=c(expression(italic("S. inaequidens")), "Natives"), lwd=2, lty=1, col=c("black", "gray50"),
       border=NA, bty="n", cex=0.9)

#Synchronous 2

newdata<-data[data$Arrival=="Synchronous2",]
newdata<-data.frame(Species=c(rep("aSenecio", 10), rep("bNatives", 10)),
                    Community=rep(newdata$Community, 2),
                    N=c(newdata$NconcSenecio, newdata$NconcNatives))

model<-lm(N*1000~Community*Species, newdata)
#plot(fitted(model), resid(model)); abline(h=0)
anova(model)

tmp<-expand.grid(Species=levels(newdata$Species), Community=levels(newdata$Community))
X<-model.matrix(~Community*Species, data=tmp)
glht(model, linfct=X)
predict(model, newdata=tmp, type="response")
Tukey<-contrMat(table(newdata$Species), "Tukey")
K1<-cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(newdata$Community)[1], rownames(K1), sep=":")
K2<-cbind(matrix(0, nrow=nrow(Tukey), ncol=ncol(Tukey)), Tukey)
rownames(K2)<-paste(levels(newdata$Community)[2], rownames(K2), sep=":")
K<-rbind(K1, K2)
colnames(K)<-c(colnames(Tukey), colnames(Tukey))

summary(glht(model, linfct=K%*%X))

par(bty="l", mar=c(5,3.5,4,2)+0.1)
beanplot(N*1000~Species*Community, newdata, side="both", what=c(1,0,1,1), border = NA, 
         col = list(c("white", "black", "black", "black"),
                    c("white", "gray50", "gray50", "gray50")), las=1, 
         ylab="", names=c(expression("Grasses", "Grasses-Legumes")), 
         ylim=c(20,80), at=c(1,2.5), xlim=c(0.5,3), main="Synchronous 2", xlab="Plant age: 39 days")
text(1, 70, "ns", cex=1.1)
text(2.5, 70, "*", cex=2)

dev.off()

###########
###########
#Figure S6
###########
###########

newdata<-data.frame(Arrival=data$Arrival, Community=as.character(data$Community), Biomass=data$BLegumes+data$BGrasses)
newdata$Arrival<-factor(newdata$Arrival,levels(newdata$Arrival)[c(2,3,4,1)])

tiff(filename="FigureS6.tif", res=1000, compression="lzw", width=8, height=8, units="cm", pointsize=7)
par(bty="l")
gf<-beanplot(Biomass~Community*Arrival, newdata, side="both", log="", border = NA,
             ylab="Total native shoot dry weight (g/m²)", las=1, what=c(1,0,1,1), 
             col = list(c("white", "black", "black", "black"),
                        c("white", "gray50", "gray50", "gray50")),
             ylim=c(0,600), at=c(1,2.5,4,5.5), xlim=c(0.5,6),
             names=c("Late", expression(Sync[1]), expression(Sync[2]), "Early"),
             xlab="Timing of arrival of the exotic species")
legend("topright", legend=c("Grasses", "Grasses-Legumes"), lwd=2, lty=1, col=c("black", "gray50"), border=NA, bty="n",
       cex=0.9)
dev.off()
