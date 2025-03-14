library(tidyverse)
library(ggpubr)
library(rstatix)
library(lme4)
library(lmerTest)
library(multcomp)
library(emmeans)
library(performance)

# load data for each condition separately
data <- read.csv('PATH/file_Music.csv') 
data$frequency = as.factor(data$frequency)
data$condition = as.factor(data$condition)
data$subject   = as.factor(data$subject)
data$region    = as.factor(data$region)

# LMERs for r values
lme0 = lmer(r ~ (1|subject), data=data, REML=FALSE)
lme1 = lmer(r ~ frequency + (1|subject), data=data, REML=FALSE)
lme2 = lmer(r ~ region    + (1|subject), data=data, REML=FALSE)
lme3 = lmer(r ~ frequency + region + (1|subject), data=data, REML=FALSE)
lme4 = lmer(r ~ frequency*region + (1|subject), data=data, REML=FALSE)
lme5 = lmer(r ~ frequency*region + (1|subject:region), data=data, REML=FALSE)
lme6 = lmer(r ~ frequency*region + (1|subject) + (1|subject:region), data=data, REML=FALSE)

model2compare1 = ####
model2compare2 = ####
anova(model2compare1,model2compare2)

bestmodel = ####

summary(bestmodel)
anova(bestmodel)

# LMERs for lag values
lme0 = lmer(lag ~ (1|subject), data=data, REML=FALSE)
lme1 = lmer(lag ~ frequency + (1|subject), data=data, REML=FALSE)
lme2 = lmer(lag ~ region    + (1|subject), data=data, REML=FALSE)
lme3 = lmer(lag ~ frequency + region + (1|subject), data=data, REML=FALSE)
lme4 = lmer(lag ~ frequency*region + (1|subject), data=data, REML=FALSE)
lme5 = lmer(lag ~ frequency*region + (1|subject:region), data=data, REML=FALSE)
lme6 = lmer(lag ~ frequency*region + (1|subject) + (1|subject:region), data=data, REML=FALSE)

model2compare1 = ####
model2compare2 = ####
anova(model2compare1,model2compare2)

bestmodel = ####
  
summary(bestmodel)
anova(bestmodel)
