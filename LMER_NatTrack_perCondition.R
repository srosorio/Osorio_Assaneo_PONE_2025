library(tidyverse)
library(ggpubr)
library(rstatix)
library(lme4)
library(lmerTest)
library(multcomp)
library(emmeans)
library(performance)

# load data for each specific condition separately
data <- read.csv('PATH/file_Speech.csv') 
data$frequency = as.factor(data$frequency)
data$condition = as.factor(data$condition)
data$subject   = as.factor(data$subject)
data$region    = as.factor(data$region)

# LMER for r values
# Best model for Speech: r ~ frequency*region + (1|subject) + (1|subject:region) 
# Best model for Music:  r ~ frequency + (1|subject)
bestmodel = lmer(r ~ frequency*region + (1|subject) + (1|subject:region), data=data, REML=FALSE)

anova(bestmodel)
shapiro.test(residuals(bestmodel))
# pairwise comparison
emmeans(bestmodel, specs =  pairwise ~ frequency|region)

# LMER for lag values
# Best model for Speech: lag ~ frequency*region + (1|subject) + (1|subject:region)
# Best model for Music:  lag ~ (1+subject)
bestmodel = lmer(lag ~  frequency*region + (1|subject) + (1|subject:region), data=data, REML=FALSE)
summary(bestmodel)
shapiro.test(residuals(bestmodel))
# pairwise comparison
emmeans(bestmodel, specs = pairwise ~ frequency|region)
