library(tidyverse)
library(ggpubr)
library(rstatix)
library(lme4)
library(lmerTest)
library(multcomp)
library(emmeans)
library(performance)

# load data
data <- read.csv('PATH/file_joint.csv') 
data$frequency = as.factor(data$frequency)
data$condition = as.factor(data$condition)
data$subject   = as.factor(data$subject)
data$region    = as.factor(data$region)

# LMER for r values
bestmodel = lmer(r ~ frequency*region*condition + (1|subject) + (1|subject:region), data=data, REML=FALSE)
summary(bestmodel)
anova(bestmodel)
# pairwise comparison
emmeans(bestmodel, specs =  ~ condition)

# LMER for lag values
bestmodel = lmer(lag ~ frequency*region*condition + (1|subject), data=data, REML=FALSE)
summary(bestmodel)
anova(bestmodel)
# pairwise comparison
emmeans(bestmodel, specs = pairwise ~ frequency:condition)
