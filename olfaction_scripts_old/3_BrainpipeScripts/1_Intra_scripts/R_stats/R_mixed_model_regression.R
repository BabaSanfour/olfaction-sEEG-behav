install.packages('lme4')
install.packages("ggplot2")
install.packages('sjPlot', dependencies = TRUE)
install.packages('sjPlot', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('MVN')

library(lme4)
library(ggplot2)
library(sjPlot)
library(sjmisc)
library(MVN)

MyData <- read.csv(file="/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/stats_R/csv/0_theta_HC_all_su_mean_trials.csv", header=TRUE, sep=",")

MyData_postulats <- data.frame('power.z'=MyData$power.z)

MyData$power.z <- scale(MyData$power, center=TRUE, scale=TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = glmer(score ~ power.z + (1|subj), data=MyData, family='poisson')
null = glmer(score ~ (1|subj), data=MyData, family='poisson')
anova(null, model)
summary(model)
plot_model(model, type='int', mdrt.values='all')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

######### POSTULATS
library(moments)
skewness(MyData$power.z)

library(psych)
describe(MyData_postulats)
