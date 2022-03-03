#Code to run ANOVAs and ANCOVA's on action prediction data

install.packages("tidyverse")
install.packages("gganimate")
install.packages("dplyr") 
install.packages("Hmisc")
install.packages("emmeans")
install.packages("afex")
install.packages("afex", dependencies = TRUE) 
install.packages("ggpubr")
install.packages("rstatix")
install.packages("broom")
install.packages("coda") #Needed for Bayes factors
install.packages("Matrix") #Needed for Bayes factors
install.packages("BayesFactor")
install.packages("bayestestR")
install.packages("naniar")
install.packages("car")
install.packages("rcompanion")# for effect sizes
install.packages("effsize")# for effect sizes
install.packages("effectsize")# for effect sizes
#install.packages("coin")# for mann whitney effect sizes
install.packages("lsr")#for cohen d effect sizes
install.packages("fitdistrplus")


library(gganimate)
library(Hmisc) #needed for correlation
library(emmeans)# needed for pairwise comparisons
library(afex) # needed for ANOVA
library(naniar) #for replacing NaN with NA
library(lsr) 

# needed for examining normal distribution etc (Levenes/sharpiro wilks)
library(ggpubr) 
library(rstatix)
library(broom)
library(coda)
library(Matrix)
library(BayesFactor)
library(bayestestR)
library(car)
library(fitdistrplus)
#library(effsize)
library(effectsize)# for effect sizes
#library(coin)


library(tidyverse)
library(dplyr) #needed for %>%


#Participant demographics################################################################################################
#Comparing autistic and non autistic groups on demographic information

remove(list = ls())
my_data_original <- read.csv(file.choose()) #loads window to choose your file. Use "Participants" spreadsheet

my_data_original_summary <- my_data_original%>% 
  group_by(Group) %>% 
  summarise(meanAge = mean(AGE, na.rm = TRUE), 
            sdAge = sd(AGE,na.rm = TRUE),
            meanELI = mean(ELI, na.rm = TRUE), 
            sdELI = sd(ELI,na.rm = TRUE),
            meanVIQ = mean(Verbal.IQ, na.rm = TRUE), 
            sdVIQ = sd(Verbal.IQ,na.rm = TRUE),
            meanPIQ = mean(PerformanceIQ, na.rm = TRUE), 
            sdPIQ = sd(PerformanceIQ,na.rm = TRUE),
            meanFIQ = mean(Full.scale.IQ, na.rm = TRUE), 
            sdFIQ = sd(Full.scale.IQ,na.rm = TRUE))


var.test(my_data_original$AGE ~ my_data_original$Group) #This tests whether variances are equal
var.test(my_data_original$Verbal.IQ ~ my_data_original$Group) #This tests whether variances are equal
var.test(my_data_original$PerformanceIQ ~ my_data_original$Group) #This tests whether variances are equal
var.test(my_data_original$Full.scale.IQ ~ my_data_original$Group) #This tests whether variances are equal

t.test(my_data_original$AGE ~ my_data_original$Group, var.equal = TRUE)
t.test(my_data_original$Verbal.IQ ~ my_data_original$Group, var.equal = TRUE)
t.test(my_data_original$PerformanceIQ ~ my_data_original$Group, var.equal = TRUE)
t.test(my_data_original$Full.scale.IQ ~ my_data_original$Group, var.equal = TRUE)



#Accuracy################################################################################################
#Action prediction accuracy spreadsheet
remove(list = ls())
my_data_original <- read.csv(file.choose()) #loads window to choose your file. Action prediction accuracy spreadsheet
colnames(my_data_original) #Checking column names of data

#changing accuracy to percentage
my_data_original$Accuracy <- my_data_original$Accuracy*100  #Turning accuracy to %


#Turning variables into factors
my_data_original$Group <- as.factor(my_data_original$Group)#change autism to a factor
levels(my_data_original$Group)#check number of levels
my_data_original$Continuation <- as.factor(my_data_original$Continuation)#change autism to a factor
levels(my_data_original$Continuation)#check number of levels

#Renaming data
my_data_original <- my_data_original %>%
  mutate(Group = recode(Group,
                         "1" = "Autistic",
                         "0" = "Non-autistic"))

my_data_original <- my_data_original %>%
  mutate(Continuation = recode(Continuation,
                        "Early_680" = "Too far behind 680",
                        "Early_840" = "Too far behind 840",
                        "In_time" = "In time",
                        "Late_680" = "Too far ahead 680",
                        "Late_840" = "Too far ahead 840"))

#Reordering factors for graph
my_data_original$Continuation<- factor(my_data_original$Continuation,levels = c("Too far behind 840", "Too far behind 680", "In time", "Too far ahead 680", "Too far ahead 840"))
levels(my_data_original$Continuation)#check number of levels

#plot the data
ggbarplot(my_data_original, x = "Continuation", y = "Accuracy", 
          add = c("mean_se",  "jitter"), # Add mean_se and jitter points
          fill="Group", color = "Group", alpha = .5, palette = c("blue", "red"),  #Adding colours, alpha changes transparency
          position = position_dodge(0.8),#separting bars
          ylab = "% correct",
          ylim = c(0, 100)) +
  
  theme( #changing fonts
    
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black")) +   #font of all axis labels apart from titles
  scale_y_continuous(breaks=seq(0, 100, 10)) # change y scale intervals


#Checking distribution
model_anova_accuracy_lm <- lm(Accuracy ~ Group + Continuation, data = my_data_original) #Create model so have residuals

hist(residuals(model_anova_accuracy_lm))
qqnorm(residuals(model_anova_accuracy_lm))
qqline(residuals(model_anova_accuracy_lm))
plot(model_anova_accuracy_lm) #checking for Homoskedasticity - needs to be blob like. 
var.test(my_data_original$Accuracy ~ my_data_original$Group)

# Checking normality
shapiro.test(model_anova_accuracy_lm$residuals)
# checking variances of residuals is homogeneous
leveneTest(Accuracy ~ Group,
           data = my_data_original)
           
# Shapiro wilks sig:  0.970 0.000431. Data is left skewed

#Transform data
model_anova_accuracy_lm_log <- lm(log ~ Group + Continuation, data = my_data_original)#Create model so have residuals
# cant do log as have meaningful zeroes in data.  
model_anova_accuracy_lm_sq <- lm(Accuracy^2 ~ Group + Continuation, data = my_data_original) #Create model so have residuals

#Square seems better (square data for left skew, square root for right skew)

hist(residuals(model_anova_accuracy_lm_sq))
qqnorm(residuals(model_anova_accuracy_lm_sq))
qqline(residuals(model_anova_accuracy_lm_sq))
plot(model_anova_accuracy_lm_sq) #checking for Homoskedasticity - needs to be blob like. 


# Checking normality
shapiro.test(model_anova_accuracy_lm_sq$residuals)
# checking variances of residuals is homogeneous
leveneTest(Accuracy ~ Group,
           data = my_data_original)
#Shapiro Wilkes still sig: 0.978 0.006

#As data still non-normal need to use non-parametric tests
#Mann Whitney
MW=wilcox.test(Accuracy ~ Group, data=my_data_original) 
MW
#W(U) = 5133.5, p-value = 0.01824

#effect size using r
#obtain z value
Za = qnorm(MW$p.value/2)

r_effect_size = abs(Za)/sqrt(36) #(z value/sqrt(N))
r_effect_size
#0.39 (medium)


#Friedman for each group

#load in separate data for each group
my_data_original_autistic <- read.csv(file.choose()) #Action_prediction_accuracy_autistic
my_data_original_non_autistic <- read.csv(file.choose()) #Action_prediction_accuracy_non_autistic

my_data_original_autistic <- my_data_original_autistic %>%
  mutate(Continuation = recode(Continuation,
                               "Early_680" = "Too far behind 680",
                               "Early_840" = "Too far behind 840",
                               "In_time" = "In time",
                               "Late_680" = "Too far ahead 680",
                               "Late_840" = "Too far ahead 840"))

my_data_original_non_autistic <- my_data_original_non_autistic %>%
  mutate(Continuation = recode(Continuation,
                               "Early_680" = "Too far behind 680",
                               "Early_840" = "Too far behind 840",
                               "In_time" = "In time",
                               "Late_680" = "Too far ahead 680",
                               "Late_840" = "Too far ahead 840"))

#changing accuracy to percentage
my_data_original_autistic$Accuracy <- my_data_original_autistic$Accuracy*100
my_data_original_non_autistic$Accuracy <- my_data_original_non_autistic$Accuracy*100


my_data_original_autistic%>%friedman_test(Accuracy ~ Continuation |Participant)
#15.7     (4) 0.00342
my_data_original_autistic %>% friedman_effsize(Accuracy ~ Continuation |Participant) #Kendall’s W uses the Cohen’s interpretation guidelines of 0.1 - < 0.3 (small effect), 0.3 - < 0.5 (moderate effect) and >= 0.5 (large effect).
#W=0.22 (small)
my_data_original_autistic %>% wilcox_test(Accuracy ~ Continuation, paired = TRUE, p.adjust.method = "bonferroni")
#These are sig: Early_680 In_time (0.003); early 840-in time (0.002)
my_data_original_autistic  %>% wilcox_effsize(Accuracy ~ Continuation, paired = TRUE) # r effect size: The effect size r is calculated as Z statistic divided by square root of the sample size 

my_data_original_non_autistic%>%friedman_test(Accuracy ~ Continuation |Participant)
#1.67     (4) 0.80
my_data_original_non_autistic %>% friedman_effsize(Accuracy ~ Continuation |Participant)
#W=0.02

#Means
my_data_original%>% 
  group_by(Group) %>% 
  summarise(meanAccuracy = mean(Accuracy, na.rm = TRUE), 
            sdAccuracy = sd(Accuracy,na.rm = TRUE))

my_data_original%>% 
  group_by(Group,Continuation) %>% 
  summarise(meanAccuracy = mean(Accuracy, na.rm = TRUE), 
            sdAccuracy = sd(Accuracy,na.rm = TRUE))



#Timing biases################################################################################
#Action_prediction_timing_bias
remove(list = ls())
my_data_original <- read.csv(file.choose()) #loads window to choose your file. Action_prediction_timing_bias
colnames(my_data_original) #Checking column names of data

#changing timing bias to percentage
my_data_original$Timing <- my_data_original$Timing*100  #Turning in time responses to %


#Turning variables into factors
my_data_original$Group <- as.factor(my_data_original$Group)#change autism to a factor
levels(my_data_original$Group)#check number of levels
my_data_original$Continuation <- as.factor(my_data_original$Continuation)#change autism to a factor
levels(my_data_original$Continuation)#check number of levels

#Renaming data
my_data_original <- my_data_original %>%
  mutate(Group = recode(Group,
                        "1" = "Autistic",
                        "0" = "Non-autistic"))

my_data_original <- my_data_original %>%
  mutate(Continuation = recode(Continuation,
                               "Early_680" = "Too far behind 680",
                               "Early_840" = "Too far behind 840",
                               "In_time" = "In time",
                               "Late_680" = "Too far ahead 680",
                               "Late_840" = "Too far ahead 840"))

#Reordering factors for graph
my_data_original$Continuation<- factor(my_data_original$Continuation,levels = c("Too far behind 840", "Too far behind 680", "In time", "Too far ahead 680", "Too far ahead 840"))
levels(my_data_original$Continuation)#check number of levels

#plot the data
ggbarplot(my_data_original, x = "Continuation", y = "Timing", 
          add = c("mean_se",  "jitter"), # Add mean_se and jitter points
          fill="Group", color = "Group", alpha = .5, palette = c("blue", "red"),  #Adding colours, alpha changes transparency
          position = position_dodge(0.8),#separting bars
          ylab = "% just in time responses",
          ylim = c(0, 100)) +
  
  theme( #changing fonts
    
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black")) +   #font of all axis labels apart from titles
  scale_y_continuous(breaks=seq(0, 100, 10)) # change y scale intervals

#Looking at distribution
model_anova_Timing_lm <- lm(Timing ~ Group * Continuation, data = my_data_original) #Create model so have residuals

hist(residuals(model_anova_Timing_lm))
qqnorm(residuals(model_anova_Timing_lm))
qqline(residuals(model_anova_Timing_lm))
plot(model_anova_Timing_lm) #checking for Homoskedasticity - needs to be blob like. 

# Checking normality
shapiro.test(model_anova_Timing_lm$residuals)
# checking variances of residuals is homogeneous
leveneTest(Timing ~ Group,
           data = my_data_original)
# Shapiro wilks not sig:  p=0.2,, data looks normal

#Run ANOVA 
model_anova_timing <- aov_4(Timing ~ Group*Continuation + (1+ Continuation | Participant), data = my_data_original)
summary(model_anova_timing)
anova(model_anova_timing) 

BF <- generalTestBF(Timing ~ Group * Continuation, my_data_original, progress = FALSE)
bayesfactor_inclusion(BF)

emmeans(model_anova_timing, pairwise ~ Continuation)
#getting effect sizes
t_to_d(t = c(-5.67,-4.47,3.19, 9.57,10.4),
       df_error = 128,
       paired = TRUE)


###################################################################################################

#Performing orthogonal polynomial contrasts

#https://blogs.uoregon.edu/rclub/2015/11/03/anova-contrasts-in-r/
#spreadsheets action_prediction_timing_bias_autistic and action_prediction_timing_bias_non_autistic
remove(list = ls())
#non-autistic data
my_data_original_non_autistic <- read.csv(file.choose()) #load action_prediction_timing_bias_non_autistic
my_data_original_non_autistic$Continuation <- as.factor(my_data_original_non_autistic$Continuation)#change Continuation to a factor
contrasts(my_data_original_non_autistic$Continuation) <- contr.poly(5) #tell r to get a polynomial contrast matrix for 5 comparisons

model_non_autistic <- aov(Timing_bias ~ Continuation, data = my_data_original_non_autistic)
summary.aov(model_non_autistic, split=list(Continuation=list("Linear"=1, "Quadratic" = 2, "Cubic" = 3)))

#autistic data
my_data_original_autistic <- read.csv(file.choose()) #load action_prediction_timing_bias_autistic
my_data_original_autistic$Continuation <- as.factor(my_data_original_autistic$Continuation)#change Continuation to a factor
contrasts(my_data_original_autistic$Continuation) <- contr.poly(5) #tell r to get a polynomial contrast matrix for 5 comparisons
model_autistic <- aov(Timing_bias ~ Continuation, data = my_data_original_autistic)
summary.aov(model_autistic, split=list(Continuation=list("Linear"=1, "Quadratic" = 2, "Cubic" = 3)))


#Questionnaire data############################################################
remove(list = ls())
#Questionnaires spreadsheet
my_data_original_questionnaires <- read.csv(file.choose()) #loads window to choose your file - Questionnaires.csv
colnames(my_data_original_questionnaires) #Checking column names of data

#changing accuracy to percentage
my_data_original_questionnaires$Mean.Accuracy <- my_data_original_questionnaires$Mean.Accuracy*100
my_data_original_questionnaires$Median.Accuracy <- my_data_original_questionnaires$Median.Accuracy*100


#Checking data for distribution
#Familiarity
hist(my_data_original_questionnaires$Familiarity)
qqnorm(my_data_original_questionnaires$Familiarity)
qqline(my_data_original_questionnaires$Familiarity)
plot(my_data_original_questionnaires$Familiarity) #checking for Homoskedasticity - needs to be blob like. 
shapiro_test(my_data_original_questionnaires$Familiarity)

#Shapriro Wilkes: 0.902 0.00336

#KVIQ.V
hist(my_data_original_questionnaires$KVIQ.V)
qqnorm(my_data_original_questionnaires$KVIQ.V)
qqline(my_data_original_questionnaires$KVIQ.V)
plot(my_data_original_questionnaires$KVIQ.V) #checking for Homoskedasticity - needs to be blob like. 
shapiro_test(my_data_original_questionnaires$KVIQ.V)

#Shapriro Wilkes:0.921  0.0189

#KVIQ.K
hist(my_data_original_questionnaires$KVIQ.K)
qqnorm(my_data_original_questionnaires$KVIQ.K)
qqline(my_data_original_questionnaires$KVIQ.K)
plot(my_data_original_questionnaires$KVIQ.K) #checking for Homoskedasticity - needs to be blob like. 
shapiro_test(my_data_original_questionnaires$KVIQ.K)

#Shapriro Wilkes: 0.97   0.33

#DCD.checklist
hist(my_data_original_questionnaires$DCD.checklist)
qqnorm(my_data_original_questionnaires$DCD.checklist)
qqline(my_data_original_questionnaires$DCD.checklist)
plot(my_data_original_questionnaires$DCD.checklist) #checking for Homoskedasticity - needs to be blob like. 
shapiro_test(my_data_original_questionnaires$DCD.checklist)

#Shapriro Wilkes:0.95   0.09

#KVIQ and Familiarity are non-normally distributed so use mann whitney for these

# Means
my_data_original_questionnaires%>% 
  group_by(Group) %>% 
  summarise(meanFam = mean(Familiarity, na.rm = TRUE), 
            sdFam = sd(Familiarity,na.rm = TRUE),
            meanKVIQ_V = mean(KVIQ.V, na.rm = TRUE),
            sdKVIQ_V = sd(KVIQ.V,na.rm = TRUE),
            meanKVIQ_K = mean(KVIQ.K, na.rm = TRUE),
            sdKVIQ_K = sd(KVIQ.K,na.rm = TRUE),
            meanDCD = mean(DCD.checklist, na.rm = TRUE),
            sdDCD = sd(DCD.checklist,na.rm = TRUE))

# Medians
my_data_original_questionnaires%>% 
  group_by(Group) %>% 
  summarise(medianFam = median(Familiarity, na.rm = TRUE), 
            medianKVIQ_V = median(KVIQ.V, na.rm = TRUE),
            medianKVIQ_K = median(KVIQ.K, na.rm = TRUE),
            medianDCD = median(DCD.checklist, na.rm = TRUE))


MW =wilcox.test(Familiarity ~ Group, correct=FALSE, data=my_data_original_questionnaires) 
MW
#W = 256.5, p-value = 0.009
#effect size using r
#obtain z value
Za = qnorm(MW$p.value/2)

r_effect_size = abs(Za)/sqrt(36) #(z value/sqrt(N))
r_effect_size
#0.43 


MW_KVIQ=wilcox.test(KVIQ.V ~ Group, correct=FALSE, data=my_data_original_questionnaires) 
MW_KVIQ

#W = 219, p-value = 0.14
#effect size using r
#obtain z value
Za_KVIQ = qnorm(MW_KVIQ$p.value/2)
r_effect_size_KVIQ = abs(Za_KVIQ)/sqrt(36) #(z value/sqrt(N))
r_effect_size_KVIQ


# t tests to compare groups on DCD/KVI-Kquestionnaires
var.test(my_data_original_questionnaires$KVIQ.K ~ my_data_original_questionnaires$Group) #This tests whether variances are equal
var.test(my_data_original_questionnaires$DCD.checklist ~ my_data_original_questionnaires$Group) #This tests whether variances are equal


t.test(my_data_original_questionnaires$KVIQ.K ~ my_data_original_questionnaires$Group)
bf = ttestBF(formula = KVIQ.K ~ Group, data = my_data_original_questionnaires)
bf
cohensD(KVIQ.K ~ Group,data = my_data_original_questionnaires)
#t = -0.001, df = 27.658, p-value = 0.9989, BF=0.32, d=<0.001

t.test(my_data_original_questionnaires$DCD.checklist ~ my_data_original_questionnaires$Group,var.equal = TRUE)
bf = ttestBF(formula = DCD.checklist ~ Group, data = my_data_original_questionnaires)
bf
cohensD(DCD.checklist ~ Group,data = my_data_original_questionnaires)
#t = -9.68, df = 34, p-value = 1.988e-11;BF= 276795044


#Correlations  

#Separate for autistic and non-autistic to perform separate correlations

my_data_original_questionnaires_non_autistic<-filter(my_data_original_questionnaires, Group == "0")
my_data_original_questionnaires_autistic<-filter(my_data_original_questionnaires, Group == "1")

# Autistic. Spearmans
my_data_correlations_autistic_spearmans <-select(my_data_original_questionnaires_autistic, Familiarity, KVIQ.V, Mean.Accuracy, Median.Accuracy,Group) 
my_data_correlations_autistic2_spearmans = rcorr(as.matrix(my_data_correlations_autistic_spearmans),type=c("spearman"))# compute correlations
my_data_correlations_autistic2_spearmans #show correlations
mydata.coeff.autistic_spearmans = round(my_data_correlations_autistic2_spearmans$r,digits=5) #show correlations in table
mydata.p.autistic_spearmans = round(my_data_correlations_autistic2_spearmans$P,digits=5)
#non sig

# Non Autistic. Spearmans
my_data_correlations_non_autistic_spearmans <-select(my_data_original_questionnaires_non_autistic, Familiarity, KVIQ.V, Mean.Accuracy, Median.Accuracy,Group) 
my_data_correlations_non_autistic2_spearmans = rcorr(as.matrix(my_data_correlations_non_autistic_spearmans),type=c("spearman"))# compute correlations
my_data_correlations_non_autistic2_spearmans #show correlations
mydata.coeff.non_autistic_spearmans = round(my_data_correlations_non_autistic2_spearmans$r,digits=5) #show correlations in table
mydata.p.non_autistic_spearmans = round(my_data_correlations_non_autistic2_spearmans$P,digits=5)
#non sig


# Autistic. Pearsons
my_data_correlations_autistic_pearsons <-select(my_data_original_questionnaires_autistic, KVIQ.K, DCD.checklist, Mean.Accuracy, Median.Accuracy,Group) 
my_data_correlations_autistic2_pearsons = rcorr(as.matrix(my_data_correlations_autistic_pearsons),type=c("spearman"))# compute correlations
my_data_correlations_autistic2_pearsons #show correlations
mydata.coeff.autistic_pearsons = round(my_data_correlations_autistic2_pearsons$r,digits=5) #show correlations in table
mydata.p.autistic_pearsons = round(my_data_correlations_autistic2_pearsons$P,digits=5)
#non sig



# Non Autistic. Pearsons
my_data_correlations_non_autistic_pearsons <-select(my_data_original_questionnaires_non_autistic, KVIQ.K, DCD.checklist, Mean.Accuracy, Median.Accuracy,Group) 
my_data_correlations_non_autistic2_pearsons = rcorr(as.matrix(my_data_correlations_non_autistic_pearsons),type=c("spearman"))# compute correlations
my_data_correlations_non_autistic2_pearsons #show correlations
mydata.coeff.non_autistic_pearsons = round(my_data_correlations_non_autistic2_pearsons$r,digits=5) #show correlations in table
mydata.p.non_autistic_pearsons = round(my_data_correlations_non_autistic2_pearsons$P,digits=5)
#non sig

#plotting KVIQ-K and accuracy for non-autistics as almost sig
ggplot(my_data_original_questionnaires_non_autistic, aes(x = KVIQ.K, y = Mean.Accuracy)) + 
  geom_point() + 
  geom_smooth(method = "lm") +  #This adds a line of best fit
  guides(colour = FALSE)



##Control task#########################################################################
#Control task spreadsheet
remove(list = ls())
my_data_original_control_task <- read.csv(file.choose()) #load control_tasks.csv 
colnames(my_data_original_control_task) #Checking column names of data

#Control task
#Checking data for distribution

hist(my_data_original_control_task$Jump.present)
qqnorm(my_data_original_control_task$Jump.present)
qqline(my_data_original_control_task$Jump.present)
plot(my_data_original_control_task$Jump.present) #checking for Homoskedasticity - needs to be blob like. 
shapiro_test(my_data_original_control_task$Jump.present)

#Not normally distributed - shaprio <0.001

MW_control_task =wilcox.test(Jump.present ~ Group, correct=FALSE, data=my_data_original_control_task) 
MW_control_task
#W = 146, p-value = 0.43
Za_control_task = qnorm(MW_control_task$p.value/2)

r_effect_size_control_task = abs(Za_control_task)/sqrt(36) #(z value/sqrt(N))
r_effect_size_control_task
#0.13

my_data_original_control_task%>% 
  group_by(Group) %>% 
  summarise(meanJump = mean(Jump.present, na.rm = TRUE), 
            sdJump = sd(Jump.present,na.rm = TRUE),
            medianJump = median(Jump.present, na.rm = TRUE),
            sdJump = sd(Jump.present, na.rm = TRUE))

#Non-autistic percentage correct:          
21.9/24
#autistic percentage correct:          
22.9/24


#Dwell time
hist(my_data_original_control_task$Percentage.Dwell.time)
qqnorm(my_data_original_control_task$Percentage.Dwell.time)
qqline(my_data_original_control_task$Percentage.Dwell.time)
plot(my_data_original_control_task$Percentage.Dwell.time) #checking for Homoskedasticity - needs to be blob like. 
shapiro_test(my_data_original_control_task$Percentage.Dwell.time)

#Not normally distributed. Shapiro <0.001

MW_Dwell =wilcox.test(Percentage.Dwell.time ~ Group, correct=FALSE, data=my_data_original_control_task) 
MW_Dwell
#W = 166, p-value = 0.66

#calculate effect size
Za_Dwell = qnorm(MW_Dwell$p.value/2)

r_effect_size_Dwell = abs(Za_Dwell)/sqrt(36) #(z value/sqrt(N))
r_effect_size_Dwell
#0.07

my_data_original_control_task%>% 
  group_by(Group) %>% 
  summarise(meanDwell = mean(Percentage.Dwell.time, na.rm = TRUE), 
            sdDwell = sd(Percentage.Dwell.time,na.rm = TRUE),
            medianDwell = median(Percentage.Dwell.time, na.rm = TRUE),
            sdDwell = sd(Percentage.Dwell.time, na.rm = TRUE))






