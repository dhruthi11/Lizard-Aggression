# Mandavilli & Brock
# August 2022

# Set working directory:
setwd("~/Desktop/Mandavilli_analysis/")

# Load packages:
library(dplyr)
library(tidyr)
library(readxl)
library(lubridate)
library(ggplot2)
library(car)
library(ggpubr)
library(FSA)
library(chisq.posthoc.test)

# Read in data:
data <- read.csv("data_2.csv")

# Clean and subset data:
data_select <- select(data, focal_morph, bite, chase, gape, lunge, int_type)
#bite counts
bite_counts <- data_select %>% group_by(int_type, focal_morph) %>% summarize(sum_bite = sum(bite))
bite_counts

#chase counts
chase_counts <- data_select %>% group_by(int_type, focal_morph) %>% summarize(sum_chase = sum(chase))
chase_counts

#gape counts
gape_counts <- data_select %>% group_by(int_type, focal_morph) %>% summarize(sum_gape = sum(gape))
gape_counts

#lunge counts
lunge_counts <- data_select %>% group_by(int_type, focal_morph) %>% summarize(sum_lunge = sum(lunge))
lunge_counts

# Data quick summary stats:
group_by(data, focal_morph) %>% 
  summarise(
    count = n(),
    mean = mean(agg_score, na.rm=T),
    sd = sd(agg_score, na.rm = T)
  )


# Aggressive Behavior Chi-Square/Bonferroni Tests --------------------------------------------------------

bite <- c(2, 25, 5)
res.bite <- chisq.test(bite, p = c(1/3,1/3,1/3))
res.bite
#X-squared = 29.312, df = 2, p-value = 4.314e-07
#There is a significant difference in one or more bite counts.

pairwise.t.test(data$bite, data$focal_morph, p.adj='bonferroni')
#o         w     
#w 0.0034  -     
#y 1.0000  0.0230
#White biting significantly differs with orange and yellow biting. Orange versus yellow biting is not significantly different.


chase <- c(1, 52, 6)
res.chase <- chisq.test(chase, p = c(1/3,1/3,1/3))
res.chase
#X-squared = 80.373, df = 2, p-value < 2.2e-16
#There is a significant difference in one or more chase counts.

pairwise.t.test(data$chase, data$focal_morph, p.adj='bonferroni')
#o          w      
#w 4.9e-12  -      
#y 0.47     4.2e-10
#White chasing significantly differs with orange and yellow chasing. Orange versus yellow chasing is not significantly different.


gape <- c(7, 16, 7)
res.gape <- chisq.test(gape, p = c(1/3,1/3,1/3))
res.gape
#X-squared = 5.4, df = 2, p-value = 0.06721
#There is no significant difference between gape counts.


lunge <- c(2, 35, 2)
res.lunge <- chisq.test(lunge, p = c(1/3,1/3,1/3))
res.lunge
#X-squared = 55.846, df = 2, p-value = 7.467e-13
#There is a significant difference in one or more lunge counts.

pairwise.t.test(data$lunge, data$focal_morph, p.adj='bonferroni')
#o          w      
#w 2.5e-06  -      
#y 1        6.5e-07
#White lunging significantly differs with orange and yellow lunging. Orange versus yellow lunging is not significantly different.


# Focal Morph - ANOVA/Kruskal-Wallis: Aggression Scores ------------------------------------------------

#Shapiro-Wilk normality test - White
white_agressionsc <- data %>% filter(focal_morph == "w") %>% pull(agg_score)
white_agressionsc
normality_W <- shapiro.test(white_agressionsc)
normality_W
#W = 0.96498, p-value = 0.4761 - Assume normality

#Shapiro-Wilk normality test - Yellow
yellow_agressionsc <- data %>% filter(focal_morph == "y") %>% pull(agg_score)
yellow_agressionsc
normality_Y <- shapiro.test(yellow_agressionsc)
normality_Y 
#W = 0.7482, p-value = 0.0008564 - Cannot assume normality

#Shapiro-Wilk normality test - Orange
orange_agressionsc <- data %>% filter(focal_morph == "o") %>% pull(agg_score)
orange_agressionsc
normality_O <- shapiro.test(orange_agressionsc)
normality_O 
#W = 0.80897, p-value = 0.008728 - Cannot assume normality

#Levene Test for variance
levene_agsc <- leveneTest(agg_score ~ focal_morph, data = data)
levene_agsc
#P-value = 0.02401
#There is a significant difference in the variance of one or more groups.

#ANNOVA - data is not normal/not same variance
#annova_agsc <- aov(agg_score ~ focal_morph, data = data)
#annova_agsc
#Residual standard error: 4.140543
#Estimated effects may be unbalanced

#Kruskal-Wallis - rank sum test
kruskal_agsc <- kruskal.test(agg_score ~ focal_morph, data = data)
kruskal_agsc
#Kruskal-Wallis chi-squared = 32.36, df = 2, p-value = 9.401e-08
#One or more morphs has a significantly different mean aggression score.

#DunnTest
dunn_agsc <- dunnTest(agg_score ~ focal_morph, data = data)
dunn_agsc
# Comparison          Z      P.unadj        P.adj
#1      o - w -4.8859472 1.029328e-06 3.087984e-06
#2      o - y -0.6103886 5.416044e-01 5.416044e-01
#3      w - y  4.4036221 1.064583e-05 2.129167e-05 

# White vs. Yellow and Orange vs. White average aggression scores are significantly different.
# Orange vs. Yellow aggression scores are not significantly different.


# Receiver Morph - ANOVA/Kruskal-Wallis: Aggression Scores ------------------------------------------------

#Shapiro-Wilk normality test - White
white_rec_agressionsc <- data %>% filter(receiver_morph == "w") %>% pull(agg_score)
white_rec_agressionsc
normality_W_rec <- shapiro.test(white_rec_agressionsc)
normality_W_rec
#W = 0.92378, p-value = 0.03364 - Cannot assume normality

#Shapiro-Wilk normality test - Yellow
yellow_rec_agressionsc <- data %>% filter(receiver_morph == "y") %>% pull(agg_score)
yellow_rec_agressionsc
normality_Y_rec <- shapiro.test(yellow_rec_agressionsc)
normality_Y_rec 
#W = 0.82502, p-value = 0.003504 - Cannot assume normality

#Shapiro-Wilk normality test - Orange
orange_rec_agressionsc <- data %>% filter(receiver_morph == "o") %>% pull(agg_score)
orange_rec_agressionsc
normality_O_rec <- shapiro.test(orange_rec_agressionsc)
normality_O_rec 
#W = 0.86597, p-value = 0.1711 - Assume normality

#Levene Test for variance
levene_agsc_rec <- leveneTest(agg_score ~ receiver_morph, data = data)
levene_agsc_rec
#P-value = 0.1459
#There is no significant difference in the variance of one or more groups.

#ANNOVA - data is not normal

#Kruskal-Wallis - rank sum test
kruskal_agsc_rec <- kruskal.test(agg_score ~ receiver_morph, data = data)
kruskal_agsc_rec
#Kruskal-Wallis chi-squared = 0.083806, df = 2, p-value = 0.959
#No morphs have a significantly different mean aggression score.

#DunnTest
dunn_agsc_rec <- dunnTest(agg_score ~ receiver_morph, data = data)
dunn_agsc_rec
# Comparison        Z      P.unadj    P.adj
#1      o - w -0.12480688 0.9006764 1.0000000
#2      o - y  0.07315734 0.9416809 0.9416809
#3      w - y  0.28501354 0.7756338 1.0000000

#No morphs have a significantly different mean aggression score.


# Graphs --------------------------------------------------------

data_bite <- read.csv("data_bite.csv")

ggboxplot(data, x = "focal_morph", y = "agg_score", 
          add = "jitter",
          fill = "focal_morph", 
          palette=c('white','gold','orange'),
          ylab = "Aggressor morph mean aggression score", xlab = "Morph")

ggboxplot(data, x = "receiver_morph", y = "agg_score",
          add = "jitter",
          fill = "receiver_morph", 
          palette=c('orange','white','gold'),
          ylab = "Receiver morph mean aggression score", xlab = "Morph")

ggbarplot(data_bite, x = "morph", y = "scars",
          add = "mean_se",
          fill = "morph", 
          palette=c('orange','white','gold'),
          ylab = "Bite scars mean", xlab = "Morph")

# Regression Model --------------------------------------------------------

model <- lm(scars ~ morph + sex + svl, data = data_bite)
summary(model)
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -4.08475    4.80514  -0.850 0.396946    
#morphw       2.80855    0.79312   3.541 0.000565 ***
#morphy      -0.81544    0.78606  -1.037 0.301614    
#sexm        -1.26796    0.69707  -1.819 0.071367 .  
#svl          0.14169    0.07419   1.910 0.058495 .  
---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  #Residual standard error: 3.432 on 122 degrees of freedom
  #Multiple R-squared: 0.1877, Adjusted R-squared: 0.1611
  #F-statistic:  7.05 on 4 and 122 DF,  p-value: 3.858e-05
  
  
  # Bite Scar Analysis ------------------------------------------------------
scar_counts <- data_bite %>% group_by(morph) %>% summarize(sums = sum(scars))
scar_counts
#"o"          170
#"w"          296
#"y"          142
# Meets chi-square assumptions.

chisq.test(c(170,296,142), p = c(1/3, 1/3, 1/3))
# X-squared = 66.408, df = 2, p-value = 3.799e-15
# There is a significant deviance in the sums of one or more groups.

pairwise.t.test(data_bite$scars, data_bite$morph, p.adj='bonferroni')
#o        w      
#w 0.0026 -      
#y 0.9562 5.6e-05