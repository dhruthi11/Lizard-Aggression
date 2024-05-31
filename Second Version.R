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
# Mean aggression score looks much higher for white morph
ggbarplot(data, x = "focal_morph", y = "agg_score",
          add = "mean_se",
          fill = "focal_morph", 
          palette=c('white','gold','orange'),
          ylab = "Focal morph mean aggression score", xlab = "Morph")

ggboxplot(data, x = "focal_morph", y = "agg_score", 
          add = "jitter",
          fill = "focal_morph", 
          palette=c('white','gold','orange'),
          ylab = "Focal morph mean aggression score", xlab = "Morph")

ggbarplot(data, x = "receiver_morph", y = "agg_score",
          add = "mean_se",
          fill = "receiver_morph", 
          palette=c('orange','white','gold'),
          ylab = "Average Agression Score", xlab = "Morph")

# Chi-square tests --------------------------------------------------------

#BITE:
biteYW <- c(5,11)
res.biteYW <- chisq.test(biteYW, p = c(1/2,1/2))
res.biteYW
# X-squared = 2.25, df = 1, p-value = 0.1336
#No difference in BITE in Y vs. W

#biteYO <- c(0,0)
#res.biteYO <- chisq.test(biteYO, p = c(1/2,1/2))
#res.biteYO

#biteWO<- c(5,2)
#res.biteWO <- chisq.test(biteWO, p = c(1/2,1/2))
#res.biteWO
# X-squared = 1.2857, df = 1, p-value = 0.2568


#CHASE:
chaseYW <- c(3,19)
res.chaseYW <- chisq.test(chaseYW, p = c(1/2,1/2))
res.chaseYW
# X-squared = 11.636, df = 1, p-value = 0.0006467
#Difference in CHASE in Y vs. W

#chaseYO <- c(0,0)
#res.chaseYO <- chisq.test(chaseYO, p = c(1/2,1/2))
#res.chaseYO

#chaseWO <- c(6,1)
#res.chaseWO <- chisq.test(chaseWO, p = c(1/2,1/2))
#res.chaseWO
# X-squared = 3.5714, df = 1, p-value = 0.05878


#GAPE:
#gapeYW <- c(5,4)
#res.gapeYW <- chisq.test(gapeYW, p = c(1/2,1/2))
#res.gapeYW
# X-squared = 0.11111, df = 1, p-value = 0.7389

#gapeYO <- c(1,0)
#res.gapeYO <- chisq.test(gapeYO, p = c(1/2,1/2))
#res.gapeYO
# X-squared = 1, df = 1, p-value = 0.3173

gapeWO <- c(3,7)
res.gapeWO <- chisq.test(gapeWO, p = c(1/2,1/2))
res.gapeWO
# X-squared = 1.6, df = 1, p-value = 0.2059
#No difference in GAPE in W vs. O


#LUNGE:
lungeYW <- c(2,15)
res.lungeYW <- chisq.test(lungeYW, p = c(1/2,1/2))
res.lungeYW
# X-squared = 9.9412, df = 1, p-value = 0.001616
#Difference in LUNGE in Y vs. W

#lungeYO <- c(0,0)
#res.lungeYO <- chisq.test(lungeYO, p = c(1/2,1/2))
#res.lungeYO

#lungeWO <- c(4,2)
#res.lungeWO <- chisq.test(lungeWO, p = c(1/2,1/2))
#res.lungeWO
# X-squared = 0.66667, df = 1, p-value = 0.4142


# Exact tests of goodness-of-fit ------------------------------------------

#BITE-WO
dat_bitewo <- data.frame(
  "biteo" = c(5,2),
  "bitew" = c(2,5),
  row.names = c("bite", "bitten"),
  stringsAsFactors = FALSE
)
colnames(dat_bitewo) <- c("orange", "white")
dat_bitewo

res.biteWO <- fisher.test(dat_bitewo)
res.biteWO
#p-value = 0.2861
#No difference in BITE in W vs. O


#CHASE-WO
dat_chasewo <- data.frame(
  "chaseo" = c(6,1),
  "chasew" = c(1,6),
  row.names = c("chasing", "chased"),
  stringsAsFactors = FALSE
)
colnames(dat_chasewo) <- c("orange", "white")
dat_chasewo

res.chaseWO <- fisher.test(dat_chasewo)
res.chaseWO
#p-value = 0.02914
#Difference in CHASE in W vs. O


#GAPE-YW
dat_gapeyw <- data.frame(
  "gapey" = c(5,4),
  "gapew" = c(4,5),
  row.names = c("gaping", "gaped"),
  stringsAsFactors = FALSE
)
colnames(dat_gapeyw) <- c("yellow", "white")
dat_gapeyw

res.gapeYW <- fisher.test(dat_gapeyw)
res.gapeYW
#p-value = 1
#No difference in GAPE in Y vs. W

#GAPE-YO
dat_gapeyo <- data.frame(
  "gapeo" = c(1,0),
  "gapey" = c(0,1),
  row.names = c("gaping", "gaped"),
  stringsAsFactors = FALSE
)
colnames(dat_gapeyo) <- c("orange", "yellow")
dat_gapeyo

res.gapeYO <- fisher.test(dat_gapeyo)
res.gapeYO
#p-value = 1
#No difference in GAPE in Y vs. O


#LUNGE-WO 
dat_lungewo <- data.frame(
  "lungeo" = c(4,2),
  "lungew" = c(2,4),
  row.names = c("lunging", "lunged"),
  stringsAsFactors = FALSE
)
colnames(dat_lungewo) <- c("orange", "white")
dat_lungewo

res.lungeWO <- fisher.test(dat_lungewo)
res.lungeWO
#p-value = 0.5671
#No difference in LUNGE in W vs. O


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


# Bite Scar Analysis ------------------------------------------------------
data_bite <- read.csv("data_bite.csv")
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

#pairwise.t.test(data_bite$scars, data_bite$morph, p.adj='holm')
#o        w      
#w 0.0017 -      
#y 0.3187 5.6e-05

#orange vs. white
#chisq.test(c(170,296), p=c(1/2,1/2))
# X-squared = 34.069, df = 1, p-value = 5.32e-09

#orange vs. yellow
#chisq.test(c(170,142), p=c(1/2,1/2))
# X-squared = 2.5128, df = 1, p-value = 0.1129

#white vs. yellow
#chisq.test(c(296,142), p=c(1/2,1/2))
# X-squared = 54.146, df = 1, p-value = 1.861e-13


# Modeling ---------------------------------------------------------
data_bite <- head(data_bite, 127)
plot(scars ~ svl, data = data_bite)
ggbarplot(data_bite, x = "sex", y = "scars",
          add = "mean_se", #I used the mean here, but the sum when doing the chi-square test.
          fill = "sex", 
          palette=c('cornsilk1','aliceblue'),
          ylab = "Bite scars mean", xlab = "Sex")
ggbarplot(data_bite, x = "morph", y = "scars",
          add = "mean_se",
          fill = "morph", 
          palette=c('orange','white','gold'),
          ylab = "Bite scars mean", xlab = "Morph")

model <- lm(scars ~ morph + sex + svl, data = data_bite)
summary(model)

# Chi-square test for sex -------------------------------------------------
scar_counts_sex <- data_bite %>% group_by(sex) %>% summarize(sums = sum(scars))
scar_counts_sex
# f       218
# m       390

chisq.test(c(218,390), p = c(1/2, 1/2))
#X-squared = 48.658, df = 1, p-value = 3.047e-12
#There is a significant difference in bite scars based on sex.


# Intramorph Contests -----------------------------------------------------

#Bite OW
dat_lungewo_intra <- data.frame(
  "biteo" = c(0,9),
  "bitew" = c(9,0),
  row.names = c("biting", "bitten"),
  stringsAsFactors = FALSE
)
colnames(dat_lungewo_intra) <- c("orange", "white")
dat_lungewo_intra

res.lungeWO_intra <- fisher.test(dat_lungewo_intra)
res.lungeWO_intra
#p-value = 4.114e-05
#Significant difference between WW and OO biting.

#Bite YW
dat_lungewy_intra <- data.frame(
  "bitey" = c(0,9),
  "bitew" = c(9,0),
  row.names = c("biting", "bitten"),
  stringsAsFactors = FALSE
)
colnames(dat_lungewy_intra) <- c("yellow", "white")
dat_lungewy_intra

res.lungeWY_intra <- fisher.test(dat_lungewy_intra)
res.lungeWY_intra
#p-value = 4.114e-05
#Significant difference between WW and YY biting.

#Bite OY
dat_lungeoy_intra <- data.frame(
  "bitey" = c(0,0),
  "biteo" = c(0,0),
  row.names = c("biting", "bitten"),
  stringsAsFactors = FALSE
)
colnames(dat_lungeoy_intra) <- c("yellow", "orange")
dat_lungeoy_intra

res.lungeOY_intra <- fisher.test(dat_lungeoy_intra)
res.lungeOY_intra
#p-value = 1
#No significant difference between OO and YY biting. values are 0!


#Chase OW
chaseOW_intra <- c(0,27)
res.chaseOW_intra <- chisq.test(chaseOW_intra, p = c(1/2,1/2))
res.chaseOW_intra
#X-squared = 27, df = 1, p-value = 2.035e-07
#Significant difference between WW and OO chasing.

#Chase YW
chaseYW_intra <- c(3,27)
res.chaseYW_intra <- chisq.test(chaseYW_intra, p = c(1/2,1/2))
res.chaseYW_intra
#X-squared = 19.2, df = 1, p-value = 1.177e-05
#Significant difference between WW and YY chasing.

#Chase OY
dat_chaseoy_intra <- data.frame(
  "chasey" = c(3,0),
  "chaseo" = c(0,3),
  row.names = c("chasing", "chased"),
  stringsAsFactors = FALSE
)
colnames(dat_chaseoy_intra) <- c("yellow", "orange")
dat_chaseoy_intra

res.chaseOY_intra <- fisher.test(dat_chaseoy_intra)
res.chaseOY_intra
#p-value = 0.1
#No significant difference between OO and YY chasing.


#Gape OW
dat_gapewo_intra <- data.frame(
  "gapeo" = c(9,0),
  "gapew" = c(0,9),
  row.names = c("gaping", "gaped"),
  stringsAsFactors = FALSE
)
colnames(dat_gapewo_intra) <- c("orange", "white")
dat_gapewo_intra

res.gapeWO_intra <- fisher.test(dat_gapewo_intra)
res.gapeWO_intra
#p-value = 4.114e-05
#Significant difference between OO and WW gaping.

#Gape YW
gapeWY_intra <- c(9,1)
res.gapeWY_intra <- chisq.test(gapeWY_intra, p = c(1/2,1/2))
res.gapeWY_intra
#X-squared = 6.4, df = 1, p-value = 0.01141
#Significant difference between YY and WW gaping.

#Gape OY
dat_gapeyo_intra <- data.frame(
  "gapeo" = c(1,0),
  "gapey" = c(0,1),
  row.names = c("gaping", "gaped"),
  stringsAsFactors = FALSE
)
colnames(dat_gapeyo_intra) <- c("orange", "yellow")
dat_gapeyo_intra

res.gapeYO_intra <- fisher.test(dat_gapeyo_intra)
res.gapeYO_intra
#p-value = 1
#No significant difference between OO and YY gaping.


#Lunge OW
dat_lungewo_intra <- data.frame(
  "lungeo" = c(16,0),
  "lungew" = c(0,16),
  row.names = c("lunging", "lunged"),
  stringsAsFactors = FALSE
)
colnames(dat_lungewo_intra) <- c("orange", "white")
dat_lungewo_intra

res.lungeWO_intra <- fisher.test(dat_lungewo_intra)
res.lungeWO_intra
#p-value = 3.327e-09
#Significant difference between OO and WW lunging.

#Lunge YW
dat_lungewy_intra <- data.frame(
  "lungey" = c(16,0),
  "lungew" = c(0,16),
  row.names = c("lunging", "lunged"),
  stringsAsFactors = FALSE
)
colnames(dat_lungewy_intra) <- c("yellow", "white")
dat_lungewy_intra

res.lungeWY_intra <- fisher.test(dat_lungewy_intra)
res.lungeWY_intra
#p-value = 3.327e-09
#Significant difference between YY and WW lunging.

#Lunge OY
#both are 0
#p-value = 1


# Color Morph Differences in Males and in Females -------------------------
#Females:
scar_counts_female <- data_bite %>% filter(sex == "f") %>% group_by(morph) %>% summarize(sums = sum(scars))
scar_counts_female
#"o"          33
#"w"          111
#"y"          74
# Meets chi-square assumptions.

chisq.test(c(33,111,74), p = c(1/3, 1/3, 1/3))
# X-squared = 41.899, df = 2, p-value = 7.975e-10
# There is a significant deviance in the sums of one or more groups in females.

data_bite_female <- data_bite %>% filter(sex == "f")
pairwise.t.test(data_bite_female$scars, data_bite_female$morph, p.adj='bonferroni')
   #o     w      
#w 0.0026 -      
#y 0.9562 5.6e-05


#Males:
scar_counts_male <- data_bite %>% filter(sex == "m") %>% group_by(morph) %>% summarize(sums = sum(scars))
scar_counts_male
#"o"          137
#"w"          185
#"y"          68
# Meets chi-square assumptions.

chisq.test(c(137,185,68), p = c(1/3, 1/3, 1/3))
# X-squared = 53.215, df = 2, p-value = 2.782e-12
# There is a significant deviance in the sums of one or more groups in males.

data_bite_male <- data_bite %>% filter(sex == "m")
pairwise.t.test(data_bite_male$scars, data_bite_male$morph, p.adj='bonferroni')
   #o      w      
#w 0.00503 -      
#y 0.64768 0.00024


# Updated Chi-Square/Bonferroni Tests --------------------------------------------------------

lunge <- c(2, 35, 2)
res.lunge <- chisq.test(lunge, p = c(1/3,1/3,1/3))
res.lunge
#X-squared = 55.846, df = 2, p-value = 7.467e-13

pairwise.t.test(data$lunge, data$focal_morph, p.adj='bonferroni')
#o          w      
#w 2.5e-06  -      
#y 1        6.5e-07


gape <- c(7, 16, 7)
res.gape <- chisq.test(gape, p = c(1/3,1/3,1/3))
res.gape
#X-squared = 5.4, df = 2, p-value = 0.06721


chase <- c(1, 52, 6)
res.chase <- chisq.test(chase, p = c(1/3,1/3,1/3))
res.chase
#X-squared = 80.373, df = 2, p-value < 2.2e-16

pairwise.t.test(data$chase, data$focal_morph, p.adj='bonferroni')
#o          w      
#w 4.9e-12  -      
#y 0.47     4.2e-10


bite <- c(2, 25, 5)
res.bite <- chisq.test(bite, p = c(1/3,1/3,1/3))
res.bite
#X-squared = 29.312, df = 2, p-value = 4.314e-07

pairwise.t.test(data$bite, data$focal_morph, p.adj='bonferroni')
#o         w     
#w 0.0034  -     
#y 1.0000  0.0230