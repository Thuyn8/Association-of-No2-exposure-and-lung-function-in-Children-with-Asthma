
library(reshape2)
library(ggplot2)
#library(tidyr)
library(dplyr)
library(corrplot)
library(nlme)
library(geepack)
library(multcomp)
library(table1)
library(tidyselect)
library(kableExtra)
library(ggpubr)
library(VIM)
library(mice)

# Load data
asth<- read.csv("AsthmaNO2.csv")
asth <- asth[, !(names(asth) %in% c("X","asthma"))]
asth$no2high <- as.factor(asth$no2high)

# Create age at baseline
asth$age_base <- rep(0, nrow(asth))
for (i in unique(asth$id)){
  asth$age_base[asth$id ==i]<- min(asth$age[asth$id==i])
}

# create height at baseline
asth$baseheight <- rep(0, nrow(asth))
for (i in unique(asth$id)){
  asth$baseheight[asth$id ==i]<- min(asth$height[asth$id==i])
}

# Number of follow up
Time <- asth %>%group_by(id) %>% summarise(Time = 1:n()) 
asth <- cbind(asth, Time[,2])

asth$male <- factor(asth$male, levels= c(0,1), labels= c("Female", "Male"))
asth$no2high <- factor(asth$no2high, levels = c(0,1), labels = c("Low N02", "High N02") )

#  Create table 1a:
## data frame of baseline measurements and number of follow-ups
tb1 <- data.frame(id = unique(asth$id),
                  age_base = asth$age_base[asth$Time ==1],
                  baseheight = asth$baseheight[asth$Time ==1],
                  Sex = asth$male[asth$Time == 1],
                  no2high = asth$no2high[asth$Time == 1])
## The maximum number of follow-ups for each individual
follow <-asth %>% group_by(id) %>% summarise(follow = max(Time)) %>% as.data.frame()
## factorize follow covariate
tb1$follow <- follow[,2] %>% as.factor()

## label covariates' name 
label(tb1$age_base) <- "Baseline Age"
label(tb1$baseheight) <- "Baseline Height"
label(tb1$follow) <- "No of Follow-ups"
units(tb1$age_base) <- "Year"
units(tb1$baseheight) <- "Inch"


## table 1a
t1a <-table1(~ . |no2high , data = tb1[,-1], caption = "Distribution of baseline covariates and number of follow-up")


# Create table 1b:
## reshape data from long format to wide format
asth_wide <- reshape(asth, direction = "wide",v.names = "fev1", idvar = "id", timevar = "Time" )

## log (fev) of five mesures
asth_wide[,8:12] <- asth_wide %>% dplyr::select(starts_with("fev1.")) %>% log()

## label covariates' name
label(asth_wide$fev1.1) <- "Baseline"    
label(asth_wide$fev1.2) <- "1 Year"                                       
label(asth_wide$fev1.3) <- "2 Year"
label(asth_wide$fev1.4) <- "3 Year"  
label(asth_wide$fev1.5) <- "4 Year"                                                        


## table 1
t1b <-table1(~ . |no2high, data= asth_wide[,c(5,8:12)], caption = "Summary of distriution of log(FEV) at different timepoints") 

## align two tables
kable(list(t1a, t1b)) %>% kable_styling(bootstrap_options = c("hover", "striped"),
                                        latex_options = "basic",
                                        position = "center",
                                        font_size = 8)
# Plot
## plots of log(fev) vs age, height
## plots of fev vs age, height
cols <- c("blue", "red")
ggarrange(ggplot(data = asth, aes(age, fev1)) + 
            geom_point() +
            geom_smooth() +labs(caption = "(a)"),
          ggplot(data = asth,aes(height, fev1)) + geom_point() +geom_smooth() +labs(caption = "(b)"),
          ggplot(data = asth, aes(age, log(fev1), col = no2high))  + 
            geom_point() +
            geom_smooth() +
            scale_color_manual(values = cols) + labs(caption = "(c)"),
          ggplot(data = asth,aes(height, log(fev1), col =no2high)) + 
            geom_point() +
            geom_smooth()+
            scale_color_manual(values = cols) +labs(caption = "(d)"), 
          common.legend = T) %>%
  annotate_figure(bottom = 
                    text_grob("Figure 1: 
          (a): Log FEV1 over age with LOESS smoothing curve.
(b): Same as (a), but for log FEV1 over height. 
(c): Same as (a), but stratified by average annual NO2 concentration.
(d): Same as (b), but stratified by average annual NO2 concentration."))

## plot individuals log(fev1) vs age
ggplot(data = asth, aes(age, log(fev1))) + 
  geom_line(aes(group=id, col= no2high), alpha=0.4) +
  facet_wrap(~no2high) +
  scale_color_manual(values = cols) +labs(caption = "Figure 2: Individual series of longitudinal log(FEV1) for all subjects, stratified by the high
and low NO2 groups.")

## correlation plots
par(mfrow = c(1,3))
corrplot(cor(asth[,2:4] , use = "pairwise.complete.obs"), 
         method = "square", 
         type = "upper", 
         diag = F, 
         addCoef.col = T)
corrplot(cor(asth[,2:4] %>% filter(asth$no2high == "High N02") , use = "pairwise.complete.obs"), 
         method = "square", 
         type = "upper", 
         diag = F, 
         addCoef.col = T)
corrplot(cor(asth[,2:4] %>% filter(asth$no2high =="Low N02") , use = "pairwise.complete.obs"), 
         method = "square", 
         type = "upper", 
         diag = F, 
         addCoef.col = T) 
mtext("Figure 3: Pairwise correlations among height, age, and log(FEV1).
  (left): on all available observations. 
      (center): For subjects in the high NO2 concentration area.
      (right): For subjects in the low NO2 concentration area.", side = 1,line = -4, outer =T)


## missing pattern
aggr_plot <- aggr(asth_wide, numbers = T, gap = 2,
                  cex.axis = .8,
                  labels = names(asth_wide),
                  sortVars = T,
                  ylab = c("Histogram of missing dat", "Pattern of missing data"), only.miss = T)
mtext("Figure 4: Pattern of missing data", side =1, line = 7, outer= T )

## impute data using MICE package
