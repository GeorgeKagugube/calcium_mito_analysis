## Set a global seed here for reproducibility
set.seed(101)

## Clear the workspce here 
rm(list = ls())

# Set the working directory here
setwd('/Users/gwk/Desktop/Bioinformatics/calcium_mito_analysis/data/Extracted_datasets')

## Load the modules that are needed for the analysis here 
# Load required packages
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(tidyverse)
library(ggpubr)

## Functions needed to analyse the data are written here
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

## Load the data to be analysed here
cal_dynamics <- read.csv("neuronal_calcium_dynamics.csv")

## explore the data here 
head(cal_dynamics)
tail(cal_dynamics)
str(cal_dynamics)
summary(cal_dynamics)

## Convert treatment to factor
cal_dynamics$Genotype <- as.factor(cal_dynamics$Genotype)
cal_dynamics$Stimulant <- as.factor(cal_dynamics$Stimulant)

## Organise teh levels such that WT alway appears before the Homs
cal_dynamics$Genotype <- ordered(cal_dynamics$Genotype,
                                 levels = c('WT', 'Hom'))


df2 <- data_summary(cal_dynamics, varname="Peak", 
                    groupnames=c("Genotype", "Stimulant"))

# Finished bar plot
p+labs(title="WT vs MUT", x="Genotype", y = "Counts (ct)")+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))


## Split the data to see how it looks when seperated here
cal <- cal_dynamics %>%
  filter(Stimulant == 'ATP')

## View a quick summary of the filtred data
summary(cal)

## Check for distribution and nromality here
head(cal)

# Run statistics here 
# Histogram
hist(cal$Repolarisation_Slope, breaks = 50)

## qqplot 
ggqqplot(cal$Repolarisation_Slope)

## Test for normality of dataset
# A pvalue > 0.05 == normal distribution (Parametric hypothesis test),
# A pvalue < 0.05 ==> other distribution (Non parametric hypothesis test)
shapiro.test(cal$Repolarisation_Slope)

## Add p-value comparing groups
## Specify the comparisons to be made
my_comparisons <- list(c('WT', "Hom"))

## Data visualisation using a dotplot
# Option 1
#p <- ggplot(cal, aes(x=Genotype, y=Repolarisation_Slope, colour = Stimulant)) + 
#  geom_dotplot(binaxis='y', stackdir='center') +
#  theme_classic() +
#  theme(axis.text.y = element_text(size = 15, face = 'bold'))

#p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#  stat_compare_means(label.y = 0.020, ref.group = 'WT')                   # Add global p-value

# Option 2
ggviolin(cal, x = "Genotype", y = "Repolarisation_Slope", fill = "Genotype",
         palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "jitter", add.params = list(fill = "white"),
         xlab = 'Genotype', ylab = 'Fura-2 AM Ratio (ABU/ms)')+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 0.020)                                      # Add global the p-value 

## Run a one way anova here 
# Compute the analysis of variance
res.aov <- aov(Delta_Ca ~ Genotype, data = cal)

# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

## Computing a test
t.test(Repolarisation_Slope ~ Genotype, data = cal)

## Run a wilcoxin test (Mann Whitney test)
wilc_test <- wilcox.test(Repolarisation_Slope ~ Genotype, data = cal)
wilc_test

