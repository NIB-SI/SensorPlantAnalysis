#######################################################################################
# EXPERIMENT 2 
# RYWALNAHG
# exposure of roGFP sensor lines to PVY virus
# mixed effect modeling analysis of confocal images, originally analyzed with
# https://github.com/NIB-SI/SensorPlantAnalysis 
# 
# Author: Anze Zupanic
# Date 7.12.2020
#
# LOADING LIBRARIES AND CHANGING SOME SETTINGS
# order of loading llibraries is IMPORTANT
library(readxl) # reading excel files
library(tidyverse) # data management and plotting (stringr, tidyr, ggplot, dplyr)
library(lmerTest) # mixed models package
library(emmeans) # statistics of mixed models
library(rmarkdown)


# defining pallette for plotting
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#####################################################################################
# READ the DATA
mydata <- read_excel("C:/projekti/projekti NIB/lezije/CompleteAnalysis/Experiment_4/MixedEffectsModeling/Rywal_normalizedratio_exp4NT.xlsx")

#####################################################################################
# QUESTION 1: is there a difference in signal between ROI1, ROI2 measurements.
# No figure, only statistical analysis (paired samples)

## MIXED EFFECTS, with leaves and plants as potential random effects
mydata2 = mydata %>% dplyr::filter(Treatment == 'PVY-WILGA') %>% 
dplyr::select(-c(Filename, Experiment, Genotype, Treatment)) %>%
  dplyr::filter(ROI == "ROI1" | ROI == "ROI2")

# Analyze the most complex model 
final_model = lmer(Measurement ~ dpi*ROI*Line +(1|Plant)+(1|Leaf) , data = mydata2)
anova(final_model) # Type III test, Satterthwaite's degrees of freedom method used
rand(final_model) # significance for random effects

# Confidence intervals and posthoc pairwise comparison within factors with 
# emmeans package: Kendall-Roger degrees of freedom used
# only factors that were deemed significant can be used
# Significance analysis for ROI
emmeans(final_model, pairwise~ROI)
# Significance analysis for Line
emmeans(final_model, pairwise~Line)

# conditional analysis for dpi 
emmeans(final_model, pairwise~ROI|dpi) 

# conditional analysis for Line (interpret with caution, interaction not significant) 
emmeans(final_model, pairwise~ROI|Line)

#####################################################################################
# QUESTION 2: is there a difference in signal between ROI1, ROI2 and CTR measurements.
# No figure, only statistical analysis.

## MIXED EFFECTS, with leaves and plants as potential random effects
mydata2 = mydata %>% dplyr::filter(Treatment == 'PVY-WILGA') %>% 
  dplyr::select(-c(Filename, Experiment, Genotype, Treatment)) %>%
  dplyr::filter(dpi == "d7")

# again issues with singularity, the most complex model used instead of moedl selection
final_model = lmer(Measurement ~ Line*ROI +(1|Plant) + (1|Leaf) , data = mydata2)
anova(final_model) # Type III test, Satterthwaite's degrees of freedom method used
rand(final_model) # significance of random effets cannot be calculated, too few values

# Confidence intervals and posthoc pairwise comparison within factors with 
# emmeans package: Kendall-Roger degrees of freedom used
# Significance of ROI
emmeans(final_model, pairwise~ROI)
# conditional analysis for Line (interpret with caution, interaction not significant)
emmeans(final_model, pairwise~ROI|Line)


#####################################################################################
# QUESTION 3: is there a difference in signal in MOCK plants and plants infected with PVY

# code for making the combined figure of all data
mydata2 = dplyr::filter(mydata, Treatment != 'H2O2')
mydata2 = dplyr::filter(mydata2, Treatment != 'DTT')
mydata2 = tidyr::unite(mydata2, ComplexTreatment, c("Treatment", "ROI"))
mydata2 = dplyr::select(mydata2,-c(Filename, Experiment, Genotype, Plant, Leaf))
# Here variable ComplexTreatment has four values for four different "treatments": 
# ROI1 and ROI2 and PVY-control and MOCK
# PLOT THE ENTIRE DATASET, FACETED
# the order in the plot is ROI1, ROI2, CTR, MOCK
# some labels are given new names
dpiLabels <- c(d3 = "3 dpi", d5 = "5 dpi", d7 = "7 dpi")
g = ggplot(mydata2,aes(x = factor(ComplexTreatment, 
                              level = c('PVY-WILGA_ROI1', 'PVY-WILGA_ROI2', 
                                        'PVY-WILGA_CTR','MOCK_NA'),
                              labels = c("ROI1", "ROI2", "CTR","MOCK")), 
                   y = Measurement, fill=factor(ComplexTreatment), colour=factor(ComplexTreatment))) +
  geom_boxplot(aes(fill = ComplexTreatment), alpha = 0.5, col = "black") +
  geom_point() +
  facet_grid(rows = vars(Line), cols = vars(dpi), labeller = labeller(dpi = dpiLabels)) +
  theme_classic() +
  geom_line(aes(group = Lesion), size=1, alpha=0.4, colour = "grey") +
  scale_color_manual(values = cbPalette, labels = c("MOCK", "CTR","ROI1", "ROI2")) +
  scale_fill_manual(values=cbPalette, labels = c("MOCK", "CTR","ROI1", "ROI2")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(strip.text.x = element_text(size=8, face="bold"),
        strip.text.y = element_text(size=8, face="bold")) +
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  theme(axis.title.x = element_blank()) +
  ylab("405/488 ratio") +
  scale_x_discrete(labels=c("PVY-WILGA_ROI1" = "ROI1", "PVY-WILGA_ROI2" = "ROI2",
                            "PVY-WILGA-control_CTR" = "CTR","MOCK_NA" = "MOCK"))
# saving the image in high res (you need first to name the image g)
# ggsave("Experiment4_complete.tiff", g, width=14, height=15, units="cm", dpi=300)

## MIXED EFFECTS, with leaves and plants as potential random effects
# prepare data
mydata2 = dplyr::filter(mydata, Treatment != 'H2O2')
mydata2 = dplyr::filter(mydata2, Treatment != 'DTT')
mydata2 = dplyr::filter(mydata2, dpi != 'd5')
mydata2 = dplyr::filter(mydata2, is.na(ROI) | ROI == c('ROI1','ROI2'))
mydata2 = tidyr::unite(mydata2, ComplexTreatment, c("Treatment", "ROI"))
mydata2 = dplyr::select(mydata2,-c(Filename, Experiment, Genotype))
# pick the best mixed effects model, starting with most complex design 
# (there are to few datapoints to look for interactions between fixed and random effects)
# model selection does not work, because of singularity, so most complex model used for analysis
final_model = lmer(Measurement ~ Line*dpi*ComplexTreatment  + (1|Leaf) + (1|Plant) , data = mydata2)
anova(final_model) # Type III test, Satterthwaite's degrees of freedom method used
rand(final_model) # significance for random effects

# Confidence intervals and posthoc pairwise comparison within factors with 
# emmeans package: Kendall-Roger degrees of freedom used
# only factors that were deemed significant can be used
# Significance for Complextreatment
emmeans(final_model, pairwise~ComplexTreatment)
# Significance for Line
emmeans(final_model, pairwise~Line)

# conditional analysis for dpi(interpret with caution as the interaction is not significant)
emmeans(final_model, pairwise~ComplexTreatment|dpi)
# conditional analysis for Line (interpret with caution as the interaction is not significant)
emmeans(final_model, pairwise~ComplexTreatment|Line)


######################################################################################
# QUESTION 4: can we confirm that the control DDT, H2O2 has an effect on the plants #

# plotting the data
# prepare data
mydata2 = dplyr::filter(mydata, Treatment %in% c("DTT", "ddH2O",'H2O2'))
mydata2 = dplyr::select(mydata2,-c(Filename, Experiment, Genotype))
# variable Treatment has three values, according to the treatment used: 
# DTT and ddH20 and H2O2
# plotting
cbPalette2 <- c("#0072B2", "#D55E00", "#CC79A7")
dpiLabels <- c(d3 = "3 dpi", d5 = "5 dpi", d7 = "7 dpi", d4 = "4 dpi")
g = ggplot(mydata2,aes(x=Treatment, y=Measurement, fill=factor(Treatment), colour=factor(Treatment))) +
  geom_boxplot(aes(fill = Treatment), alpha = 0.5, col = "black") +
  geom_point() +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  facet_grid(rows = vars(Line), cols = vars(dpi), labeller = labeller(dpi = dpiLabels)) +
  scale_fill_manual(values=cbPalette2)+ # Boxplot fill color
  scale_color_manual(values = cbPalette2) + # Jitter color palette
  theme_classic() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(strip.text.x = element_text(size=8, face="bold"),
        strip.text.y = element_text(size=8, face="bold")) +
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 65, vjust=0.6)) +
  theme(axis.title.x = element_blank()) +
  ylab("405/488 ratio")
# saving the image in high res (you need first to name the image g)
# ggsave("Experiment4_controls.tiff", g, width=10, height=10.5, units="cm", dpi=300)

## MIXED EFFECTS, with leaves and plants as potential random effects
# pick the best mixed effects model, starting with most complex design 
# (there are to few datapoints to look for interactions between fixed and random effects)
final_model = lmer(Measurement ~ Line*Treatment  + (1|Plant) + (1|Leaf) , data = mydata2)
anova(final_model) # Type III test, Satterthwaite's degrees of freedom method used
rand(final_model) # significance for random effects

# Confidence intervals and posthoc pairwise comparison within factors with 
# emmeans package: Kendall-Roger degrees of freedom used
# only factors that were deemed significant can be used
# If non-significant factors need to be analyzed, then they need to be included 
# into the final_model manually
# Significance analysis for treatment
emmeans(final_model, pairwise~Treatment)
# conditional analysis for Line (interpret with caution since interaction is not significant)
emmeans(final_model, pairwise~Treatment|Line)










