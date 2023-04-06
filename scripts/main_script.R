
# loading and wrangling ####

## load required packages
pkgs <- c("nlme", "lme4", "MuMIn", "lattice", "glmmTMB",  "dplyr", "car", "visreg", "ggeffects",
          "cowplot", "dotwhisker", "tidyverse", "devtools", "ggplot2", "ggpubr", "DHARMa", "emmeans", "ggarrange", "ggrepel") 

lapply(pkgs, library, character.only = TRUE)


## set a plotting theme and load color palettes
theme_set(theme_classic())

#Color palettes
# custom colors for manuscript:
Comps.colors <- c(n = "#e78c51", y = "#349481") #orange and turquoise instead = much better complementary contrast
type.colors <- c(greenhouse = "#8dcfaa", field = "#f6be8b") #green and burnt orange, verified works for all colorblindness and greyscale

#add code here to source("env var.R") - want to run everything in that script to calculate PCA values and produce the biplot for the supplement.

outcome <- read_csv("data/persistence_project_outcome_duplicate8y_removed_newcol.csv") 
original <- read_csv("data/persistence_project_original.csv")

names(outcome)
names(original)

all <- merge(outcome, original, by = "Planting_ID")
names(all)

abiotic <- read_csv("data/env data with site 11.csv")
#This csv has all the ln-transformed variables that were used in the PCA biplot, and the pca1 and pca2 columns are correctly based on that. These pca columns are identical to the pca1.new and pca2.new columns that we added at some point to env.means.all.old.csv, which otherwise has the untransformed data in addition to just three unique variables not currently included in this csv: productivity, lat, and long.

#need to rename col 1 to "population":
colnames(abiotic)[1] = "population"
#need to add a V to every cell in this column, so it matches with the column in all:
abiotic$population<-sub("^","V",abiotic$population)

all_abiotic <- merge(all, abiotic, by = "population") #merging our dataset with the abiotic site data from 2014. THIS MERGE REMOVES THE THREE VS SITES, because they don't have abiotic data!
names(all_abiotic)
head(all_abiotic)

data1 <- all_abiotic
data1$Comps <- factor(data1$Comps, levels=c("n","y"))

data1$Block_unique <- paste(data1$Block, data1$Comps, sep = "")  

data1$Block <- as.factor(data1$Block)
data1$Block_unique <- as.factor(data1$Block_unique)
levels(data1$Block_unique)

unique(data1$population) #DOES NOT INCLUDE VS2, VS31, or VS33.
data1$population<-factor(data1$population, levels=c("V1","V3","V4","V5","V6","V8","V9","V10","V11","V13","V15","V16","V17","V19","V20","V22","V25","V26","V27","V29"))

count(data1, plant_present, sort = T) #says there are 29 instances where the tube was missing ("x")

data1_sub_present <- subset(data1, plant_present !="x") #this removes the instances where the tube was missing

#include all transplanted individuals, even if failed to germinate
data1_sub_present$seeds <- ifelse(is.na(data1_sub_present$seeds)==TRUE,0,data1_sub_present$seeds)
data1_sub_present$seeds<-ifelse(data1_sub_present$Block==1, data1_sub_present$seeds+0.5,data1_sub_present$seeds)

head(data1_sub_present)
str(data1_sub_present)
summary(subset(data1_sub_present,Block==1)) #to check if +0.5 seeds worked


# "data23" ####
#repeating these wrangling steps but for all 23 populations, without the abiotic data included. This df will only be used in analyses that did not include the environment as a predictor (e.g. figure 3, and some supplementary figures)
data23 <- all
data23$Comps <- factor(data23$Comps, levels=c("n","y"))

data23$Block_unique <- paste(data23$Block, data23$Comps, sep = "")  

data23$Block <- as.factor(data23$Block)
data23$Block_unique <- as.factor(data23$Block_unique)
levels(data23$Block_unique)

unique(data23$population) 
data23$population<-factor(data23$population, levels=c("V1","V3","V4","V5","V6","V8","V9","V10","V11","V13","V15","V16","V17","V19","V20","V22","V25","V26","V27","V29","VS2","VS31","VS33"))

count(data23, plant_present, sort = T) #says there are 29 instances where the tube was missing ("x")

data23 <- subset(data23, plant_present !="x") #this removes the instances where the tube was missing

#include all transplanted individuals, even if failed to germinate
data23$seeds <- ifelse(is.na(data23$seeds)==TRUE,0,data23$seeds)
data23$seeds<-ifelse(data23$Block==1, data23$seeds+0.5,data23$seeds)
data23<-subset(data23,type=="greenhouse")

head(data23)
str(data23)
summary(subset(data23,Block==1))

# environmental history - main analysis  ####

data1_sub_present2<-data1_sub_present %>% select(-"seed weight",-notes,-glumes,-undeveloped,-Notes,-present_with_seeds)
data_env<-na.omit(data1_sub_present2) 
data_env<-subset(data_env,type=="greenhouse")

data<-subset(data_env,population!="V6") #remove the local population


#lm<-glmmTMB(seeds ~ Comps * poly(pca1,2)  + Comps * poly(pca2,2)  + (1|Block_unique) + (1|population), data = data, family = poisson, ziformula=~1); Anova(lm)
#lm1<-glmmTMB(seeds ~ Comps * poly(pca1,2)  + Comps * poly(pca2,2)  + (1|Block_unique) + (1|population), data = data, family = poisson, ziformula=~Comps); Anova(lm1)
lm2<-glmmTMB(seeds ~ Comps * poly(pca1,2)  + Comps * poly(pca2,2)  + (1|Block_unique) + (1|population), data = data, family = poisson, ziformula=~.)
#lm3<-glmmTMB(seeds ~ Comps * poly(pca1,2)  + Comps * poly(pca2,2)  + (1|Block_unique) + (1|population), data = data, family = poisson, ziformula=~pca1+pca2); Anova(lm3)

#anova(lm,lm1,lm2,lm3)

lm<-lm2 #lm2 is the best model, so using that as "lm" from now forward.

summary(lm)

#use DHARMa to test goodness of fit
library(DHARMa)
simulationOutput<-simulateResiduals(lm)
plot(simulationOutput)
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)

Anova(lm,type = 3, component="cond")
Anova(lm,type = 3, component="zi")

##emtrends 
emt1 <- emtrends(lm, ~ Comps | pca1, "pca1", max.degree = 2)
emt1

emt2 <- emtrends(lm, ~ Comps | pca2, "pca2", max.degree = 2) 
emt2





#### figure 2 - main results ####
#Figure 2A = fe PCA1 predictions only -> reproductive fitness
#Figure 2B = fe PCA2 predictions only -> reproductive fitness
#Figure 2C = fe.zi PCA1 predictions and means+error bars -> population growth
#Figure 2D = fe.zi PCA2 predictions and means+error bars -> population growth

#fig 2B (PCA2 reproductive fitness) ####
vis<-ggpredict(lm, terms = c("pca2 [all]","Comps"), type="fe") #pca2 fe - no zeros. This is reproductive fitness only for the plants that grew and produced seeds
#basic plot:
plot1<-plot(vis) + geom_vline(xintercept=0.55167, linetype="dashed"); plot1

#formatted fig for top panel (B):
#adjust the axis limits
min(data$pca2) #-2.069
max(data$pca2) #4.236

PCA2_repro <- ggplot(vis, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymax=conf.high, ymin=conf.low), alpha = 0.2, color = NA) +
  scale_fill_manual(values=Comps.colors, guide = FALSE) +
  scale_color_manual(values=Comps.colors, name = "Community treatment", labels = c("cleared", "intact")) +
  geom_vline(xintercept=0.5517, linetype="dashed") + #this is the PCA2 value for V6 (the local site)
  labs(x = "PCA axis 2", y = element_blank()) +
  #theme(legend.position="top") 
  theme(
  legend.position = c(1, 1),
legend.justification = c("right", "top"),
legend.box.just = "right",
legend.margin = margin(6, 6, 6, 6), 
plot.margin = unit(c(1,0.5,0.5,0.1), "cm"))+ 
  scale_x_continuous(limits = c(-1.5,4.55), expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 8, limits = c(0,120), expand = c(0, 0))
PCA2_repro


#add data means 
#calculate means for plants that reproduced only
means_rm0s <- data_env %>%
  group_by(population, Comps) %>%
  filter(seeds>0.5) %>%
  mutate(mean_seeds = mean(seeds))

#calculate means for all
means_all <- data_env %>%
  group_by(population) %>%
  mutate(mean_seeds = mean(seeds, na.rm = TRUE),
         sd_seeds = sd(seeds, na.rm = TRUE))

PCA2_repro_means <- PCA2_repro + 
  geom_point(data=means_rm0s, aes(x = pca2, y = mean_seeds, color = Comps), size = 2, alpha = 1, inherit.aes = FALSE) +
  geom_point(data=subset(means_rm0s, population == "V6"), aes(x = pca2, y = mean_seeds, color = Comps), size = 4, shape = 18, inherit.aes = FALSE) +
  stat_summary(mapping = aes(x = pca2, y = seeds, color = Comps), data = means_rm0s, inherit.aes = F, fun.data = mean_se, geom = "errorbar", alpha = 0.5, width = 0.1) #SE bars 

PCA2_repro_means
#this plots the means for each population, from ONLY PLANTS THAT PRODUCED SEEDS (all zeros removed before calculating the means). Might not be a valuable visual for this reason.


  




  

#fig 2A (PCA1 reproductive fitness) ####
vis<-ggpredict(lm, terms = c("pca1 [all]","Comps"), type="fe") #pca1 fe
plot2<-plot(vis) + geom_vline(xintercept=0.4659278, linetype="dashed"); plot2

#formatted fig for top panel (A):
#adjust the axis limits
min(data$pca1) #-2.582
max(data$pca1) #4.904

#formatted: 
PCA1_repro <- ggplot(vis, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymax=conf.high, ymin=conf.low), alpha = 0.2, color = NA) +
  scale_fill_manual(values=Comps.colors, guide = FALSE) +
  scale_color_manual(values=Comps.colors, name = "Community treatment", labels = c("cleared", "intact")) +
  geom_vline(xintercept=0.4659, linetype="dashed") + #this is the PCA1 value for V6 (the local site)
  labs(x = "PCA axis 1", y = "Mean reproductive fitness \n(number of seeds)") +
  theme(legend.position="none",plot.margin = unit(c(1,0.5,0.5,1), "cm"));PCA1_repro
#scale_fill_discrete(name = "Community treatment", labels = c("cleared", "intact"))
#legend(location = "topright", inset = .05, title = "Community treatment", c("cleared", "intact"), fill = group)

PCA1_repro_means <- PCA1_repro + 
  geom_point(data=means_rm0s, aes(x = pca1, y = mean_seeds, color = Comps), size = 2, alpha = 1, inherit.aes = FALSE) +
  geom_point(data=subset(means_rm0s, population == "V6"), aes(x = pca1, y = mean_seeds, color = Comps), size = 4, shape = 18, inherit.aes = FALSE)+
  stat_summary(mapping = aes(x = pca1, y = seeds, color = Comps), data = means_rm0s, inherit.aes = F, fun.data = mean_se, geom = "errorbar", alpha = 0.5, width = 0.1) +  #foreign SE bars
  scale_x_continuous(limits = c(-2.7,5), expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 8, limits = c(0,120), expand = c(0, 0))
PCA1_repro_means
#this plots the means for each population, from ONLY PLANTS THAT PRODUCED SEEDS (all zeros removed before calculating the means). Might not be a valuable visual for this reason.



fig2A_B <- plot_grid(PCA1_repro_means, PCA2_repro_means, labels = c('A', 'B'), align = "h", hjust = .25, rel_widths = c(1, 1));fig2A_B
#ggsave("fig2A-B.png", plot = fig2A_B, width = 8, height = 4)




#for bottom panel (fig C and D), use fe.zi - includes 0s in the prediction curves, so = population growth.

#fig 2D (PCA2 pop growth) ####

vis<-ggpredict(lm, terms = c("pca2[all]", "pca1[all]","Comps"), type="fe.zi")

PCA2_all <- ggplot(subset(vis,group==-0.186), aes(x = x, y = predicted, color = facet, fill = facet)) + #group ==-0.186 is the median of PCA1
  geom_line() +
  geom_ribbon(aes(ymax=conf.high, ymin=conf.low), alpha = 0.2, color = NA) +
  scale_fill_manual(values=Comps.colors, guide = FALSE) +
  scale_color_manual(values=Comps.colors, name = "Community treatment", labels = c("cleared", "intact")) +
  geom_vline(xintercept=0.5517, linetype="dashed") + #this is the PCA2 value for V6 (the local site)
  labs(x = "PCA axis 2", y = element_blank()) +
  theme(legend.position="none", plot.margin = unit(c(0.1,0.5,0.5,0.1), "cm")); PCA2_all


#code for the means and error bars
PCA2_all_means <- PCA2_all + 
  geom_hline(yintercept=1, linetype='dotted', col = 'black') +
  stat_summary(mapping = aes(x = pca2, y = seeds, color = Comps), data = subset(data_env, population!= "V6"), inherit.aes = F, fun.data = mean_se, geom = "errorbar", alpha = 0.5, width = 0.1) + #foreign SE bars
  stat_summary(mapping = aes(x = pca2, y = seeds, color = Comps), data = subset(data_env, population!= "V6"), inherit.aes = F, fun.data = mean_se, geom = "point", alpha = 1, size = 2) + #foreign means
  stat_summary(mapping = aes(x = pca2, y = seeds, color = Comps), data = subset(data_env, population == "V6"), inherit.aes = F, fun.data = mean_se, geom = "errorbar", alpha = 0.9, width = .1) + #local SE bars
  stat_summary(mapping = aes(x = pca2, y = seeds, color = Comps), data = subset(data_env, population == "V6"), inherit.aes = F, fun.data = mean_se, geom = "point", alpha = 1, size = 4, shape = 18) + #local means
  scale_x_continuous(limits = c(-1.5,4.55), expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 10, expand = c(0.01, 0)); 
PCA2_all_means



#fig 2C (PCA1 pop growth) ####

vis<-ggpredict(lm, terms = c("pca1[all]", "pca2[all]","Comps"), type="fe.zi")

PCA1_all <- ggplot(subset(vis,group==-0.12), aes(x = x, y = predicted, color = facet, fill = facet)) + #group ==-0.12 is the median of PCA2
  geom_line() +
  geom_ribbon(aes(ymax=conf.high, ymin=conf.low), alpha = 0.2, color = NA) +
  scale_fill_manual(values=Comps.colors, guide = FALSE) +
  scale_color_manual(values=Comps.colors, name = "Community treatment", labels = c("cleared", "intact")) +
  geom_vline(xintercept=0.4659, linetype="dashed") + #this is the PCA1 value for V6 (the local site)
  labs(x = "PCA axis 1", y = "Per capita population growth \n(lambda)") +
  theme(legend.position="none",
        plot.margin = unit(c(0.1,0.5,0.5,1), "cm")); PCA1_all


#code for the means and error bars
PCA1_all_means <- PCA1_all + 
  geom_hline(yintercept=1, linetype='dotted', col = 'black') +
  stat_summary(mapping = aes(x = pca1, y = seeds, color = Comps), data = subset(data_env, population!= "V6"), inherit.aes = F, fun.data = mean_se, geom = "errorbar", alpha = 0.5, width = 0.1) + #foreign SE bars
  stat_summary(mapping = aes(x = pca1, y = seeds, color = Comps), data = subset(data_env, population!= "V6"), inherit.aes = F, fun.data = mean_se, geom = "point", alpha = 1, size = 2) + #foreign means
  stat_summary(mapping = aes(x = pca1, y = seeds, color = Comps), data = subset(data_env, population == "V6"), inherit.aes = F, fun.data = mean_se, geom = "errorbar", alpha = 0.9, width = .1) + #local SE bars
  stat_summary(mapping = aes(x = pca1, y = seeds, color = Comps), data = subset(data_env, population == "V6"), inherit.aes = F, fun.data = mean_se, geom = "point", alpha = 1, size = 4, shape = 18)+  #local means 
  scale_x_continuous(limits = c(-2.7,5), expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 10, expand = c(0.01, 0))
PCA1_all_means

fig2C_D <- plot_grid(PCA1_all_means, PCA2_all_means, align = "h", rel_widths = c(1, 1));fig2C_D
#ggsave("fig2C-D.png", plot = fig2C_D, width = 8, height = 4)


#WHOLE THING!!
#fig2 <- plot_grid(fig2A_B, fig2C_D, nrow = 2);fig2
#ggsave("fig2.png", plot = fig2, width = 14, height = 8)

fig2_labeled <- plot_grid(PCA1_repro_means, PCA2_repro_means, PCA1_all_means, PCA2_all_means, nrow = 2, align = "h", rel_widths = c(1, 1), labels = "AUTO", label_x = c(.15, .01, .15,.01));fig2_labeled

#, label_x = c(-.035, -.035, -.035,-.035), label_y = c(1.05, 1.05, 1.05, 1.05)

#ggsave("fig2_labeled.png", plot = fig2_labeled, width = 8, height = 5)









#figure 3 - rank/abundance and evenness ####

df<-data23 %>% group_by(population, Comps) %>%
  mutate(rank = as.integer(rank(desc(seeds), ties.method = "first"))) 

V6 <- df %>%
  filter(population == "V6")

examples <- df %>%
  filter(population == c("V1", "V13", "V19"))


#fig 3A - rank abundance style plot ####
ra_indivs <- ggplot(df, aes(x = rank, y = seeds, color = Comps, linetype=type,group=interaction(population, Comps, type))) + 
  geom_line(size = 1, alpha = .7) + 
  labs(x = "Individual (N=10 per population per treatment)", y = "Individual fitness (number of seeds)", color = "Community treatment") + 
  geom_line(data=V6, 
             aes(x=rank,y=seeds), 
             color='black',
             size=1, alpha = .5,
            linetype = "dotted") + #add a dashed/slightly colored line for V6 in both treatments. 
  theme_classic() +
  #theme(aspect.ratio = 0.80) +
  theme(legend.position="none") +
  #geom_text(aes(label=ifelse(seeds>55,as.character(population),'')),hjust=-.3,vjust=0, size = 2.5, color = "gray30") +
  #geom_text(aes(label=ifelse(population==c("V1", "V13") & Comps== c("n"),as.character(population),'')),hjust=-.6,vjust=0, size = 2.5, check_overlap = T) +
  scale_color_manual(labels = c("cleared", "intact"), values=Comps.colors) +
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(-1, 120)) +
  scale_x_continuous(breaks = seq(1, 10, by = 1), expand = c(0, 0))
ra_indivs    

#NOTE - this includes the three VS populations that were excluded from the lm because they didn't have any associated abiotic data... and two of them have extreme values shown here (VS33 and VS2).

#ggsave("ra_indivs.png", plot = ra_indivs, width = 5, height = 5)
#ggsave("ra_indivs_V6dot.png", plot = ra_indivs, width = 5, height = 5)



#fig 3B - evenness plot ####

summary<-data23 %>% group_by(population, Comps) %>%
  summarise(n=n(),sum=sum(seeds)) 

#plotting as function of mean fitness
even <- data23 %>% left_join(summary) %>%
  mutate(prop=(seeds+1)/(sum+n)) %>%
  mutate(intermediate=prop*log(prop)) %>%
  group_by(population,Comps) %>%
  summarize(n=mean(n),sum=(sum(intermediate)*-1), mean_seeds=mean(seeds)) %>%
  mutate(evenness=sum/log(n)) %>%
  
  ggplot(aes(x=mean_seeds,y=evenness, color=Comps, group=population)) +
  geom_line(color="black", alpha = 0.2) +
  geom_point(size=3, alpha = .7) +
  #geom_jitter(width = 0.5, height = 0.5) +
  scale_color_manual(labels = c("cleared", "intact"), values=Comps.colors) +
  #geom_text(aes(label=ifelse(evenness<1,as.character(population),'')),hjust=1.6,vjust=.5, size = 2, check_overlap = T) +
  #geom_text_repel(aes(label=ifelse(evenness<1,as.character(population),'')), size = 3, check_overlap = F) +
  #geom_text_repel(aes(label=ifelse(population == "V6",as.character(population),'')), size = 3, color = "black", check_overlap = T) +
  #geom_text(aes(label=ifelse(population==c("V1", "V13", "V19"),as.character(population),'')),hjust=-.6,vjust=0, size = 2.5, check_overlap = T, color = "gray30") +
  labs(x = "Mean Individual fitness", y = "Evenness", color = "Community treatment") +
  scale_x_continuous(limits= c(0,20), expand = c(0, 0)) +
  theme(legend.position = c(1, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(0, 0, 0, 0),);even


#ggsave("evenness.png", plot = even, width = 5, height = 5)

fig3 <- plot_grid(ra_indivs, even, align = "v", labels = c('A', 'B'));fig3

fig3_labels <- plot_grid(ra_indivs, even, align = "v", labels = c('A', 'B'));fig3_labels
#ggsave("fig3.png", plot = fig3, width = 8, height = 4)
















#SI figure: zi.prob visuals ####
# zi.prob is just the zero-inflated probs --> SI figures

vis<-ggpredict(lm, terms = c("pca1 [all]","Comps"), type="zi_prob") #pca1 zi_prob
plot1<-plot(vis) + geom_vline(xintercept=0.4659278, linetype="dashed")
plot1

PCA1_zi <- ggplot(vis, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymax=conf.high, ymin=conf.low), alpha = 0.2, color = NA) +
  scale_fill_manual(values=Comps.colors, guide = FALSE) +
  scale_color_manual(values=Comps.colors, name = "Community treatment", labels = c("cleared", "intact")) +
  geom_vline(xintercept=0.4659, linetype="dashed") + #this is the PCA1 value for V6 (the local site)
  labs(x = "PCA axis 1", y = "Predicted zero-inflated probability") +
  theme(legend.position="none",plot.margin = unit(c(1,0.5,0.5,1), "cm"));PCA1_zi

vis<-ggpredict(lm, terms = c("pca2 [all]","Comps"), type="zi_prob") #pca2 zi_prob
plot2<-plot(vis) + geom_vline(xintercept=0.55167, linetype="dashed")
plot2

PCA2_zi <- ggplot(vis, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymax=conf.high, ymin=conf.low), alpha = 0.2, color = NA) +
  scale_fill_manual(values=Comps.colors, guide = FALSE) +
  scale_color_manual(values=Comps.colors, name = "Community \ntreatment", labels = c("cleared", "intact")) +
  geom_vline(xintercept=0.5517, linetype="dashed") + #this is the PCA2 value for V6 (the local site)
  labs(x = "PCA axis 2", y = element_blank()) +
  #theme(legend.position="top") 
  theme(
    legend.position = "right")+ 
  scale_x_continuous(limits = c(-1.5,4.55), expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 8, expand = c(0, 0))
PCA2_zi



zi_plots <- plot_grid(
  PCA1_zi + theme(legend.position = "none"), 
  PCA2_zi + theme(legend.position = "none"), 
  align = "vh", labels = c('A', 'B'),
  hjust = -2
  );zi_plots

legend <- get_legend(
  PCA2_zi + theme(legend.box.margin = margin(0,0,0,12))
)

zi_plots <- plot_grid(zi_plots, legend, rel_widths = c(3, .4)); zi_plots
#ggsave("zi_plots_formatted.png", plot = zi_plots, width = 8, height = 4)













#Seeds per pop boxplot: Useful SI visual ####
vis<-ggpredict(lm, terms = c("population","Comps"), type="re") #predicted counts of seeds per pop, colored by comps.
plot1<-plot(vis); plot1

#with actual data:
#seeds per population 
#options for which pops to include:
  #(data = only the 19 pops with abiotic data; included in the model)
  #(data1_sub_present = all 23 pops used in my experiment)
seeds_pop <- ggplot(data = data23, aes(x= population, y = seeds, color = Comps)) +
  geom_boxplot() +
  scale_fill_manual(values=Comps.colors, guide = "none") +
  scale_color_manual(values=Comps.colors, name = "Community treatment", labels = c("cleared", "intact")) +
  labs(x = "Population", y = "Number of seeds produced \nper individual planted") +
  scale_y_continuous(n.breaks = 10)+
  #theme(legend.position = c(1, 1), legend.justification = c("right", "top"), legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))
  theme(legend.position = "top"); seeds_pop

#ggsave("seeds_pop.png", plot = seeds_pop, width = 8, height = 5)









# 3D graph for SI ---------------------####
#need to look at each environment separately to get 3d graphs

data.n<-subset(data_env, population!="V6" & type=="greenhouse" & Comps == "n")
lm.n<-glmmTMB(seeds ~  poly(pca1,2)  + poly(pca2,2)  + (1|Block_unique) + (1|population), data = data.n, family = poisson, ziformula = ~.); Anova(lm.n)
summary(lm.n,type=3)
visreg2d(lm.n,x="pca2",y="pca1",plot.type="persp", scale="response") #surface for the cleared treatment (Comps = n)


vis1<-ggpredict(lm.n, terms = c("pca2[all]"), type="fe")
plot(vis1) 
vis1<-ggpredict(lm.n, terms = c("pca2[all]"), type="zi_prob")
plot(vis1) 

#best from high water and high nutrients, though not significantly so
data.y<-subset(data_env, population!="V6" & type=="greenhouse" & Comps == "y")
lm.y<-glmmTMB(seeds ~  poly(pca1,2)  + poly(pca2,2)  + (1|Block_unique) + (1|population), data = data.y, family = poisson, ziformula = ~.); Anova(lm.y)
summary(lm.y,type=3)
visreg2d(lm.y,x="pca2",y="pca1",plot.type="persp", scale="response") #surface for the intact treatment (Comps = y)


vis1<-ggpredict(lm.y, terms = c("pca2[all]"), type="fe")
plot(vis1) 
vis1<-ggpredict(lm.y, terms = c("pca2[all]"), type="zi_prob")
plot(vis1) 
















# pie charts for SI ####
#pie charts of no germination, germination but no seed attempt, seed attempt but none viable, successful seed attempt - in four conditions

#all planted data = all
#split into Comps = y/n dfs
comps_y <- all %>% 
  filter(Comps =='y')

comps_n <- all %>% 
  filter(Comps =='n')
#want total number of rows in each = total # planted
total_compsy <- nrow(comps_y);total_compsy #459 rows total (I lost one individual while planting)
total_compsn <- nrow(comps_n);total_compsn #460 rows total

#want number that did not sprout; plant_present = n
nosprout_compsy <- sum(comps_y$plant_present == "n") #340
nosprout_compsn <- sum(comps_n$plant_present == "n")

#want number that were lost; plant_present = x
lost_compsy <- sum(comps_y$plant_present == "x") #5
lost_compsn <- sum(comps_n$plant_present == "x")
#want number sprouted, no seeds; present_with_seeds = NA
comps_y[is.na(comps_y)] <- 0 #note - this replaced all NAs in entire df with 0. 
comps_n[is.na(comps_n)] <- 0
#sum for all rows where plant_present = y and present_with_seeds = 0....
sprout1 <- sum(comps_y$plant_present == "y" & comps_y$present_with_seeds == "0"); sprout1

sprout2<- sum(comps_y$plant_present == "y" & comps_y$present_with_seeds == "n");sprout2

sprout_noseeds_compsy <- sprout1 + sprout2; sprout_noseeds_compsy #this for Compsy

sprout3 <- sum(comps_n$plant_present == "y" & comps_n$present_with_seeds == "0"); sprout3

sprout4<- sum(comps_n$plant_present == "y" & comps_n$present_with_seeds == "n");sprout4

sprout_noseeds_compsn <- sprout3 + sprout4; sprout_noseeds_compsn #this for Compsn

#want number sprouted, with seeds; present_with_seeds = y
withseeds_compsy <- sum(comps_y$present_with_seeds == "y") #101
withseeds_compsn <- sum(comps_n$present_with_seeds == "y")

#Convert all of those into proportions of the total # planted
group <-  c("Did not grow", "Lost", "Present with no seeds", "Present with seeds")
value <-  c(nosprout_compsy, lost_compsy, sprout_noseeds_compsy, withseeds_compsy)
comps_y_Ns <- data.frame(group, value
  ) #total numbers of each, in the intact community treatment

Group <-  c("Did not grow", "Lost", "Present with no seeds", "Present with seeds")
Value <-  c(nosprout_compsn, lost_compsn, sprout_noseeds_compsn, withseeds_compsn)
comps_n_Ns <- data.frame(Group, Value
) #total numbers of each, in the cleared community treatment


#then plot pies for Comps =y/n separately
#comps = y; intact community treatment:

pie_y <- ggplot(comps_y_Ns, aes(x = "", y = value, fill = group)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = value),
            position = position_stack(vjust = 0.5), show.legend = F) +
  theme_void() +
  scale_fill_brewer(palette = "Accent", name = "Planted in intact \ncommunity")
pie_y


pie_n <- ggplot(comps_n_Ns, aes(x = "", y = Value, fill = Group)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = Value),
                   position = position_stack(vjust = 0.5), show.legend = F) +
  theme_void() +
  scale_fill_brewer(palette = "Accent", name = "")
pie_n

#plots together
pie_plots <- plot_grid(
  pie_n + theme(legend.position = "none"), 
  pie_y + theme(legend.position = "none"), 
  align = "vh", labels = c('A) Cleared community', 'B) Intact community'),
  hjust = -0.1
);pie_plots


#ggsave("pie_plots.png", plot = pie_plots, width = 6, height = 4)


#then plot bar charts for each population - UNFINISHED

Compsn_pops <- data23 %>% group_by(population) %>%
  subset(Comps=="n")

Compsy_pops <- data23 %>% group_by(population) %>%
  subset(Comps=="y")

#continue working on this later










# map of McLaughlin and my sites ####
library(sf)
mclaughlin_shapefile <- st_read("data/UC_Natural_Reserve_System_Boundaries/UC_Natural_Reserve_System_Boundaries.shp") 

McL_map <- ggplot() +
  geom_sf(data = subset(mclaughlin_shapefile, objectid == "21"),
               fill = "seashell2", color = "black", alpha = 0.3) +
  coord_sf() + #plots the outline of the reserve on a white background! How to get more map details?
  theme(
    panel.background = element_rect(fill = 'transparent')
  )
McL_map
##add my sites:
#load datasheet with lat/long for each site
env_data_2 <- read_csv("data/env.means.all.old.csv")
#create new df with just lat and long
sites <- subset(env_data_2, select=c("population", "lat", "long", "pca1.new", "pca2.new"))
#add sites to the map
library(RColorBrewer)
map_sites_pca2 <- McL_map + geom_point(data=sites, aes(x = long, y = lat, color = pca2.new), size = 3, shape = 16, alpha = .9, inherit.aes = FALSE) + scale_color_viridis_c(option = "plasma") +
  labs(x = "Longitude", y = "Latitude", color = "PCA 2") +
  geom_text_repel(data=subset(sites, population == "V6"), aes(long, lat, label = "garden"), nudge_x = -.005, nudge_y = .005, size = 3) +
  geom_text_repel(data = sites, aes(x= long, y = lat, label=ifelse(lat<40,as.character(population),'')),  size = 3, color = "gray30", check_overlap = F ) +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  );map_sites_pca2 

#ggsave("map_sites_pca2.png", plot = map_sites_pca2, width = 8, height = 7)

map_sites_pca1 <- McL_map + geom_point(data=sites, aes(x = long, y = lat, color = pca1.new), size = 3, shape = 16, alpha = .9, inherit.aes = FALSE) + scale_color_viridis_c(option = "plasma") +
  labs(x = "Longitude", y = "Latitude", color = "PCA 1") +
  geom_text_repel(data=subset(sites, population == "V6"), aes(long, lat, label = "garden"), nudge_x = -.005, nudge_y = .005, size = 3) +
  geom_text_repel(data = sites, aes(x= long, y = lat, label=ifelse(lat<40,as.character(population),'')),  size = 3, color = "gray30", check_overlap = F ) +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  );map_sites_pca1

#ggsave("map_sites_pca1.png", plot = map_sites_pca1, width = 8, height = 7)



#was able to separately get a geographic map this way:
mclaughlin_geo <- get_stamenmap(bbox = c(left = -122.46, bottom = 38.80, 
                                  right = -122.30, top = 38.9), 
                         zoom = 11)

geo <- ggmap(mclaughlin_geo);geo #not the nicest map, but fine

#can't yet figure out how to add this as a layer though:
geo + 
  geom_sf(data = subset(mclaughlin_shapefile, objectid == "21"),
             fill = "seashell2", color = "black", alpha = 0.3) +
  coord_sf()


