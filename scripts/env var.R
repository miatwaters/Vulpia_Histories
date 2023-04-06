
require(vegan); require(doBy)
require(lme4); require(car)
require(MuMIn);require(glmmADMB)


#### Data prep ---------------------------####

env.data<-read.csv("data/env data with site 11.csv",header=TRUE)
env.data$site<-as.factor(env.data$site)


#don't log xK,cec, ph, om, sunlight
#env.data$tot.slope<-log(env.data$tot.slope);env.data$soil.depth<-log(env.data$soil.depth);env.data$N<-log(env.data$N);   env.data$Na<-log(env.data$Na);   env.data$xNA<-log(env.data$xNA)
#env.data$P<-log(env.data$P);env.data$Ca<-log(env.data$Ca);env.data$Mg<-log(env.data$Mg);env.data$ca_mg<-log(env.data$ca_mg)
#subset(env.data,xNA>.05) #in the original data, site 25 had a value of 0.23, which is about an order of magnitude larger than the close-by sites (0.03-0.04). we have changed the data sheet to 0.023

dim(env.data)
lab<-names(env.data)

#### PCA of abiotic factors ------------------------------ ####

#vars from cca
pairs(env.data[,c(2,4:6,8:10,12,14,17:20)+1])
pc.env<-prcomp(env.data[,c(2,4:6,8:10,12,14,17:20)+1],scale.=TRUE,center=TRUE)

biplot(pc.env)

#formatted nicer:
library("factoextra")
fviz_pca_biplot <- fviz_pca_biplot(pc.env, repel = TRUE,
                                   col.var = "contrib", # Variables color
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                   col.ind = "#696969",  # Individuals color,
                                   labels = colnames(env.data$site), label.size= 8
) ;fviz_pca_biplot
#ggsave("new_biplot.png", plot = fviz_pca_biplot, width = 5, height = 5)


# Results for Variables
res.var <- get_pca_var(pc.env)
res.var$coord          # Coordinates
res.var$contrib #my favorite thing to look at - % contributions of each variable to each axis
fviz_eig(pc.env) #scree plot

summary(pc.env)
env.data$pca1<-pc.env$x[,1]
env.data$pca2<-pc.env$x[,2]

screeplot(pc.env,bstick=TRUE)

#write.csv(env.data, "env data with site 11.csv")

#### Productivity differences among sites ------------------------- ####

alt.env.data<-read.csv("data/env.means.all.old.csv",header=TRUE)
#This is another csv file that contains several additional environmental variables not relevant for the PCA above: productivity, lat, and long. Nothing else in this datasheet should be used, because the columns are not ln-transformed. 
alt.env.data$site<-as.factor(alt.env.data$site)
names(alt.env.data)
str(alt.env.data)

#Find the max and min productivity values. These represent the (dry?) biomass in grams collected from 0.75m x 0.75m plots from each Vulpia site at McLaughlin reserve in 2017.
prod_min <- min(alt.env.data$productivity, na.rm = T)
prod_max <- max(alt.env.data$productivity, na.rm = T)
#those are g/0.5625m^2....need to convert to g/1m^2 to more clearly report the biomass:
prod_min <- prod_min/0.5625;prod_min
prod_max <- prod_max/0.5625;prod_max

#which sites were they?
alt.env.data[which.min(alt.env.data$productivity),] #site 4
alt.env.data[which.max(alt.env.data$productivity),] #site 21
