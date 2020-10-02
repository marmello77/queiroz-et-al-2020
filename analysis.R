################################################################################
################################################################################
#Supplement to the paper Queiroz et al. (2020, Biotropica).

#Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com.

#Authors: Joel A. Queiroz, Ugo M. Diniz, Diego P. Vázquez, Zelma M. Quirino,
#Francisco A.R. Santos, Marco A.R. Mello, Isabel C. Machado.

#See README for further info:
#https://github.com/marmello77/queiroz_et_al_2020/blob/main/README.md
################################################################################
################################################################################



################################################################################
##### SET THE STAGE
################################################################################


#Set the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Delete all previous objects
rm(list= ls())

#Load the required packages
library(igraph)
library(bipartite)
library(Rmisc)
library(vegan)
library(gdata)
library(ggplot2)

#Load some custom-made functions
source("RestNullModel.R")
source("PosteriorProb.R")
source("MyTriangle.R")



################################################################################
##### PROCESS THE NETWORK
################################################################################


#Import the network
data <- as.matrix(read.delim("data/network.txt", row.names=1))

#Inspect the network
class(data)
data
dim(data)
min(data)
max(data)

#Plot the matrix
visweb(data)

#Convert the network to igraph format
data2 <- graph_from_incidence_matrix(data, directed = F, weighted = TRUE) 

#Inspect object
class(data2)
data2
E(data2)
V(data2)$name

#Inform which nodes represent which taxonomic groups
V(data2)$set[1:nrow(data)] = c("Moths",	"Moths", "Bats",	"Bats",	"Moths",	"Moths",
                       "Moths",	"Moths",	"Moths",	"Bats",	"Bats",	"Moths",
                       "Moths",	"Moths",	"Moths",	"Moths",	"Moths",
                       "Moths",	"Moths")
V(data2)$set[(nrow(data)+1):(nrow(data)+ncol(data))] = "Plants"



################################################################################
##### DRAW THE NETWORK
################################################################################


#Set layout
lay1 <- layout_nicely(data2)

#Set edge curvatures
curves1 = curve_multiple(data2)

#Set edge mode and width
E(data2)$arrow.mode = 0
E(data2)$width = E(data2)$weight/5+1

#Calculate Louvain modularity (resolution = 1.0, similar to DIRT_LPA+)
data2.lou = cluster_louvain(data2)

#Import "diamond" vertex shape
source("MyTriangle.R")

#Set vertex shapes
V(data2)$shape = V(data2)$set
V(data2)$shape = gsub("Bats","diamond",V(data2)$shape)
V(data2)$shape = gsub("Moths","square",V(data2)$shape)
V(data2)$shape = gsub("Plants","circle",V(data2)$shape)

##Set node and cloud colors by modularity
colrs <- rainbow(length(data2.lou), alpha = 1.0, s = 1, v = 0.8)
V(data2)$color <- colrs[data2.lou$membership]
clouds = rainbow(length(data2.lou), alpha = 0.1)

#Plot and export the graph
tiff(filename= "figures/network.tif", res= 300, height= 3000, width= 3100)
par(mfrow=c(1,1),mar=c(1,1,1,5))
plot(data2.lou,
     data2,
     col = V(data2)$color,
     mark.border="lightgrey", 
     mark.col=clouds,
     vertex.size=7.5,
     vertex.label=V(data2)$name,
     vertex.label.color="white",
     vertex.label.cex=.3,
     edge.color = adjustcolor("grey", alpha.f = .5),
     edge.curved=0.3,
     edge.width = 3,
     layout=lay1)
legend(x = 0.9,y = 1.0, legend = c("Bats", "Moths", "Plants"),
       pch = c(18,15,19),  title="Taxon",
       text.col = "gray20", title.col = "black",
       box.lwd = 0, cex = 2, col=c("grey", "grey", "grey"))
par(mfrow=c(1,1))
dev.off()



################################################################################
##### NETWORK LEVEL ANALYSIS
################################################################################


#Set seed
set.seed(14)

#Set the number of permutations to be used in all null model analyses
permutations <- 1000

#Generate randomized matrices
nulls <- nullmodel(data, N=permutations, method="vaznull")


##### MODULARITY

#Calculate metric for the original network
Mod <- computeModules(data, method = "Beckett")
Mod@likelihood

#Extract module membership
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(data)]
col.Part <- Part[(nrow(data)+1):(nrow(data)+ncol(data))]

#Calculate metric for the randomized networks
nullmod <- sapply(nulls, computeModules, method = "Beckett")
modnull <- sapply(nullmod, function(x) x@likelihood)
(Mod@likelihood - mean(modnull))/sd(modnull) # Z value
Mod.sig <- sum(modnull>(Mod@likelihood)) / length(modnull) # p value
Mod.sig

#Plot the observed value against the distribution of randomized values
plot(density(modnull), main="Observed vs. randomized",
     xlim=c(min((Mod@likelihood), min(modnull)), 
            max((Mod@likelihood), max(modnull))))
abline(v=Mod@likelihood, col="red", lwd=2, xlab="")


Mod@likelihood #observed
mean(modnull) #randomized mean
sd(modnull) #randomized SD
(Mod@likelihood - mean(modnull))/sd(modnull) # Z-value
sum(modnull>(Mod@likelihood)) / length(modnull) #randomized > observed
sum(modnull<(Mod@likelihood)) / length(modnull) #randomized < observed


##### SPECIALIZATION 

#Calculate metric for the original network
Spec <- networklevel(data, index="H2")
class(Spec)

#Calculate metric for the randomized networks
randomized.Spec <- unlist(sapply(nulls, networklevel, index="H2"))
(Spec - mean(randomized.Spec))/sd(randomized.Spec) # Z value
Spec.sig <- sum(randomized.Spec>Spec)/length(randomized.Spec) # p value
Spec.sig

#Plot the observed value against the distribution of randomized values
plot(density(randomized.Spec), main="Observed vs. randomized",
     xlim=c(min((Spec), min(randomized.Spec)), 
            max((Spec), max(randomized.Spec))))
abline(v=Spec, col="red", lwd=2, xlab="")

Spec #observed
mean(randomized.Spec) #randomized mean
sd(randomized.Spec) #randomized SD
(Spec - mean(randomized.Spec))/sd(randomized.Spec) # Z-value
sum(randomized.Spec>(Spec)) / length(randomized.Spec) #randomized > observed
sum(randomized.Spec<(Spec)) / length(randomized.Spec) #randomized < observed


##### NESTEDNESS

#Calculate metric for the original network
Nest <- networklevel(data, index="weighted NODF")

#Calculate metric for the randomized networks
randomized.Nest <- unlist(sapply(nulls, networklevel, index="weighted NODF"))
(Nest - mean(randomized.Nest))/sd(randomized.Nest) # Z value
Nest.sig <- sum(randomized.Nest>Nest)/length(randomized.Nest) # p value
Nest.sig

#Plot the observed value against the distribution of randomized values
plot(density(randomized.Nest), main="Observed vs. randomized",
     xlim=c(min((Nest), min(randomized.Nest)), 
            max((Nest), max(randomized.Nest))))
abline(v=Nest, col="red", lwd=2, xlab="")

Nest #observed
mean(randomized.Nest) #randomized mean
sd(randomized.Nest) #randomized SD
(Nest - mean(randomized.Nest))/sd(randomized.Nest) # Z-value
sum(randomized.Nest>(Nest)) / length(randomized.Nest) #randomized > observed
sum(randomized.Nest<(Nest)) / length(randomized.Nest) #randomized < observed


##### COMPOUND TOPOLOGY 

#Calculate the desired nestedness metric (here WNODA) for the original network.
obs.com <- unlist(bipartite::nest.smdm(x = data, 
                                       constraints = Part, #Input the modular structured recovered from step 2
                                       weighted = T, #By considering the edge weights, you are choosing WNODA
                                       decreasing = "abund"))

#Check the scores
obs.com

#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = data, 
                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Check what those probabilities look like
Pij
dim(Pij)

#Generate randomized networks with the null model of your choice, considering the interaction probabilities calculated before. 
nulls.com <- RestNullModel(M = data, 
                           Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                           Numbernulls = permutations, #This step may take long, so start experimenting with low values
                           Print.null = T, 
                           allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                           return.nonrm.species = F, 
                           connectance = T, byarea = T, 
                           R.partitions = row.Part, 
                           C.partitions = col.Part)

#Calculate the nestedness within and between modules
rest.nest <- nest.smdm(data, constraints = Part, 
                       weighted = T, 
                       decreasing = "abund", 
                       sort = T)

unlist(rest.nest)

null.com <- sapply(nulls.com, 
                   function(x) bipartite::nest.smdm(x = x,
                                                    constraints = Part, 
                                                    weighted = T, 
                                                    decreasing = "abund"))
WNODA.null.com <- unlist(null.com[3,])
WNODAsm.null.com <- unlist(null.com[8,])
WNODAdm.null.com <- unlist(null.com[9,])

#Plot the observed nestedness value against the distribution of randomized values
par(mfrow = c(1,3))

plot(density(WNODA.null.com), xlim=c(min(obs.com[3], min(WNODA.null.com)),
                                     max(obs.com[3], max(WNODA.null.com))), 
     main="observed vs. randomized", xlab = "WNODA matrix")
abline(v=obs.com[3], col="red", lwd=2)

plot(density(WNODAsm.null.com), xlim=c(min(obs.com[8], min(WNODAsm.null.com)),
                                       max(obs.com[8], max(WNODAsm.null.com))), 
     main="observed vs. randomized", xlab = "WNODAsm matrix")
abline(v=obs.com[8], col="red", lwd=2)    

plot(density(WNODAdm.null.com), xlim=c(min(obs.com[9], min(WNODAdm.null.com)),
                                       max(obs.com[9], max(WNODAdm.null.com))), 
     main="observed vs. randomized", xlab = "WNODAdm matrix")
abline(v=obs.com[9], col="red", lwd=2)    

par(mfrow = c(1,1))

#Estimate the p-values

#Nestedness in th entire network
praw.WNODA <- sum(WNODA.null.com>obs.com[3]) / length(WNODA.null.com)
p.WNODA <- ifelse(praw.WNODA > 0.5, 1- praw.WNODA, praw.WNODA)    # P-value
p.WNODA

#Nestedness within the modules
praw.WNODAsm <- sum(WNODAsm.null.com>obs.com[8]) / length(WNODAsm.null.com)
p.WNODAsm <- ifelse(praw.WNODAsm > 0.5, 1- praw.WNODAsm, praw.WNODAsm)    # P-value
p.WNODAsm

#Nestedness between the modules
praw.WNODAdm <- sum(WNODAdm.null.com>obs.com[9]) / length(WNODAdm.null.com)
p.WNODAdm <- ifelse(praw.WNODAdm > 0.5, 1- praw.WNODAdm, praw.WNODAdm)    # P-value
p.WNODAdm


##### PLOT THE COMPOUND TOPOLOGY

par(mfrow = c(1,1))

#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp <- bipartite::sortmatrix(matrix = data, topology = "compound", 
                                   sort_by = "weights", 
                                   row_partitions = row.Part, 
                                   col_partitions = col.Part)

#Assign colors for the modules
modcol <- rainbow((length(unique(Part))), alpha=1, s = 1, v = 1)

#Plot the matrix
png("figures/compound.png", width = 3000, height = 3000, res = 300)
plotmatrix(data.comp$matrix, 
           row_partitions = data.comp$row_partitions, 
           col_partitions = data.comp$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol,
           within_color = modcol, 
           between_color = "lightgrey")
dev.off()


##### EXPORT A SUMMARY OF THE TOPOLOGICAL RESULTS

sink(file = "results/results_topology.txt")

paste("Topological analysis of the nocturnal pollination network")
paste("Queiroz et al. 2020, Biotropica")
cat("\n")
paste("The network has", nrow(data), "rows and", ncol(data), "columns.")
cat("\n")
paste("The network's specialization (H2) is", round(Spec, 2),",", "P =", round(Spec.sig, 2))
cat("\n")
paste("The network's modularity (DIRT_LPA+) is", round(Mod@likelihood, 2), ",", "P =", round(Mod.sig, 2), ",", "and it contains", length(unique(Part)),"modules.")
cat("\n")
paste("The network's nestedness (WNODF) is", round(Nest/100, 2),",", "P =", round(Nest.sig, 2))
cat("\n")
paste("The network shows the following scores of nestedness (WNODA):")
cat("\n")
paste("Entire network =", round(rest.nest$WNODAmatrix/100, 2), ",", "P =", round(p.WNODA, 2))
cat("\n")
paste("Between the modules =", round(rest.nest$WNODA_DM_matrix/100, 2), ",", "P =", round(p.WNODAdm, 2))
cat("\n")
paste("Within the modules =", round(rest.nest$WNODA_SM_matrix/100, 2), ",", "P =", round(p.WNODAsm, 2))
cat("\n")

sink(file = NULL, )



################################################################################
##### SPECIES LEVEL ANALYSIS (CENTRALITY)
################################################################################


# Specialization (d') 
d <- specieslevel(data,index="d")
dplants <- d$`higher level`
write.csv(dplants, "results/dplants.csv") # writing a separate csv file   

# Betweenness centrality (BC)
BC <- specieslevel(data, index="betweenness")
BCplants <- BC$higher
write.csv(BCplants, "results/BCplants.csv") # writing a separate csv file

# Normalized degree (nk)
ND <-ND(data, normalised=T)
NDplants <- ND$higher
write.csv(NDplants, "results/NDplants.csv") # writing a separate csv file


##### Plot centrality metrics by syndrome

#Import data
plants <- read.xls("data/plants.xlsx", h=T) # reading compiled spreadsheet with species & metrics classified by guild

# Change reference level for GLMs
ord <-  ordered(plants$Guild, levels = c("sphin", "chiro", "other"))

# BC

tiff("figures/bc.tiff", width = 12, height = 15, units = "cm", res = 300)
ggplot(plants, aes(x=ord, y=bc, fill=Guild)) + 
  ylab("BC")+ xlab("")+ ylim(0, 0.15) +
  scale_fill_manual(values=c("darkolivegreen1", "sandybrown", "orchid1"))+
  geom_boxplot(width=0.5, color="black") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5) ,
        axis.title.y = element_text(color="black", face ="italic", size =23),
        axis.text= element_text(color="black", size=19),
        legend.position = "none")
dev.off()

# d'

tiff("figures/d.tiff", width = 12, height = 15, units = "cm", res = 300)

ggplot(plants, aes(x=ord, y=d, fill=Guild)) + 
  ylab("d'")+ xlab("")+ ylim(0, 0.8) +
  scale_fill_manual(values=c("darkolivegreen1", "sandybrown", "orchid1"))+
  geom_boxplot(width=0.5, color="black") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5) ,
        axis.title.y = element_text(color="black", face ="italic", size =23),
        axis.text= element_text(color="black", size=19),
        legend.position = "none")

dev.off()

# ND

tiff("figures/nk.tiff", width = 12, height = 15, units = "cm", res = 300)

ggplot(plants, aes(x=ord, y=nk, fill=Guild)) + 
  ylab("nk")+ xlab("")+ ylim(0, 1.1) +
  scale_fill_manual(values=c("darkolivegreen1", "sandybrown", "orchid1"))+
  geom_boxplot(width=0.5, color="black") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5) ,
        axis.title.y = element_text(color="black", face ="italic", size =23),
        axis.text= element_text(color="black", size=19),
        legend.position = "none")

dev.off()


##### Run GLMs to compare centrality scores

table(plants$Guild)
plants$Guild <- relevel(plants$Guild, ref="chiro") #changing reference level for GLMs
plants$Guild <- relevel(plants$Guild, ref="sphin")
plants$Guild <- relevel(plants$Guild, ref="other")

# d'

glmd <- glm(plants$d ~ plants$Guild, family=quasibinomial("logit"))
summary(glmd)
glm_d <- anova(glmd, test = "Chisq")
glm_d

## BC

glmbc <- glm(plants$bc ~ plants$Guild, family=quasibinomial("logit"))
summary(glmbc)
glm_bc <- anova(glmbc, test = "Chisq")
glm_bc

## nk

glmnk <- glm(plants$nk ~ plants$Guild, family=quasibinomial("logit"))
summary(glmnk)
glm_k <- anova(glmnk, test = "Chisq")
glm_k

#Export the results
sink(file = "results/results_centrality.txt")

paste("Comparison of centrality by syndrome")
paste("Queiroz et al. 2020, Biotropica")
cat("\n")
paste(capture.output(glmd))
paste(capture.output(glm_d))
cat("\n")
paste(capture.output(glmbc))
paste(capture.output(glm_bc))
cat("\n")
paste(capture.output(glmnk))
paste(capture.output(glm_k))

sink(file = NULL, )



################################################################################
##### MORPHOMETRIC ANALYSIS
################################################################################


##### Import morphology data
morph_plants<-read.xls("data/morph_pla.xlsx", h=T)
morph_pol <- read.xls("data/morph_pol.xlsx", h=T)

str(morph_plants)
str(morph_pol)

# Change reference level for GLMs
morph_plants$module <- relevel(morph_plants$module, ref="bat")
morph_pol$module <- relevel(morph_pol$module, ref="hawk1")


##### Run GLMs to compare modules

# Pollinator tongues
glm_pol <- glm(morph_pol$length_pol~morph_pol$module, family=gaussian())
summary(glm_pol)
anova(glm_pol, test = "Chisq")

# Floral width (w) and length (l)
glm_pla_l <- glm(morph_plants$length_pla~morph_plants$module, family=gaussian())
glm_pla_w <- glm(morph_plants$width_pla~morph_plants$module, family=gaussian())
summary(glm_pla_l)
anova(glm_pla_l, test = "Chisq")
summary(glm_pla_w)
anova(glm_pla_w, test = "Chisq")

### Plot
ggmorph<-read.xls("data/morph_graph.xlsx",h=T) 

tiff("figures/morph.tiff", width = 24, height = 10, units = "cm", res = 300)
ggplot(ggmorph, aes(x=module, y=measure, fill=variable)) +
  ylab("Measure (mm)")+ xlab("Modules")+ ylim(0, 150) +
  scale_fill_manual(values=c("slateblue", "goldenrod1", "coral1"))+
  geom_boxplot(width=0.5, color="black", position = position_dodge(width=0.5)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5) ,
        axis.title.y = element_text(color="black", size =20),
        axis.title.x = element_text(color="black", size =20),
        axis.text= element_text(color="black", size=19), legend.position = "none")
dev.off()



################################################################################
##### SAMPLING COMPLETENESS ANALYSIS
################################################################################


# Loading interaction data for Chao1 estimator
sampbat<- read.xls("data/sampbat.xlsx", h=T)
estimateR(sampbat, index =c("chao")) 
str(sampbat)

samphawk <- read.xls("data/samphawk.xlsx", h=T)
estimateR(samphawk, index =c("chao"))
str(samphawk)

# Loading interaction data for rarefaction curve drawing
sampling_bats <- read.xls("data/sampling_bats.xlsx", h=T)
curve_bat<- specaccum(sampling_bats, method="rarefaction")

sampling_hawkmoths <- read.xls("data/sampling_hawkmoths.xlsx", h=T)
curve_hawk<- specaccum(sampling_hawkmoths, method="rarefaction")

# Plot curves
tiff("figures/sampling.tiff", width = 20, height = 20, units = "cm", res = 600)
par(mfrow=c(1,2))

plot(curve_hawk, ci.type = "poly", xvar = "individuals", ci.lty=0, ylab = NA,
     xlab = NA,
     ci.col=rgb(0.7, 0, 0.2, 0.3), ylim=c(0,30))
abline(h=22.2, lty=1, col=rgb(0.7, 0, 0.2, 0.3), lwd=2.5)
abline(h=(22.2+0.6195203), lty=3, col=rgb(0.7, 0, 0.2, 0.3), lwd=2)
abline(h=(22.2-0.6195203), lty=3, col=rgb(0.7, 0, 0.2, 0.3), lwd=2)

plot(curve_bat, ci.type = "poly", xvar = "individuals", ci.lty=0,
     ci.col=rgb(0, 0, 0.5, 0.3), ylim=c(0,30), ylab=NA, xlab=NA)
abline(h=14, lty=1, col=rgb(0, 0, 0.5, 0.3), lwd=2.5)
abline(h=(14-2.283481), lty=3, col=rgb(0, 0, 0.5, 0.3), lwd=2)
abline(h=(14+2.283481), lty=3, col=rgb(0, 0, 0.5, 0.3), lwd=2)

dev.off()