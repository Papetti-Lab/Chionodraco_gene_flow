### --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ###
### This file contains the command lines used in R to run the analysed presented in the paper "Limited interspecific gene flow in the evolutionary history of the icefish Chionodraco spp." ###
### --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ###


# Set the working directory
setwd("path/to/your/working/direcotory")

	### PCoA ###

# Load required packages
library(vcfR)
library(adegenet)
library(ggplot2)

# Load data
vcf <- read.vcfR("dataset2.vcf")
# Convert vcfR object to genlight object
x <- vcfR2genlight(vcf)

# Define species groups
x$pop <- as.factor(rep(c("hamatus", "myersi", "rastrospinosus", "wilsoni"), c(14, 13, 38, 20)))

## Define colours for each group
species.colours <- c("red1", "gold", "dodgerblue4", "springgreen")

# Compute PCoA
tab.x <- tab(x, freq=TRUE, NA.method="asis")
pco.x <- dudi.pco(dist(tab.x), scannf=TRUE, ful=FALSE)


# Calculate contributions of the displayed principal components to the variance
var.axis.1 <- pco.x$eig[1]/sum(pco.x$eig)*100
var.axis.1
var.axis.2 <- pco.x$eig[2]/sum(pco.x$eig)*100
var.axis.2

# Plot
ggplot(pco.x$li, aes(x=A1, y=A2)) +
  geom_point(aes(fill=x$pop), size=5, shape=21, alpha=0.8) +
  theme_classic() +
  theme_bw(base_size=18) +
  xlab(paste("PC1 ", "(", round(var.axis.1, 2), " %", ")", sep="")) + 
  ylab(paste("PC2 ", "(", round(var.axis.2, 2), " %", ")", sep="")) +
  scale_fill_manual(name="Species:",labels=c(bquote(italic("C. hamatus")),
                               bquote(italic("C. myersi")),
                               bquote(italic("C. rastrospinosus")),
                               bquote(italic("C. wilsoni"))), values=species.colours) +
  theme(legend.position=c(0.8, 0.3), legend.background=element_rect(linetype="solid",linewidth=0.5,colour ="black"))





  ### PCoA only for C. rastrospinosus populations ###

# Load data
rastro.vcf <- read.vcfR("dataset2-rastro.vcf")
# Convert vcfR object to genlight object
y <- vcfR2genlight(rastro.vcf)

# Define populations groups
y$pop <- as.factor(rep(c("SO", "AP", "AS", "WS"), c(11, 11, 7, 9)))

## Define colours for each group
rastro.colours <- c("#FF69B4", "#FFA500", "#333333", "#008B8B")

# Compute PCoA
tab.y <- tab(y, freq=TRUE, NA.method="asis")
pco.y <- dudi.pco(dist(tab.y), scannf=TRUE, ful=FALSE)


# Calculate contributions of the displayed principal components to the variance
y.var.axis.1 <- pco.y$eig[1]/sum(pco.y$eig)*100
y.var.axis.1
y.var.axis.2 <- pco.y$eig[2]/sum(pco.y$eig)*100
y.var.axis.2

# Plot
ggplot(pco.y$li, aes(x=A1, y=A2)) +
  geom_point(aes(fill=y$pop), size=5, shape=21, alpha=0.8) +
  theme_classic() +
  theme_bw(base_size=18) +
  xlab(paste("PC1 ", "(", round(y.var.axis.1, 2), " %", ")", sep="")) + 
  ylab(paste("PC2 ", "(", round(y.var.axis.2, 2), " %", ")", sep="")) +
  scale_fill_manual(name="Populations:",
                    labels=c("Antarctic Peninsula","Antarctic Sound","South Orkneys","Weddel Sea"),
                    values=rastro.colours) +
  theme(legend.position="bottom")





	### Fst ###

# Load required package
library(StAMPP)

# Load data
vcf <- read.vcfR("dataset2.vcf")
x <- vcfR2genlight(vcf)

# Define groups
x$pop <- as.factor(rep(c("hamatus", "myersi", "rastrospinosus", "wilsoni"), c(14, 13, 38, 20)))

# Compute Fst with associated p-values
fst.stampp <- stamppFst(x, nboots = 1000, percent = 95, nclusters = 8)

# Look at the results
fst.stampp$Fsts
fst.stampp$Pvalues





  ### Plot Admixture results ###

# Required package
library(pophelper)

# Import cross-validation errors
cv <- read.table("dataset2-cv_errors.txt", quote="\"")
# Plot cross-validation errors
ggplot(data=cv, aes(x=V1, y=V2)) +
  geom_line(linetype = "dashed") +
  geom_point(size=4) +
  theme_bw(base_size=20) +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(name="K", breaks=cv$V1) +
  ylab("Cross-validation error")

# Import Admixture output
qvalues <- readQ("dataset2.4.Q")
# Define species groups
cat <- data.frame(Species=rep(c("C. hamatus", "C. myersi", "C. rastrospinosus", "C. wilsoni"),
                                         c(14, 13, 38, 20)))
# Plot
plotQ(qvalues, showindlab=F, showlegend=F, showsp=T, showyaxis=T, showticks=T, 
      clustercol=c("springgreen", "red1", "gold", "dodgerblue4"),
      splab=c("K=4"), grplab=cat, grplabface="italic", grplabsize=1.2,
      exportpath=getwd(), imgtype="pdf")





  ### Plot TreeMix results ###

# Load required packages
library(RColorBrewer)
library(R.utils)
# Import plotting functions provided with TreeMix
source("path/to/treemix-1.13/src/plotting_funcs.R")
# I have modified the original functions for plotting to make them more customisable and display clearer plots.
# Place the script "custom_plotting_funcs.R" in the same folders where "plotting_funcs.R" is and import it.
source("path/to/treemix-1.13/src/custom_plotting_funcs.R")

# Define the prefix of the TreeMix outputs (it is the same name used for the input file)
prefix="dataset2"


# Plot with one migration edge allowed
pdf("plot.pdf",width=21, height=15)
par(mar=c(9,0,8,18))
custom_plot_tree("dataset2.1", font=3, cex=5, lwd=5, arrow=0.6, lwd.arrow=15, cex.colorbar=4, cex.mse=4,
             cex.x.label=4, cex.x.values=4, line.x=7, padj.x=0.7, disp = 0.002)
title("1 edge", cex.main=5, line=3)
dev.off()


# Plot trees for all the number of edges that were allowed
for(edge in 0:5){
  pdf(paste0("path/to/output/folder/",prefix,".",edge,".pdf"), width=21, height=15)
  par(mar=c(9,0,8,18))
  custom_plot_tree(paste0(prefix,".",edge), font=3, cex=5, lwd=5, arrow=0.6, lwd.arrow=15, cex.colorbar=4, cex.mse=4,
               cex.x.label=4, cex.x.values=4, line.x=7, padj.x=0.7, disp = 0.002)
  title(paste(edge,"edges"), cex.main=5, line=3)
  dev.off()
}


# Plot residuals for all the number of edges that were allowed
for(edge in 0:5){
  pdf(paste0("path/to/output/folder/",prefix,".",edge,".residuals.pdf"),
      width=8.98, height=5.74)
  par(mar=c(9,13,2,0))
  custom_plot_resid(paste0(prefix,".",edge), pop_order="species-order.txt", cex=2, cex.bar=1.6, font=3, y.pos=-0.07)
  title(paste(edge,"edges"), cex.main=2)
  dev.off()
}



### Plot TreeMix results when considering separated populations ###
# Unlike the previous case, the plots differ depending on the number of migration edges allowed. They must therefore be adapted on a case-by-case basis.
prefix="dataset2_pops"

# 1 edge
pdf("path/to/output/folder/dataset2_pops.1.pdf", width=25, height=20)
par(mar=c(9,0,7,0))
custom_plot_tree("dataset2_pops.1", font=1, cex=5, lwd=5, arrow=0.6, lwd.arrow=8, cex.colorbar=4, cex.mse=4,
                 cex.x.label=4, cex.x.values=4, line.x=7, padj.x=0.7, disp = 0.001, se.y=-0.05, cb.y=-0.1)
title(paste("1 edge"), cex.main=5, line=3)
dev.off()


# 2 edges
pdf("path/to/output/folder/dataset2_pops.2.pdf", width=25, height=20)
par(mar=c(9,0,7,0))
custom_plot_tree("dataset2_pops.2", font=1, cex=5, lwd=5, arrow=0.6, lwd.arrow=8, cex.colorbar=4, cex.mse=4,
                 cex.x.label=4, cex.x.values=4, line.x=7, padj.x=0.7, disp = 0.001, se.y=-0.05)
title(paste("2 edges"), cex.main=5, line=3)
dev.off()


# 3 edges
pdf("path/to/output/folder/dataset2_pops.3.pdf", width=25, height=20)
par(mar=c(9,0,7,0))
custom_plot_tree("dataset2_pops.3", font=1, cex=5, lwd=5, arrow=0.6, lwd.arrow=8, cex.colorbar=4, cex.mse=4,
                 cex.x.label=4, cex.x.values=4, line.x=7, padj.x=0.7, disp = 0.001, se.y=-0.05, cb.y=-0.1)
title(paste("3 edges"), cex.main=5, line=3)
dev.off()


# 4 edges
pdf("path/to/output/folder/dataset2_pops.4.pdf", width=25, height=20)
par(mar=c(9,0,7,0))
custom_plot_tree("dataset2_pops.4", font=1, cex=5, lwd=5, arrow=0.6, lwd.arrow=8, cex.colorbar=4, cex.mse=4,
                 cex.x.label=4, cex.x.values=4, line.x=7, padj.x=0.7, disp = 0.001, se.y=-0.05)
title(paste("4 edges"), cex.main=5, line=3)
dev.off()


# 5-10 edges
for(edge in 5:10){
pdf(paste0("path/to/output/folder/","dataset2_pops.",edge,".pdf"), width=25, height=20)
par(mar=c(9,0,7,0))
custom_plot_tree(paste0(prefix,".",edge), font=1, cex=5, lwd=5, arrow=0.6, lwd.arrow=8, cex.colorbar=4, cex.mse=4,
                 cex.x.label=4, cex.x.values=4, line.x=7, padj.x=0.7, disp = 0.001, se.y=-0.05, cb.y=0.25)
title(paste(edge,"edges"), cex.main=5, line=3)
dev.off()
}



# Plot residuals 
# One migration edge allowed
pdf("path/to/output/folder/dataset2_pops.1.residuals.pdf", width=8.98, height=5.74)
par(mar=c(7,9,2,0))
custom_plot_resid("dataset2_pops.1", pop_order="pop-order.txt", cex=2, cex.bar=1.6, y.pos=-0.05)
title(paste("1 edge"), cex.main=2)
dev.off()

# Two to ten migration edges allowed
for(edge in 2:10){
  pdf(paste0("path/to/output/folder/","dataset2_pops",".",edge,".residuals.pdf"), width=8.98, height=5.74)
  par(mar=c(7,9,2,0))
  custom_plot_resid(paste0(prefix,".",edge), pop_order="pop-order.txt", cex=2, cex.bar=1.6, y.pos=-0.05)
  title(paste(edge,"edges"), cex.main=2)
  dev.off()
}





	### Plot site frequency spectra ###

# Import needed functions (downloaded from http://cmpg.unibe.ch/software/fastsimcoal27/additionalScripts.html)

source("/path/to/ParFileInterpreter_VS.r")
source("/path/to/utilFscOutput.r")

# Define the tag used to identify the model for which you want to plot the site frequency spectra
scenario <- "scenario_name"

# Comparison 1_0
obssfs <- read.table(paste(scenario, "_jointMAFpop1_0.obs", sep=""), skip=1, stringsAsFactors = F, header=T)
expsfs <- read.table(paste(scenario, "_jointMAFpop1_0.txt", sep=""), stringsAsFactors = F, header=T)
plot2dSFS(obsSFS=obssfs, expSFS=expsfs, xtag="hamatus", ytag="rastrospinosus", minentry=1)
plot_relDiff2dSFS(obsSFS=obssfs, expSFS=expsfs, xtag="hamatus", ytag="rastrospinosus", minentry=1)
# Comparison 2_0
obssfs <- read.table(paste(scenario, "_jointMAFpop2_0.obs", sep=""), skip=1, stringsAsFactors = F, header=T)
expsfs <- read.table(paste(scenario, "_jointMAFpop2_0.txt", sep=""), stringsAsFactors = F, header=T)
plot2dSFS(obsSFS=obssfs, expSFS=expsfs, xtag="myersi", ytag="rastrospinosus", minentry=1)
plot_relDiff2dSFS(obsSFS=obssfs, expSFS=expsfs, xtag="myersi", ytag="rastrospinosus", minentry=1)
# Comparison 2_1
obssfs <- read.table(paste(scenario, "_jointMAFpop2_1.obs", sep=""), skip=1, stringsAsFactors = F, header=T)
expsfs <- read.table(paste(scenario, "_jointMAFpop2_1.txt", sep=""), stringsAsFactors = F, header=T)
plot2dSFS(obsSFS=obssfs, expSFS=expsfs, xtag="myersi", ytag="hamatus", minentry=1)
plot_relDiff2dSFS(obsSFS=obssfs, expSFS=expsfs, xtag="myersi", ytag="hamatus", minentry=1)






