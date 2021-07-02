setwd("~/Dropbox/aims2015/tagseq_28june21/deseq2_host/")

library(DESeq2) # for differential gene expression analysis
library(tidyverse) # for wrangling and plotting
library(VennDiagram) # for venn diagram
library(pheatmap) # for pretty heatmaps
library(arrayQualityMetrics) # for sample outlier analysis

# sequencing stats --------
countPreTrim <- read.delim("../countsPreTrim.txt", header=F)
head(countPreTrim)
tail(countPreTrim)
# remove the first two rows of description
countPreTrim <- countPreTrim[-c(1,2),]
countPreTrim$V1 <- gsub(".fq", "", countPreTrim$V1)
names(countPreTrim)[2] <- "preTrim"
head(countPreTrim)

countPostTrim <- read.delim("../countsTrimmed.txt", header=F)
head(countPostTrim)
countPostTrim <- countPostTrim[-c(1,2),]
countPostTrim$V1 <- gsub(".fq.trim", "", countPostTrim$V1)
names(countPostTrim)[2] <- "postTrim"
head(countPostTrim)

mappedCounts <- read.delim("../AlignmentRate.txt", header=F)
head(mappedCounts)
mapped <- as.numeric(substr(mappedCounts$V1,1,5))

allStats <- merge(countPreTrim,countPostTrim,by="V1")
allStats <- cbind(allStats,mapped)
head(allStats)

allStats %>% summarise(meanPreTrim = mean(preTrim),
                       meanPostTrim = mean(postTrim),
                       meanTrimPercRemaining = mean(postTrim/preTrim),
                       meanMappedPerc = mean(mapped))

# meanPreTrim meanPostTrim meanTrimPercRemaining meanMappedPerc
# 1,859,609     95,3188.8             0.5154209       84.24156

# Load expression table --------
countdata0 <- read.delim("../allcounts.coral.txt")
dim(countdata0)

# look at the corners
countdata0[1:4,1:3] # looks good
countdata0[(nrow(countdata0)-6):nrow(countdata0),1:3] # has non-gene rows
countdata0 <- countdata0 %>% filter(grepl("Amillepora",countdata0$X)) # get rid of non-gene rows
countdata0[(nrow(countdata0)-6):nrow(countdata0),1:3] # Looks good now
countdata0[(nrow(countdata0)-4):nrow(countdata0),(length(countdata0)-4):length(countdata0)] # has "NA" column
countdata0$X.1 <- NULL
countdata0[(nrow(countdata0)-4):nrow(countdata0),(length(countdata0)-4):length(countdata0)] # looks good now

# look at the corners again
countdata0[1:4,1:3]
countdata0[(nrow(countdata0)-4):nrow(countdata0),1:3]
countdata0[(nrow(countdata0)-4):nrow(countdata0),(length(countdata0)-4):length(countdata0)]
countdata0[1:4,(length(countdata0)-4):length(countdata0)]

# make genes the row names
row.names(countdata0) <- countdata0$X
countdata0$X <- NULL

# make column names nicer
head(names(countdata0))

tank <- as.factor(substr(names(countdata0),2,3))
summary(tank)
# 01 02 03 04 05 06 09 11 18 20 21 22 23 25 26 27 28 29 30 32 35 36 37 38 
# 25  3  9  7  8  3  4 21 16 11  7 16 17  2 24 21  8 26 10 19  6 22  2  7 

genet <- as.factor(substr(names(countdata0),5,6))
summary(genet)
# 01 04 05 06 07 09 11 12 14 15 16 18 19 20 22 23 24 25 29 30 32 33 35 38 39 40 50 53 54 55 
# 10 10 10 10 10 10  9 10 10 10  9 10 10 10 10 10 10 10 10  8 10 10 10 10 10 10 10  8 10 10

names(countdata0) <- paste("t",tank,"_g",genet,sep="")
head(names(countdata0))

# fix the clones
# 11 and 12 = 60 (Rib)
# 24 and 32 = 61 (Pandora)
# 18 and 20 = 62 (Rib)
# add "a" and "b" to indicate the replicates
names(countdata0) <- gsub("g11", "g60a", names(countdata0))
names(countdata0) <- gsub("g12", "g60b", names(countdata0))
names(countdata0) <- gsub("g24", "g61a", names(countdata0))
names(countdata0) <- gsub("g32", "g61b", names(countdata0))
names(countdata0) <- gsub("g18", "g62a", names(countdata0))
names(countdata0) <- gsub("g20", "g62b", names(countdata0))
head(names(countdata0))

# now pad the other names with "a" to indicate the single replicate for that genet/treat
names(countdata0) <- str_pad(names(countdata0), 8, 
                            side = "right", pad = "a")
head(names(countdata0))

# count data ready to go!
countdata <- countdata0
head(countdata)

##################
# get trait data
##################

load("../../AllTraits_May2018/traits/allTraitData_Survival.RData")
head(d)

# summarize by tank and genet (average biological replicates)
d_cond_data <- d %>% ungroup() %>% 
  dplyr::select(tank,geno,tank,reef,treat) %>% 
  rename(genet = geno) %>%
  mutate(sam = as.factor(paste("t",tank,"_g",genet,"a",sep=""))) %>%
  distinct(sam, .keep_all=T)
head(d_cond_data)

# make rows with the other replicate
dup_rows <- d_cond_data %>% filter (genet %in% c("60","61","62"))
dup_rows$sam <- gsub("a","b",dup_rows$sam)
head(dup_rows)

# add to the other data frame
conds <- rbind(d_cond_data,dup_rows)
head(conds)

conds <- conds[conds$sam %in% names(countdata),]
head(conds)

# sort to match the order in the counts file
conds <- conds[match(names(countdata), conds$sam),]
table(conds$sam == names(countdata))

# summarize samples
conds %>% group_by(treat) %>% summarize(n = n())

# Construct data object ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = conds,
  design = ~ treat + genet)

save(countdata, conds, dds, file="dds.Rdata")
load("dds.Rdata")

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("treat", "genet"), force=T)
save(vsd, e, file="VSD_arrayQM.Rdata")
# no samples failed more than one test. Keep all samples.

# load all data --------------
load("dds.Rdata")
load("VSD_arrayQM.Rdata")

# Set base mean minimum ---------------------------------------------------
means <- apply(countdata,1,mean)
table(means>3)
# FALSE  TRUE 
# 17987 10201 

means3 <- names(means[means>3])
head(means3)
length(means3)
# 10201

countFilt <- countdata[row.names(countdata) %in% means3,]
head(countFilt)

totalCountsFilt <- colSums(countFilt)
totalCountsFilt
min(totalCountsFilt) # 91368
max(totalCountsFilt) # 698852
mean(totalCountsFilt) # 314409.1

# check sample order
table(names(countFilt) == as.vector(conds$sam))

# Reconstruct data object (filtered) ------------------------------

ddsFilt <- DESeqDataSetFromMatrix(
  countData = countFilt,
  colData = conds,
  design = ~ treat + genet)

save(countdata, conds, dds, countFilt, ddsFilt, file="ddsFilt.Rdata")

# DESeq -------------------------------------------------------------------

load("ddsFilt.Rdata")

# new VSD ------
vsd.filt <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFilt)
# this may take awhile

#---Results
resBac <- results(deds, independentFiltering = F, contrast = c("treat","b","c"))
resBac
summary(resBac, alpha = 0.05)
# LFC > 0 (up)       : 2317, 23%
# LFC < 0 (down)     : 2369, 23%
# outliers [1]       : 1, 0.0098%

resComb <- results(deds, independentFiltering = F, contrast = c("treat","a","c"))
resComb
summary(resComb, alpha = 0.05)
# LFC > 0 (up)       : 1750, 17%
# LFC < 0 (down)     : 1963, 19%
# outliers [1]       : 1, 0.0098%

resHeat <- results(deds, independentFiltering = F, contrast = c("treat","h","c"))
resHeat
summary(resHeat, alpha = 0.05)
# LFC > 0 (up)       : 456, 4.5%
# LFC < 0 (down)     : 417, 4.1%
# outliers [1]       : 1, 0.0098%

resPco2 <- results(deds, independentFiltering = F, contrast = c("treat","p","c"))
resPco2
summary(resPco2, alpha = 0.05)
# LFC > 0 (up)       : 141, 1.4%
# LFC < 0 (down)     : 76, 0.75%
# outliers [1]       : 1, 0.0098%

resCombvsBac <- results(deds, independentFiltering = F, contrast = c("treat","a","b"))
resCombvsBac
# log2 fold change (MLE): treat a vs b
summary(resCombvsBac, alpha = 0.05)
# LFC > 0 (up)       : 1133, 11%
# LFC < 0 (down)     : 1127, 11%
# outliers [1]       : 1, 0.0098%

resCombvsHeat <- results(deds, independentFiltering = F, contrast = c("treat","a","h"))
resCombvsHeat 
summary(resCombvsHeat, alpha = 0.05)
# LFC > 0 (up)       : 1692, 17%
# LFC < 0 (down)     : 2013, 20%

resCombvsPCO2 <- results(deds, independentFiltering = F, contrast = c("treat","a","p"))
resCombvsPCO2
summary(resCombvsPCO2, alpha = 0.05)
# LFC > 0 (up)       : 1941, 19%
# LFC < 0 (down)     : 2215, 22%

# plot summary of DEGs
alpha = 0.05
bacUP <- resBac %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, padj) %>% 
  filter(padj<alpha & log2FoldChange>0) %>%
  pull(geneID)
bacDOWN <- resBac %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, padj) %>% 
  filter(padj<alpha & log2FoldChange<0) %>%
  pull(geneID)
combUP <- resComb %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, padj) %>% 
  filter(padj<alpha & log2FoldChange>0) %>%
  pull(geneID)
combDOWN <- resComb %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, padj) %>% 
  filter(padj<alpha & log2FoldChange<0)  %>%
  pull(geneID)
heatUP <- resHeat %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, padj) %>% 
  filter(padj<alpha & log2FoldChange>0) %>%
  pull(geneID)
heatDOWN <- resHeat %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, padj) %>% 
  filter(padj<alpha & log2FoldChange<0)  %>%
  pull(geneID)
pco2UP <- resPco2 %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, padj) %>% 
  filter(padj<alpha & log2FoldChange>0) %>%
  pull(geneID)
pco2DOWN <- resPco2 %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, padj) %>% 
  filter(padj<alpha & log2FoldChange<0)  %>%
  pull(geneID)

candidatesUP <- list("Bacteria"=bacUP,"Combined"=combUP,"Heat"=heatUP, "pCO2"=pco2UP)

pdf(file="fig_vennDEGSUP.pdf", height=5,width=5)
prettyvennUP <- venn.diagram(
  x = candidatesUP,
  filename = NULL,
  col = "transparent",
  fill = c("yellow", "darkorange", "red", "blue"),
  alpha = 0.5,
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08, 0.08, 0.08),
  cat.pos = 1
);
grid.draw(prettyvennUP)
dev.off()

candidatesDOWN <- list("Bacteria"=bacDOWN,"Combined"=combDOWN,"Heat"=heatDOWN, "pCO2"=pco2DOWN)

pdf(file="fig_vennDEGSDOWN.pdf", height=5,width=5)
prettyvennDOWN <- venn.diagram(
  x = candidatesDOWN,
  filename = NULL,
  col = "transparent",
  fill = c("yellow", "darkorange", "red", "blue"),
  alpha = 0.5,
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08, 0.08, 0.08),
  cat.pos = 1
);
grid.draw(prettyvennDOWN)
dev.off()

# what are the overlapping genes?
Reduce(intersect, list(bacDOWN, combDOWN, heatDOWN, pco2DOWN))
# [1] "Amillepora02995" "Amillepora04844" "Amillepora07197" "Amillepora20489" "Amillepora20490" "Amillepora35339"

Reduce(intersect, list(bacUP, combUP, heatUP, pco2UP))
# [1] "Amillepora12443" "Amillepora16461" "Amillepora20865"

# Save/Load Data ----------------------------------------------------------

save(vsd.filt, countFilt, conds, ddsFilt, deds, 
     resBac, resComb, resHeat, resPco2, resCombvsBac, resCombvsHeat, resCombvsPCO2,
     file="results.Rdata")
load("results.Rdata")


# Explore with plots ------------------------------------------------------
#Sample distance heatmap
pheatmap(cor(assay(vsd.filt)),border_color=NA, main="SampleHeatmap")

# Diagnostics -------------------------------------------------------------

#Dispersions plot
plotDispEsts(deds, main="Dispersion Plot Response")

#MA plot
plotMA(resBac, ylim = c(-1, 1), main="MA Plot Bacteria") 
plotMA(resComb, ylim = c(-1, 1), main="MA Plot Combined") 
plotMA(resHeat, ylim = c(-1, 1), main="MA Plot Heat") 
plotMA(resPco2, ylim = c(-1, 1), main="MA Plot Pco2") 

# Write results for making heatmaps ---------------------------------------

###--------------Get pvals & LFCs
head(resBac)
valsBac <- cbind(resBac$log2FoldChange, resBac$pvalue, resBac$padj)
head(valsBac)
colnames(valsBac)=c("LFC.b", "pval.b", "padj.b")
length(valsBac[,1])
table(complete.cases(valsBac))

head(resComb)
valsComb <- cbind(resComb$log2FoldChange, resComb$pvalue, resComb$padj)
head(valsComb)
colnames(valsComb)=c("LFC.a", "pval.a", "padj.a")
length(valsComb[,1])
table(complete.cases(valsComb))

head(resHeat)
valsHeat <- cbind(resHeat$log2FoldChange, resHeat$pvalue, resHeat$padj)
head(valsHeat)
colnames(valsHeat)=c("LFC.h","pval.h", "padj.h")
length(valsHeat[,1])
table(complete.cases(valsHeat))

head(resPco2)
valsPco2 <- cbind(resPco2$log2FoldChange, resPco2$pvalue, resPco2$padj)
head(valsPco2)
colnames(valsPco2)=c("LFC.p","pval.p", "padj.p")
length(valsPco2[,1])
table(complete.cases(valsPco2))

#Make VSD and pvals table
vsdpvals <- cbind(assay(vsd.filt), valsBac, valsComb, valsHeat, valsPco2)
vsdpvals[1:2,(ncol(vsdpvals)-12):ncol(vsdpvals)]
write.csv(vsdpvals, "29june21_WaldtoCont_VSDandPVALS.csv", quote=F)

# Write results for GO/KOG analysis -------------------------------------------

# by -log p-value
resBac %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, pvalue) %>%
  mutate(logp = round(-log(pvalue+1e-10,10),1),
         sign = ifelse(log2FoldChange>0,1,-1)) %>%
  mutate(logpsign = logp*sign) %>%
  #head # run through this line to check that everything looks OK, then comment out and move ahead
  select(geneID, logpsign) %>%
  write_csv("resBac_logpval.csv") # CHANGE THE FILE NAME!

resComb %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, pvalue) %>%
  mutate(logp = round(-log(pvalue+1e-10,10),1),
         sign = ifelse(log2FoldChange>0,1,-1)) %>%
  mutate(logpsign = logp*sign) %>%
  # head # run through this line to check that everything looks OK, then comment out and move ahead
  select(geneID, logpsign) %>%
  write_csv("resComb_logpval.csv") # CHANGE THE FILE NAME!

resHeat %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, pvalue) %>%
  mutate(logp = round(-log(pvalue+1e-10,10),1),
         sign = ifelse(log2FoldChange>0,1,-1)) %>%
  mutate(logpsign = logp*sign) %>%
  # head # run through this line to check that everything looks OK, then comment out and move ahead
  select(geneID, logpsign) %>%
  write_csv("resHeat_logpval.csv") # CHANGE THE FILE NAME!

resPco2 %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>%
  select(geneID, log2FoldChange, pvalue) %>%
  mutate(logp = round(-log(pvalue+1e-10,10),1),
         sign = ifelse(log2FoldChange>0,1,-1)) %>%
  mutate(logpsign = logp*sign) %>%
  # head # run through this line to check that everything looks OK, then comment out and move ahead
  select(geneID, logpsign) %>%
  write_csv("resPco2_logpval.csv") # CHANGE THE FILE NAME!

# Principal Coordinates Analysis ---------------
library(vegan) # for adonis
library(ape) # for pcoa

# variance stabilized expression data
exp <- data.frame(assay(vsd.filt))
head(exp)

# condition data (check to make sure they match)
head(conds)
table(conds$sam == names(exp))

# compute dissimilarity indices
dd.veg <- vegdist(t(exp), "manhattan")
div.dd.veg <- dd.veg/1000
head(div.dd.veg)

# perform PERMANOVA  
set.seed(1)
adonisRes <- adonis(t(exp)~treat+reef,data=conds,method="manhattan")
adonisRes
 
#           Df  SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# treat       4  190169497 47542374  5.7340 0.07106  0.001 ***
# reef        3  114855662 38285221  4.6175 0.04292  0.001 ***

# compute principal coordinate decomposition
dd.pcoa <- pcoa(div.dd.veg)
head(dd.pcoa)
scores <- dd.pcoa$vectors

# plotting PCoA----
margin <- 1

# play around with these numbers to see different axes
xaxis <- 1
yaxis <- 2

# PCoA for mid by site type
pdf(file = "fig_pca_host.pdf", width = 7, height = 6)
# quartz()
pca_plot <- plot(scores[,xaxis], scores[,yaxis],type="n", 
                 main = "Coral Gene Expression",
                 xlim=c(min(scores[,xaxis])-margin,max(scores[,xaxis])+margin),
                 ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
                 mgp=c(2.3,1,0),
                 xlab=paste("PCo", xaxis," (", 
                            round(dd.pcoa$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
                 ylab=paste("PCo", yaxis," (", 
                            round(dd.pcoa$values$Relative_eig[yaxis]*100,1),"%)",sep=""))
points(scores[conds$treat=="c",xaxis],
       scores[conds$treat=="c",yaxis],
       col="black", pch=15, cex=1.5) +
points(scores[conds$treat=="a",xaxis],
         scores[conds$treat=="a",yaxis],
         col="orange", pch=15, cex=1.5) +
points(scores[conds$treat=="b",xaxis],
         scores[conds$treat=="b",yaxis],
         col="yellow", pch=15, cex=1.5) +
points(scores[conds$treat=="h",xaxis],
         scores[conds$treat=="h",yaxis],
         col="red", pch=15, cex=1.5) +
points(scores[conds$treat=="p",xaxis],
         scores[conds$treat=="p",yaxis],
         col="blue", pch=15, cex=1.5)

ordiellipse(scores,conds$treat,draw="polygon",
            col = c("black","orange","yellow","red","blue"),
            label = F)

# # legend of sites 
legend("bottomleft", 
       c("Control", "Vibrio", "Combined"), 
       text.col = c("black","gold","orange"), 
       cex = 1,
       ncol = 1,
       bty = "n")
legend("bottomright", 
       c("Heat", "Low pH"), 
       text.col = c("red","blue"), 
       cex = 1,
       ncol = 1,
       bty = "n")

#insert p value
# legend("topleft",
       # paste("Treatment p = ",adonisRes$aov.tab$`Pr(>F)`[1], sep=" "),
       # cex=1, bty='n')
# legend("topright",
       # paste("Reef p = ",adonisRes$aov.tab$`Pr(>F)`[2], sep=" "),
       # cex=1, bty='n')
dev.off()

# make heatmap of differentially expressed genes ----
# Load annotations
gg <- read.delim("~/Dropbox/genomes/amil_genome/amil_gene2geneName.tab",
                sep="\t",header=T, skip=3)
head(gg)
names(gg) <- c("V1", "V2")

# Format expression data
vsdpvals <- read.csv("29june21_WaldtoCont_VSDandPVALS.csv")
head(vsdpvals)
names(as.data.frame(vsdpvals))
exp <- vsdpvals[-c(grep("pval",names(vsdpvals)),
                   grep("padj",names(vsdpvals)),
                        grep("LFC",names(vsdpvals)))]
head(exp)

#     Make p-value cut-offs
sig.heat <- row.names(vsdpvals[vsdpvals$padj.h<0.0001 & !is.na(vsdpvals$padj.h),])
length(sig.heat)

# sig expression
exp.sig <- exp[row.names(exp) %in% sig.heat,]
nrow(exp.sig)

# species heatmap ---------
# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(exp.sig)){
  s=subset(gg,V1==i)
  if (length(s[,1])>0){
    gnames=append(gnames,paste(s$V2[1],i,sep="."))
    expg=rbind(expg,exp.sig[i,])
  } 
}
row.names(expg)=gnames
expl=expg
means=apply(expl,1,mean) # means of rows
explc=expl-means # subtracting them
# explc_sort <- explc[sorted.order]
# names(explc_sort)

# shorten the really long names manually
row.names(explc_sort)[3] <- "Component of U1 snRNA.Mcavernosa01384"

# make colors 
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),
                               bias=0.8)(100)

# cluster plot
quartz()
pheatmap(as.matrix(explc_sort),
         color = heat.colors,cex=0.9,border_color=NA,
         clustering_distance_rows="correlation", cluster_cols=F)


