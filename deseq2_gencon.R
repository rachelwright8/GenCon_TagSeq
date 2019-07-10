library(DESeq2) # for differential gene expression analysis
library(tidyverse) # for data wrangling and plotting

# set to working directory on your computer
setwd("~/Documents/aims2015/tagseq_19apr19/gencon_mega_git/")

# Look at the sequencing stats --------
allStats <- read.delim("sequencing_stats.txt")
head(allStats)

# Calculate mean reads after trimming and mapping
# & percentages of reads remaining after trimming and mapping
allStats %>% mutate(trimPerc = round(postTrim/preTrim,2), 
                    mappedPerc = round(mapped/postTrim,2)) %>% 
  summarise(meanTrimReads = mean(postTrim), meanMappedReads = mean(mapped), 
            meanTrimPerc = round(mean(trimPerc),2), meanMappedPerc = round(mean(mappedPerc),2))

# meanTrimReads meanMappedReads meanTrimPerc meanMappedPerc
# 1687253         1276020         0.91           0.76

# Load expression table and meta data --------
load("counts_and_traits.RData")

# Look at the expression matrix. Rows are genes. Columns are samples.
head(countData)

# Look at the sample meta data. "sam" matches the sample name in the expression matrix.
# "tank" = factor representing the tank that coral fragment was in
# "genet" = factor representing the genetic identity of the fragment
# "reef" = factor representing where the colony was collected
# "treat" = factor representing experimental treatments
head(metaData)

table(metaData$treat)
# a = All stressors (see below)
# b = bacteria (Vibrio owensii)
# c = control (27°C, pH 8.0, no Vibrio)
# h = heat (30°C)
# p = low pH (7.8)

# a  b  c  h  p 
# 58 59 60 60 57 

# Look at the trait data
head(traitData)
# there is trait data for 939 samples, but "only" expression data for 294 ;)

# For example, plot the average change in coral color by treatment
traitData %>% group_by(treat) %>% summarise(meanColorChange = mean(color.log, na.rm=T),
                                            sdColorChange = sd(color.log, na.rm=T)) %>%
  ggplot(aes(x=treat, y=meanColorChange, fill=treat)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=meanColorChange-sdColorChange, 
                    ymax=meanColorChange+sdColorChange),
                width=.2,
                position=position_dodge(.9))+
  scale_fill_manual(values=c("grey", "orange", "yellow", "red", "blue"))+
  theme_bw()

# have fun!