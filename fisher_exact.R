# Fisher test v7 # update 29.11.13 this is the script I am using for analysis v1.1
# This includes multiple testing correction. 
# To use put all genes in mapman bins. Column name = "Mapman". It is best to modify the mapman bins to reduce complexity. 
# Based on http://cgrlucb.wikispaces.com/Functional+Enrichment+Analysis.
# Where k = c1, m-k = c2, n-k = c3, N-m-n-k = c4

# load the packages

library(scales)
library(ggplot2)
library(reshape)

# Import the data

mesophyll <- read.csv ("setariamesophyll.csv", header = TRUE, sep = "\t")
bundlesheath <- read.csv ("setariabundlesheath.csv", header = TRUE, sep = "\t")
both <- rbind(mesophyll, bundlesheath)

# set the variables 

# -log10 pvalue (0.59/0.6 (for 0.25), 0.82 (for 0.15), 1 (for 0.1), 1.3 (for 0.05), 2 (for 0.01))
# mutliple testing correction, used 'fdr', used 0.59 for setaria and 1 for maize (0.25, and 0.1)
# filter out low count groups, 

MTC <- "fdr"      
pval <- 1
title <- ""
filter <- 5 # set as 0
zz <- 9 # max log2p value

# may need to alter column names so mapman is with a capital M

# first get a list of unique mapman catagories

unique <- as.character(unique(both$Mapman))  

## loop for mesophyll

result <- data.frame(cat=character(1), count=numeric(1), p=numeric(1), stringsAsFactors=FALSE)

for(category in unique) {
  categoryname <- mesophyll[ which(mesophyll$Mapman==category), ] 
  c1 <- length(categoryname$id) # mesophyll DE genes in cat
  count <- c1
  categorynamex <- both[ which(both$Mapman==category), ] 
  c5 <- length(categorynamex$id) # genes in cat all
  c3 <- c5 - c1 # genes in cat not DE
  totalDE <- length(mesophyll$id)
  c2 <- totalDE - c1 # genes not in cat DE
  notDE <- length(both$id) - length(mesophyll$id) # total not DE
  c4 <- notDE - c3
  table <- matrix(c(c1,c3,c2,c4), nrow=2, ncol=2)
  fisher.result <- fisher.test(table, alternative='g')
  p.value <- fisher.result$p.value
  meslogp <- -log10(p.value)
  result <- rbind(result, c(category, count, p.value)) 
}

# save the raw counts to a backup file for retrieval

mesophyllcounts <- result
mesophyllcounts <- mesophyllcounts[with(mesophyllcounts, order(cat)), ]

## loop for bundle sheath

result2 <- data.frame(cat=character(1), count=numeric(1), p=numeric(1), stringsAsFactors=FALSE)

for(category in unique) {
  categoryname <- bundlesheath[ which(bundlesheath$Mapman==category), ] 
  c1 <- length(categoryname$id) # 
  count <- c1
  categorynamex <- both[ which(both$Mapman==category), ] 
  c5 <- length(categorynamex$id) # genes in cat all
  c3 <- c5 - c1 # genes in cat not DE
  totalDE <- length(bundlesheath$id)
  c2 <- totalDE - c1 # genes not in cat DE
  notDE <- length(both$id) - length(bundlesheath$id) # total not DE
  c4 <- notDE - c3
  table <- matrix(c(c1,c3,c2,c4), nrow=2, ncol=2)
  fisher.result <- fisher.test(table, alternative='g')
  p.value <- fisher.result$p.value
  meslogp <- -log10(p.value)
  result2 <- rbind(result2, c(category, count, p.value)) 
}

# save the counts to a file

bundlesheathcounts <- result2
bundlesheathcounts <- bundlesheathcounts[with(bundlesheathcounts, order(cat)), ]

# merge raw counts and p values

unsigmergedcounts <- merge(mesophyllcounts, bundlesheathcounts, by = "cat", all.x = TRUE, all.y = TRUE)
# write.csv(unsigmergedcounts, file = "unsigmergedcounts.csv", row.names = FALSE)

# filter the counts with a certain cut off

write.csv(result, file = "temp.csv", row.names = FALSE)
result <- read.csv(file = "temp.csv", header = TRUE)
write.csv(result2, file = "temp.csv", row.names = FALSE)
result2 <- read.csv(file = "temp.csv", header = TRUE)

merged1 <- merge(result, result2, by = "cat")
colnames(merged1)[2] <- "countmes"
colnames(merged1)[3] <- "pmes"
colnames(merged1)[4] <- "countbs"
colnames(merged1)[5] <- "pbs" 
merged1 <- merged1[-1,]
merged1$sum <- (merged1$countmes + merged1$countbs)
merged2 <- subset(merged1, sum >= filter)                 # filter the merged data set

# multiple testing correction and take -log10 for merged set

# for mesophyll

padjmes <- p.adjust(merged2$pmes, method=MTC) # set multiple testing correction method
merged2 <- cbind(merged2, padjmes)
meslogp <- -log10(padjmes)
merged2 <- cbind(merged2, meslogp)

# for bundlesheath

padjbs <- p.adjust(merged2$pbs, method=MTC) # set multiple testing correction method
merged2 <- cbind(merged2, padjbs)
bslogp <- -log10(padjbs)
merged2 <- cbind(merged2, bslogp)

# reorder the whole data frame

merged3 <- merged2[c("cat", "countmes", "pmes", "padjmes", "meslogp", "countbs", "pbs", "padjbs", "bslogp")]
merged7 <- merged2[c("cat", "padjmes", "padjbs")] # for entire set
colnames(unsigmergedcounts)[2] <- "countmes"
colnames(unsigmergedcounts)[3] <- "pmes"
colnames(unsigmergedcounts)[4] <- "countbs"
colnames(unsigmergedcounts)[5] <- "pbs"
merged9 <- merge(unsigmergedcounts, merged7, by = "cat", all.x = TRUE)
merged9 <- merged9[-1,]
write.csv(merged9, file = "alldata.csv", row.names = FALSE)

# subset based on a p value threshold

merged3$meslogp[ merged3$meslogp >= zz ] <- zz
merged3$bslogp[ merged3$bslogp >= zz ] <- zz

mesophyllupregulated <- subset(merged3, meslogp >= pval)
bundlesheathupregulated <- subset(merged3, bslogp >= pval)

# merge the sets

significant <- subset(merged3, meslogp >= pval | bslogp >= pval)  # set threshold (1 for 0.1, 1.3 for 0.05, 2 for 0.01)

# find and replace

significant$cat <- gsub("RNA.regulation of transcription.", "transcription factor - ", significant$cat)
significant$cat <- gsub("transcription factor family", "", significant$cat)

#write a copy of just mes and bs logps to merge with maize later

tomergeforboth <- significant[c("cat", "meslogp", "bslogp")]
write.csv(tomergeforboth, file = "tomergeforboth.csv", row.names = FALSE)

#write the significant results

significant[is.na(significant)] <- 0

write.csv(significant, file = "significant.csv", row.names = FALSE)
significant <- read.csv ("significant.csv", header = TRUE, row.names = 1)

# pull out just the cat and logp values

significant <- significant[c("meslogp", "bslogp")]

## transform the data

rows <- rownames(significant)
significant$catagory <- rows

significant$catagory <- gsub("\\."," ", significant$catagory)

rownames(significant) <- NULL
significant <- significant[c(3,1,2)]
colnames(significant)[2] <- "mesophyll"
colnames(significant)[3] <- "bundle sheath"
significant.m <- melt(significant) # just from here three lines if you are running both together

significant.m <- ddply(significant.m, .(variable), transform, log2p = rescale(value)) # rescale function around value removed

colnames(significant.m)[2] <- "cell"

## draw the plot

q <- ggplot(significant.m, aes(cell, catagory)) + geom_tile(aes(fill = log2p), colour = "white") + 
  scale_fill_gradient(low = "white", high = "red") + 
  labs(title = title) + theme(axis.text.x = element_text(size = 12, face = "bold", colour = "black"),
      axis.text.y = element_text(colour = "black"), plot.title = element_text(size = 14, face = "bold"))
