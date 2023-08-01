similarityModules1 <- as.matrix(read.table(paste0('KRAS G12/output/S1/moduleSimilarity.tsv'), sep='\t'))
similarityModules2 <- as.matrix(read.table(paste0('KRAS G12/output/S2/moduleSimilarity.tsv'), sep='\t'))
similarityModules3 <- as.matrix(read.table(paste0('KRAS G12/output/S3/moduleSimilarity.tsv'), sep='\t'))

averageSimilarityModules1 <- median(apply(similarityModules1,2,max))
averageSimilarityModules2 <- median(apply(similarityModules2,2,max))
averageSimilarityModules3 <- median(apply(similarityModules3,2,max))

similaritySubnetworks1 <- as.matrix(read.table(paste0('KRAS G12/output/S1/subnetworkSimilarity.tsv'), sep='\t'))
similaritySubnetworks2 <- as.matrix(read.table(paste0('KRAS G12/output/S2/subnetworkSimilarity.tsv'), sep='\t'))
similaritySubnetworks3 <- as.matrix(read.table(paste0('KRAS G12/output/S3/subnetworkSimilarity.tsv'), sep='\t'))

averageSimilaritySubnetworks1 <- median(apply(similaritySubnetworks1,2,max))
averageSimilaritySubnetworks2 <- median(apply(similaritySubnetworks2,2,max))
averageSimilaritySubnetworks3 <- median(apply(similaritySubnetworks3,2,max))

similaritySSubnetworks1 <- as.matrix(read.table(paste0('KRAS G12/output/S1/sSubnetworkSimilarity.tsv'), sep='\t'))
colorsA <- as.vector(t(read.table('KRAS G12/output/S1/training/subnetworks.tsv', sep='\t', header = FALSE)))
colorsB <- as.vector(t(read.table('KRAS G12/output/S1/validation/subnetworks.tsv', sep='\t', header = FALSE)))
exclusiveA <- colorsA[!(colorsA == 'grey')]
exclusiveB <- colorsB[!(colorsB == 'grey')]
similaritySSubnetworks1 <- similaritySSubnetworks1[colorsB, colorsA]
similaritySSubnetworks2 <- as.matrix(read.table(paste0('KRAS G12/output/S2/sSubnetworkSimilarity.tsv'), sep='\t'))
colorsA <- as.vector(t(read.table('KRAS G12/output/S2/training/subnetworks.tsv', sep='\t', header = FALSE)))
colorsB <- as.vector(t(read.table('KRAS G12/output/S2/validation/subnetworks.tsv', sep='\t', header = FALSE)))
exclusiveA <- colorsA[!(colorsA == 'grey')]
exclusiveB <- colorsB[!(colorsB == 'grey')]
similaritySSubnetworks2 <- similaritySSubnetworks2[colorsB, colorsA]
similaritySSubnetworks3 <- as.matrix(read.table(paste0('KRAS G12/output/S3/sSubnetworkSimilarity.tsv'), sep='\t'))
colorsA <- as.vector(t(read.table('KRAS G12/output/S3/training/subnetworks.tsv', sep='\t', header = FALSE)))
colorsB <- as.vector(t(read.table('KRAS G12/output/S3/validation/subnetworks.tsv', sep='\t', header = FALSE)))
exclusiveA <- colorsA[!(colorsA == 'grey')]
exclusiveB <- colorsB[!(colorsB == 'grey')]
similaritySSubnetworks3 <- similaritySSubnetworks3[colorsB, colorsA]

averageSimilaritySSubnetworks1 <- median(apply(similaritySSubnetworks1,2,max))
averageSimilaritySSubnetworks2 <- median(apply(similaritySSubnetworks2,2,max))
averageSimilaritySSubnetworks3 <- median(apply(similaritySSubnetworks3,2,max))

resultSimilarity1 <- as.matrix(read.table(paste0('KRAS G12/output/S1/similarity.tsv'), sep='\t'))
resultSimilarity2 <- as.matrix(read.table(paste0('KRAS G12/output/S2/similarity.tsv'), sep='\t'))
resultSimilarity3 <- as.matrix(read.table(paste0('KRAS G12/output/S3/similarity.tsv'), sep='\t'))

GSEASimilarity1 <- as.matrix(read.table(paste0('KRAS G12/output/S1/gsea.tsv'), sep='\t'))
GSEASimilarity2 <- as.matrix(read.table(paste0('KRAS G12/output/S2/gsea.tsv'), sep='\t'))
GSEASimilarity3 <- as.matrix(read.table(paste0('KRAS G12/output/S3/gsea.tsv'), sep='\t'))

end <- data.frame('Modules' = c(averageSimilarityModules1, averageSimilarityModules2, averageSimilarityModules3),
                  'Subnetworks' = c(averageSimilaritySubnetworks1, averageSimilaritySubnetworks2, averageSimilaritySubnetworks3),
                  'Significant subnetworks' = c(averageSimilaritySSubnetworks1, averageSimilaritySSubnetworks2, averageSimilaritySSubnetworks3),
                  'GO genes' = c(resultSimilarity1['genes',1], resultSimilarity2['genes',1], resultSimilarity3['genes',1]),
                  'BP' = c(resultSimilarity1['BP',1], resultSimilarity2
                           ['BP',1], resultSimilarity3['BP',1]),
                  'MF' = c(resultSimilarity1['MF',1], resultSimilarity2['MF',1], resultSimilarity3['MF',1]),
                  'GSEA BP' = c(GSEASimilarity1['BP',1], GSEASimilarity2['BP',1], GSEASimilarity3['BP',1]),
                  'GSEA MF' = c(GSEASimilarity1['MF',1], GSEASimilarity2['MF',1], GSEASimilarity3['MF',1]))

rownames(end) <- c('Run 1', 'Run 2', 'Run 3')

png(file=paste0('KRAS G12/output/similarity.png'), width=1920, height=1020)
boxplot(end, main = 'Average similarity')
dev.off()

R1 <- data.frame(Type=character(0), Similarity=numeric(0))
for (l in apply(similarityModules1,2,max)){
  R1[nrow(R1) + 1,] <-  c('Module', as.numeric(l))
}
for (l in apply(similaritySubnetworks1,2,max)){
  R1[nrow(R1) + 1,] <- c('Subnetwork', as.numeric(l))
}
for (l in apply(similaritySSubnetworks1,2,max)){
  R1[nrow(R1) + 1,] <- c('Significant subnetwork', as.numeric(l))
}

R1[nrow(R1) + 1,] <- c('GO genes', as.numeric(resultSimilarity1['genes',1]))
R1[nrow(R1) + 1,] <- c('GO BP', as.numeric(resultSimilarity1['BP',1]))
R1[nrow(R1) + 1,] <- c('GO MF', as.numeric(resultSimilarity1['MF',1]))
R1[nrow(R1) + 1,] <- c('GSEA BP', as.numeric(GSEASimilarity1['BP',1]))
R1[nrow(R1) + 1,] <- c('GSEA MF', as.numeric(GSEASimilarity1['MF',1]))

R1$Similarity <- as.numeric(R1$Similarity)

png(file=paste0('KRAS G12/output/S1/similarity.png'), width=1920, height=1020)
ggplot(R1, aes(x=Type, y=Similarity, fill=Type)) + 
  geom_violin() +
  geom_boxplot(width=0.05, color="white") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.75, color = 'black', binwidth = 1/100) +
  scale_y_continuous(name = "Similarity", limits=c(0,1), breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1)) +
  scale_x_discrete(limits=c("Module", "Subnetwork", "Significant subnetwork", 'GO genes', 'GO BP', 'GO MF', 'GSEA BP', 'GSEA MF'))
dev.off()

R2 <- data.frame(Type=character(0), Similarity=numeric(0))
for (l in apply(similarityModules2,2,max)){
  R2[nrow(R2) + 1,] <-  c('Module', as.numeric(l))
}
for (l in apply(similaritySubnetworks2,2,max)){
  R2[nrow(R2) + 1,] <- c('Subnetwork', as.numeric(l))
}
for (l in apply(similaritySSubnetworks1,2,max)){
  R2[nrow(R2) + 1,] <- c('Significant subnetwork', as.numeric(l))
}

R2[nrow(R2) + 1,] <- c('GO genes', as.numeric(resultSimilarity1['genes',1]))
R2[nrow(R2) + 1,] <- c('GO BP', as.numeric(resultSimilarity1['BP',1]))
R2[nrow(R2) + 1,] <- c('GO MF', as.numeric(resultSimilarity1['MF',1]))
R2[nrow(R2) + 1,] <- c('GSEA BP', as.numeric(GSEASimilarity1['BP',1]))
R2[nrow(R2) + 1,] <- c('GSEA MF', as.numeric(GSEASimilarity1['MF',1]))

R2$Similarity <- as.numeric(R2$Similarity)

png(file=paste0('KRAS G12/output/S2/similarity.png'), width=1920, height=1020)
ggplot(R2, aes(x=Type, y=Similarity, fill=Type)) + 
  geom_violin() + 
  geom_boxplot(width=0.05, color="white") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.75, color = 'black', binwidth = 1/100) +
  scale_y_continuous(name = "Similarity", limits=c(0,1), breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1)) +
  scale_x_discrete(limits=c("Module", "Subnetwork", "Significant subnetwork", 'GO genes', 'GO BP', 'GO MF', 'GSEA BP', 'GSEA MF'))
dev.off()

R3 <- data.frame(Type=character(0), Similarity=numeric(0))
for (l in apply(similarityModules3,2,max)){
  R3[nrow(R3) + 1,] <-  c('Module', as.numeric(l))
}
for (l in apply(similaritySubnetworks3,2,max)){
  R3[nrow(R3) + 1,] <- c('Subnetwork', as.numeric(l))
}
for (l in apply(similaritySSubnetworks3,2,max)){
  R3[nrow(R3) + 1,] <- c('Significant subnetwork', as.numeric(l))
}

R3[nrow(R3) + 1,] <- c('GO genes', as.numeric(resultSimilarity3['genes',1]))
R3[nrow(R3) + 1,] <- c('GO BP', as.numeric(resultSimilarity3['BP',1]))
R3[nrow(R3) + 1,] <- c('GO MF', as.numeric(resultSimilarity3['MF',1]))
R3[nrow(R3) + 1,] <- c('GSEA BP', as.numeric(GSEASimilarity3['BP',1]))
R3[nrow(R3) + 1,] <- c('GSEA MF', as.numeric(GSEASimilarity3['MF',1]))

R3$Similarity <- as.numeric(R3$Similarity)

png(file=paste0('KRAS G12/output/S3/similarity.png'), width=1920, height=1020)
ggplot(R3, aes(x=Type, y=Similarity, fill=Type)) + 
  geom_violin() + 
  geom_boxplot(width=0.05, color="white") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.75, color = 'black', binwidth = 1/100) +
  scale_y_continuous(name = "Similarity", limits=c(0,1), breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1)) +
  scale_x_discrete(limits=c("Module", "Subnetwork", "Significant subnetwork", 'GO genes', 'GO BP', 'GO MF', 'GSEA BP', 'GSEA MF'))
dev.off()

modules1 <- data.frame(matrix(nrow = length(colnames(similarityModules1)), ncol = 5))
rownames(modules1) <- colnames(similarityModules1)
colnames(modules1) <- c('Name1', 'Name2', 'Similarity', 'Size1', 'Size2')
modules1[, 'Name1'] <- rownames(modules1)
modules1[, 'Name2'] <- rownames(similarityModules1)[apply(similarityModules1,2,which.max)]
modules1[, 'Similarity'] <- apply(similarityModules1,2,max)
for (i in 1:nrow(modules1)){
  modules1[i, 'Size1'] <- length(as.vector(t(read.table(paste0('KRAS G12/output/S1/training/coexpression/', row.names(modules1[i,]), '.tsv'), sep='\t', header = FALSE))))
  modules1[i, 'Size2'] <- length(as.vector(t(read.table(paste0('KRAS G12/output/S1/validation/coexpression/', modules1[i, 'Name2'], '.tsv'), sep='\t', header = FALSE))))
}
png(file=paste0('KRAS G12/output/S1/modules.png'), width=1920, height=1020)
ggplot(modules1, aes(x=Size1,y=Similarity)) + 
  geom_point(aes(size=Size2, colour = Name1))
dev.off()

modules2 <- data.frame(matrix(nrow = length(colnames(similarityModules2)), ncol = 5))
rownames(modules2) <- colnames(similarityModules2)
colnames(modules2) <- c('Name1', 'Name2', 'Similarity', 'Size1', 'Size2')
modules2[, 'Name1'] <- rownames(modules2)
modules2[, 'Name2'] <- rownames(similarityModules2)[apply(similarityModules2,2,which.max)]
modules2[, 'Similarity'] <- apply(similarityModules2,2,max)
for (i in 1:nrow(modules2)){
  modules2[i, 'Size1'] <- length(as.vector(t(read.table(paste0('KRAS G12/output/S2/training/coexpression/', row.names(modules2[i,]), '.tsv'), sep='\t', header = FALSE))))
  modules2[i, 'Size2'] <- length(as.vector(t(read.table(paste0('KRAS G12/output/S2/validation/coexpression/', modules2[i, 'Name2'], '.tsv'), sep='\t', header = FALSE))))
}
png(file=paste0('KRAS G12/output/S2/modules.png'), width=1920, height=1020)
ggplot(modules2, aes(x=Size1,y=Similarity)) + 
  geom_point(aes(size=Size2, colour = Name1))
dev.off()

modules3 <- data.frame(matrix(nrow = length(colnames(similarityModules3)), ncol = 5))
rownames(modules3) <- colnames(similarityModules3)
colnames(modules3) <- c('Name1', 'Name2', 'Similarity', 'Size1', 'Size2')
modules3[, 'Name1'] <- rownames(modules3)
modules3[, 'Name2'] <- rownames(similarityModules3)[apply(similarityModules3,2,which.max)]
modules3[, 'Similarity'] <- apply(similarityModules3,2,max)
for (i in 1:nrow(modules3)){
  modules3[i, 'Size1'] <- length(as.vector(t(read.table(paste0('KRAS G12/output/S3/training/coexpression/', row.names(modules3[i,]), '.tsv'), sep='\t', header = FALSE))))
  modules3[i, 'Size2'] <- length(as.vector(t(read.table(paste0('KRAS G12/output/S3/validation/coexpression/', modules3[i, 'Name2'], '.tsv'), sep='\t', header = FALSE))))
}
png(file=paste0('KRAS G12/output/S3/modules.png'), width=1920, height=1020)
ggplot(modules3, aes(x=Size1,y=Similarity)) + 
  geom_point(aes(size=Size2, colour = Name1))
dev.off()

end <- data.frame('Modules' = c(averageSimilarityModules1, averageSimilarityModules2, averageSimilarityModules3),
                  'Subnetworks' = c(averageSimilaritySubnetworks1, averageSimilaritySubnetworks2, averageSimilaritySubnetworks3),
                  'Significant subnetworks' = c(averageSimilaritySSubnetworks1, averageSimilaritySSubnetworks2, averageSimilaritySSubnetworks3),
                  'GO genes' = c(resultSimilarity1['genes',1], resultSimilarity2['genes',1], resultSimilarity3['genes',1]),
                  'BP' = c(resultSimilarity1['BP',1], resultSimilarity2['BP',1], resultSimilarity3['BP',1]),
                  'MF' = c(resultSimilarity1['MF',1], resultSimilarity2['MF',1], resultSimilarity3['MF',1]),
                  'GSEA BP' = c(GSEASimilarity1['BP',1], GSEASimilarity2['BP',1], GSEASimilarity3['BP',1]),
                  'GSEA MF' = c(GSEASimilarity1['MF',1], GSEASimilarity2['MF',1], GSEASimilarity3['MF',1]))

rownames(end) <- c('Run 1', 'Run 2', 'Run 3')

png(file=paste0('KRAS G12/output/similarity.png'), width=1920, height=1020)
boxplot(end, main = 'Average similarity')
dev.off()

print(paste('Similarity for modules across runs:', mean(c(averageSimilarityModules1, averageSimilarityModules2, averageSimilarityModules3))))
print(paste('Similarity for subnetworks across runs:', mean(c(averageSimilaritySubnetworks1, averageSimilaritySubnetworks2, averageSimilaritySubnetworks3))))
print(paste('Similarity for significant subnetworks across runs:', mean(c(averageSimilaritySSubnetworks1, averageSimilaritySSubnetworks2, averageSimilaritySSubnetworks3))))
print(paste('Similarity for significant genes in subnetworks across runs:', mean(resultSimilarity1['genes',1], resultSimilarity2['genes',1], resultSimilarity3['genes',1])))
print(paste('Similarity for significant BP terms in subnetworks across runs:', mean(resultSimilarity1['BP',1], resultSimilarity2['BP',1], resultSimilarity3['BP',1])))
print(paste('Similarity for significant MF terms in subnetworks across runs:', mean(resultSimilarity1['MF',1], resultSimilarity2['MF',1], resultSimilarity3['MF',1])))
print(paste('Similarity for GSEA BP across runs:', mean(GSEASimilarity1['BP',1], GSEASimilarity2['BP',1], GSEASimilarity3['BP',1])))
print(paste('Similarity for GSEA MF across runs:', mean(GSEASimilarity1['MF',1], GSEASimilarity2['MF',1], GSEASimilarity3['MF',1])))