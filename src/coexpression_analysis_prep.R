library(WGCNA)

# Load data
load(file='rdata/normalized_data.RData')
load(file='rdata/differential_expression.RData')
d.test <- t(d.norm)

# Check for missing values in data
check.missing <- goodSamplesGenes(d.test, verbose = 3)
check.missing$allOK

# Plot data to check for sample outliers
sampleDend <- hclust(dist(d.test), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleDend, main="Sample clustering", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)

# Check if some clinical variables correlate with the clustering
dim(meta.data)
names(meta.data)
# remove columns that hold information we do not need
allTraits <- meta.data[, c(2,4,5,8,9,10,11,12,13)]
allTraits.num <- allTraits
allTraits.num$group <- as.numeric(factor(allTraits.num$group), levels=c('NDC', 'FTD_TDP', 'FTD_TAU'))
allTraits.num$type <- as.numeric(factor(allTraits.num$type), levels=c('0', 'sporadic', 'genetic'))
allTraits.num$mutation <- as.numeric(factor(allTraits.num$mutation,
                                            levels=c('no', 'to do', 'GRN (NBB)', 'GRN (new mut)',
                                                     'MAPT (NBB)', 'C9ORF (NBB)')))
allTraits.num$gender <- as.numeric(factor(allTraits.num$gender, levels=c('m', 'f')))
allTraits.num$TAU..AT8. <- as.numeric(factor(allTraits.num$TAU..AT8., levels=c('0', 'pos but weak', 'pos')))
allTraits.num$TDP43 <- as.numeric(factor(allTraits.num$TDP43, levels=c('0', '?', 'pos but weak', 'pos (threads, few)', 'pos')))
allTraits.num$p62 <- as.numeric(factor(allTraits.num$p62, levels=c('0', 'repeat because background', 'ND', 'pos but weak', 'pos (threads, few)', 'pos')))
allTraits.num$Amyloid.beta.plaques <- as.numeric(factor(allTraits.num$Amyloid.beta.plaques, levels=c('', '0', 'ND', 'pos (few)', 'pos')))
allTraits.num$age <- as.numeric(allTraits.num$age)
dim(allTraits.num)
names(allTraits.num)
collectGarbage()

# Re-cluster samples
sampleDend2 <- hclust(dist(d.test), method="average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(allTraits.num, signed=FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleDend2, traitColors, groupLabels=names(allTraits.num),
                    main="Sample relations to traits")
