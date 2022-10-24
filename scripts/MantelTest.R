########## Mantel test for phylogenetic signal of sars-cov-2 mutation spectrum evolution ######
# Annabel Beichman

########## load packages and set seed ##############
require(vegan) # for mantel test 
require(dplyr)
require(tidyverse)
require(ggplot2)
require(scales)
require(broom) # for tidying distance matrix
set.seed(42)

npermutations=9999999 # 10M permutations for Mantel test
todaysdate=format(Sys.Date(),"%Y%m%d")



### input files ###


# 1. pairwise matrix of phylogenetic distances 

resultsdir="/Users/annabelbeichman/Documents/UW/sars_cov_2_spectrum/SARS2-mut-spectrum/results/"

outdir=paste0(resultsdir,"mantel_test/")
dir.create(outdir,showWarnings = F)

# read in, skipping row 1 which is empty and setting row names equal to column 1 
phylodist <- read.table(paste0(resultsdir,"/clade_founder_tree/clade_founders.mldist"),header=F,skip=1,row.names=1)
# make column names the same as row names:
colnames(phylodist) <- rownames(phylodist)
dim(phylodist) # 19x19 matrix with 0s along the diagonal, good

head(phylodist)

# note there are some clades in this tree that aren't in the spectrum distances , is that correct?

# 2. pairwise (Euclidean) spectrum distances
spectrumdist_E <- read.table(paste0(resultsdir,"/synonymous_mut_rates/clade_rate_distances.csv"),sep=",",header=T)

# need to turn this into square matrix with filled upper and lower triangles
# so need to duplicate the values, switching what is in clade 1 and clade 2

# need to add self-self comparisons for diagonal
clades <- unique(c(spectrumdist_E$clade_1,spectrumdist_E$clade_2)) # gets all the clade names
# want the self-self distances (0) for the diagonal
self_distances_diagonal <- data.frame(clade_1=clades,clade_2=clades,mut_rate_distance=0)
# add this in
spectrumdist_E_plusDiagonal <- bind_rows(spectrumdist_E,self_distances_diagonal)

# pivot wider
spectrumdist_E_matrix <- spectrumdist_E_plusDiagonal %>%
  pivot_wider(names_from = clade_2, values_from = mut_rate_distance) %>%
  column_to_rownames('clade_1') %>%
  as.matrix()

# reorder columns to match order of rownames
spectrumdist_E_matrix <- spectrumdist_E_matrix[,rownames(spectrumdist_E_matrix)]


# this works to fill in the matrix (note both are lower.tri )
spectrumdist_E_matrix[lower.tri(spectrumdist_E_matrix)] <- t(spectrumdist_E_matrix)[lower.tri(spectrumdist_E_matrix)] # this tranposes matrix so upper tri becomes lower tri, then fills in the original matrix with the new lower tri 

# spot check:
spectrumdist_E_matrix
# diagonal is 0 -- check
# is symmetrical -- check
# spot checks match  original file from Jesse -- (21C x 22C = .753990; 22A x 21J = 1.02600) -- check.

# 3. Spectrum dist with Aitchison distances <- to be added

###### order the two matrices in the same way; also removes clades from the phylo matrix that aren't in spectrum distance matrix #######

phylodist_reordered_subset <- phylodist[rownames(spectrumdist_E_matrix),colnames(spectrumdist_E_matrix)] # reorder phylo matrix to match spectrum matrix row/column names
# make sure it looks good with original file
# phylo dist: phylodist["21C","22C"] == phylodist_reordered_subset["21C","22C"] TRUE. good.

# check dimensions and look visually at matrices to make sure it all looks good
dim(phylodist_reordered_subset)
dim(spectrumdist_E_matrix)
head(phylodist_reordered_subset)
head(spectrumdist_E_matrix)
# make sure row/col names are the same 
if(sum(rownames(phylodist_reordered_subset)!=rownames(spectrumdist_E_matrix))!=0){
  stop("something is wrong with your matrix setup")
}

if(sum(colnames(phylodist_reordered_subset)!=colnames(spectrumdist_E_matrix))!=0){
  stop("something is wrong with your matrix setup")
}


########## run the mantel test #########
mantelTestResults_E <- vegan::mantel(phylodist_reordered_subset,spectrumdist_E_matrix,permutations = npermutations) # default method = pearson. 

mantelOutFile_E <- paste0(outdir,todaysdate,".mantelResults.NOTsqrtdist.Euclidean.txt")
sink(file=mantelOutFile_E)
print(mantelTestResults_E)
sink()

# make a little dataframe with the results
mantelTestResult_E_df <- data.frame(method=as.character(mantelTestResults_E$method),statistic=as.numeric(mantelTestResults_E$statistic),permCount=as.numeric(mantelTestResults_E$permutations),signif=as.numeric(mantelTestResults_E$signif))


######### TO DO: mantel test with sqrt distance ############
mantelTestResults_E_sqrtDistance <- vegan::mantel(sqrt(phylodist_reordered_subset),spectrumdist_E_matrix,permutations = npermutations) # default method = pearson. 

mantelOutFileSqrtDistance_E <- paste0(outdir,todaysdate,".mantelResults.sqrtdistance.Euclidean.txt")
sink(file=mantelOutFileSqrtDistance_E)
print(mantelTestResults_E_sqrtDistance)
sink()

# make a little dataframe with the results

mantelTestResultSqrtDistance_E_df <- data.frame(method=as.character(mantelTestResults_E_sqrtDistance$method),statistic=as.numeric(mantelTestResults_E_sqrtDistance$statistic),permCount=as.numeric(mantelTestResults_E_sqrtDistance$permutations),signif=as.numeric(mantelTestResults_E_sqrtDistance$signif))
#### merge distances for plotting #####
phylodist_reordered_subset_table <- tidy(as.dist(phylodist_reordered_subset)) # turn into a table
# spot check
head(phylodist_reordered_subset_table)
head(phylodist_reordered_subset)
# this got rid of self-self distances and duplicates -- good. ends up with dim 91 which is 14 choose 2. check. 
colnames(phylodist_reordered_subset_table) <- c("clade_1","clade_2","phylogenetic_distance")


## combine:
merged_distances_E <- merge(spectrumdist_E,phylodist_reordered_subset_table,by=c("clade_1","clade_2"))
# add labels:

merged_distances_E$comparison_label <- paste0(merged_distances_E$clade_1,"x",merged_distances_E$clade_2)
####### plot results #####


mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_E <- ggplot(merged_distances_E,aes(x=sqrt(phylogenetic_distance),y=mut_rate_distance))+
  geom_point(size=1,alpha=0.85)+
  geom_text(data=mantelTestResultSqrtDistance_E_df,aes(x=0.1,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\nEuclidean spectrum distances; ",unique(mantelTestResultSqrtDistance_E_df$permCount)," permutations; SQRT of phylo distance"))+
  theme_bw()+
  xlab("square root of phylogenetic distance (shared branch length)")+
  ylab("Euclidean distance between mutation spectra")
# making strips not gray

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_E
ggsave(paste0(outdir,todaysdate,"mantelTest.CorrelationPlot.SQRTDISTANCE.Euclidean.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_E,height=5,width=7)

# add labels:
# color based on whether comparison involves omicron?

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_E_labels <- ggplot(merged_distances_E,aes(x=sqrt(phylogenetic_distance),y=mut_rate_distance))+
  #geom_point(size=1,alpha=0.85)+
  geom_text(aes(label=comparison_label),size=2)+
  geom_text(data=mantelTestResultSqrtDistance_E_df,aes(x=0.1,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\nEuclidean spectrum distances; ",unique(mantelTestResultSqrtDistance_E_df$permCount)," permutations; SQRT of phylo distance"))+
  theme_bw()+
  xlab("square root of phylogenetic distance (shared branch length)")+
  ylab("Euclidean distance between mutation spectra")

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_E_labels
ggsave(paste0(outdir,todaysdate,"mantelTest.CorrelationPlot.SQRTDISTANCE.Euclidean.labels.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_E_labels,height=5,width=7)
