require(data.table)
require(ggplot2)
require(DESeq2)

#Import dataframes containing read counts
#Column 2 is for human read counts, column 4 is for lemur read counts
#Column 18,19,20 are for genomic locations
##Import CPD replicates
human_lemur_cpd_rep1 = fread("human_lemur_cpd.tsv",select = c(2,4,18,19,20),header=F)
human_lemur_cpd_rep2 = fread("human_lemur_cpd_rep2.tsv",select = c(2,4,18,19,20),header=F)
##Import 6-4 replicates
human_lemur_64_rep1 = fread("human_lemur_64.tsv",select = c(2,4,18,19,20),header=F)
human_lemur_64_rep2 = fread("human_lemur_64_rep2.tsv",select = c(2,4,18,19,20),header=F)
##Import GM12878 cells, column 2 is rep1, column 4 is rep2 read counts
gm = fread("gm12878.tsv",select = c(2,4,18,19,20),header=F)

#Paste columns 18,19,20 together to identify homologies and merge them into unified dataframe
# and define a function to trim unwanted columns and create a homology dataframe on position
trimmer=function(df,col1,col2) {
    df$pos = paste0(df$V18,df$V19,df$V20)
    df= df[,c(1,2,6)]
    colnames(df) = c(col1,col2,"pos")
    return(df)
}

gm = trimmer(gm,"gm1","gm2")
cpd1 = trimmer(human_lemur_cpd_rep1,"human_cpd1","lemur_cpd1")
cpd2 = trimmer(human_lemur_cpd_rep2,"human_cpd2","lemur_cpd2")
six4_1 = trimmer(human_lemur_64_rep1,"human_six4_1","lemur_six4_1") 
six4_2 = trimmer(human_lemur_64_rep2,"human_six4_2","lemur_six4_2") 

#merge all datasets based on position
cpds = merge(cpd1,cpd2,by="pos",all=T)
six4 = merge(six4_1,six4_2,by="pos",all=T)
merged = merge(cpds,six4,by="pos",all=T)
merged = merge(merged,gm,by="pos",all=T)

#delete column "pos" for numerical operations since its class is character 
merged = merged[,-c("pos")]

#coerce NAs to Zeros
merged[is.na(merged)] <- 0

#Keep rows where there are more than 10 reads
keep <- rowSums(merged) >= 10
merged <- merged[keep,]

#Perform Variance stabilizing transformation from DESeq2
vsd <- vst(as.matrix(merged), blind=FALSE)

#Perform PCA
pca= prcomp(vsd, center=T, scale=T)
plot(pca$rotation[,1],pca$rotation[,2], xlab = "PC1", ylab = "PC2")
text(pca$rotation[,1],pca$rotation[,2], row.names(pca$rotation), cex=0.5, pos=4)

#
