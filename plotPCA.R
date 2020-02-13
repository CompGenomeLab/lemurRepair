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
#Summary PCA
summary(pca)
#Label rows
pca_data = data.frame(pca$rotation[,1],pca$rotation[,2],row.names = rownames(pca$rotation))
colnames(pca_data)[1] = "pc1"
colnames(pca_data)[2] = "pc2"
pca_data$CellLine = c("Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","B Lymphocytes","B Lymphocytes") 
pca_data$Organism = c("Human","Mouse Lemur","Human","Mouse Lemur","Human","Mouse Lemur","Human","Mouse Lemur","Human","Human")
pca_data$DamageType = c("CPD","CPD","CPD","CPD","6-4","6-4","6-4","6-4","CPD","CPD")

#draw the plot
ggplot(pca_data,aes(x=pc1,y=pc2,color=DamageType,shape=Organism))+geom_point(size=7,alpha=0.7) + 
    ggforce::geom_mark_ellipse(aes(label = CellLine, group = CellLine),color="Gray",con.colour = "gray",con.type="straight",label.colour = "Black",expand=0.03) + 
    theme_minimal() + 
    xlab("PC1 (%63.96)") +
    ylab("PC2 (%13.95)") + 
    xlim(0.20,0.38) +
    scale_y_continuous() + 
    theme(
        axis.title.x = element_text(size = 14, vjust=-5),
        axis.text.x = element_text(size = 14, vjust=-1),
        axis.title.y = element_text(size = 14, vjust=5),
        axis.text.y = element_text(size = 14, vjust=-1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        plot.margin=unit(c(1,1,1.5,1.2),"cm")) 