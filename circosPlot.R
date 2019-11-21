#!/usr/bin/Rscript


#load necessay packages
library(circlize)
library(data.table)
library(viridis)
library(scales)
library(svglite)
library(ComplexHeatmap)
library(yarrr)

#Load chr lengths for circos plot to initialize to represent two genomes as one
chrLengths = fread("chrLengths.txt")
#         chr start       end
# 1: LemurChr1     1 114275543
# 2: LemurChr2     1 141296404
# 3: LemurChr3     1 119751181
#               ...
#54:  HumanChr3     1 198022430
#55:  HumanChr2     1 243199373
#56:  HumanChr1     1 249250621

# Load links for overlapping regions between two genomes as the result of whole genome-genome alignment
links = dplyr::select(fread("bedf.bed"),c(1,2,3,4,5,6,9))
# Its a bed file where 1,2,3. columns represent lemur and 4,5,6. columns represent aligned regions of human:
#chr1	17703769	17709999	HumanChr3	25469161	25475344	0	+	88.472443	+	5738
#chr1	67213689	67219159	HumanChr3	173513208	173518761	0	+	84.097366	+	5470
#chr1	76715866	76721084	HumanChr3	183880975	183886195	0	+	85.458992	+	5218

#Choose the alignment percentages for colorization
perc = data.frame(links$V9)

#choose and correct the column names of the links starts and ends (link1 (lemur) and link2 (human))
link1 = links[,c(1,2,3)]
link2 = links[,c(4,5,6)]
#reformat the links chr names to match with chr names in chrLengths
link1$V1 = gsub("chr","LemurChr",link1$V1)
link2$V4 = gsub("chr","HumanChr",link2$V4)
colnames(link1) = c("chr","start","end")
colnames(link2) = c("chr","start","end")
#manually define sector labels to print on circos plot
sector_label=c(1:32,"x","x",22:1)

#either start a png session
#png("test_heatmap.png", width = 16, height = 16, units = 'in', res = 1200)

#or start an svg session for plotting the circos plot
svglite("test_svg.svg", width = 16, height = 16)

#clear previous circos session if present
circos.clear()
#define the graphical parameters for circos session
circos.par(gap.after = c(rep(0.5,32),15, rep(0.5,22),15), start.degree = 90,cell.padding = c(0, 0, 0, 0))
#initialize a null genome with defined 56 chr (sum of lemur chrs and human chrs)
circos.genomicInitialize(chrLengths, plotType = NULL)

#Load the RPKM values of each analysis
res = fread("./counts/final_go/res.tsv")
colnames(res)=c("el.HumanXrRPKM", "el.HumanXrCount", "el.LemurXrRPKM", "el.LemurXrCount", "el.HumanSimRPKM", "el.HumanSimCount", "el.LemurSimRPKM", "el.LemurSimCount", "el.HumanRNARPKM", "el.HumanRNACount", "el.LemurRNARPKM", "el.LemurRNACount", "el.Genic", "el.Identity", "el.LemurChr", "el.LemurStart", "el.LemurEnd", "el.HumanChr", "el.HumanStart", "el.HumanEnd")

######################################################
# Human side of the circos plot shows xr-seq log2(rpkm) values of human and lemur, both divided by theirs simulated XR-seq representatives
# Lemur side of the circos plot shows rna-seq log2(rpkm) values

lemur_rna_lemur = cbind(data.frame(res$el.LemurChr),res$el.LemurStart,res$el.LemurEnd,res$el.LemurRNARPKM)
colnames(lemur_rna_lemur) = c("chr","start","end","value1")
lemur_rna_lemur = subset(lemur_rna_lemur, value1 != 0)
lemur_rna_lemur$value1 = log2(lemur_rna_lemur$value1)
lemur_rna_lemur$chr = gsub("chr","LemurChr",lemur_rna_lemur$chr)

lemur_rna_human = cbind(data.frame(res$el.LemurChr),res$el.LemurStart,res$el.LemurEnd,res$el.HumanRNARPKM)
colnames(lemur_rna_human) = c("chr","start","end","value1")
lemur_rna_human = subset(lemur_rna_human, value1 != 0)
lemur_rna_human$value1 = log2(lemur_rna_human$value1)
lemur_rna_human$chr = gsub("chr","LemurChr",lemur_rna_human$chr)

human_xr_over_sim_lemur = cbind(data.frame(res$el.HumanChr),res$el.HumanStart,res$el.HumanEnd,(res$el.LemurXrRPKM/res$el.LemurSimRPKM))
colnames(human_xr_over_sim_lemur) = c("chr","start","end","value1")
human_xr_over_sim_lemur=subset(human_xr_over_sim_lemur[Reduce(`&`, lapply(human_xr_over_sim_lemur, function(x) !is.na(x)  & is.finite(x))),], value1!=0)
human_xr_over_sim_lemur$value1 = log2(human_xr_over_sim_lemur$value1)
human_xr_over_sim_lemur$chr = gsub("chr","HumanChr",human_xr_over_sim_lemur$chr)

human_xr_over_sim_human = cbind(data.frame(res$el.HumanChr),res$el.HumanStart,res$el.HumanEnd,(res$el.HumanXrRPKM/res$el.HumanSimRPKM))
colnames(human_xr_over_sim_human) = c("chr","start","end","value1")
human_xr_over_sim_human=subset(human_xr_over_sim_human[Reduce(`&`, lapply(human_xr_over_sim_human, function(x) !is.na(x)  & is.finite(x))),],value1!=0)
human_xr_over_sim_human$value1 = log2(human_xr_over_sim_human$value1)
human_xr_over_sim_human$chr = gsub("chr","HumanChr",human_xr_over_sim_human$chr)

######################################################

#define a color scale for all xr-seq log2(rpkm) values
xr_all = rbind(human_xr_over_sim_human,human_xr_over_sim_lemur)
xr_scale = colorRamp2(unique(sort(xr_all$value1)),viridis(length(unique(xr_all$value1))))

#draw the genomics tracks of xr-seq analysis
circos.genomicTrack(human_xr_over_sim_human, stack = TRUE, 
                    panel.fun = function(region, value, ...) {
                        col = xr_scale(value[[1]])
                        circos.genomicRect(region, value, border = col, ...)
                    },bg.border=NA, track.index = 1)
circos.genomicTrack(human_xr_over_sim_lemur, stack = TRUE, 
                    panel.fun = function(region, value, ...) {
                        col = xr_scale(value[[1]])
                        circos.genomicRect(region, value, border = col, ...)
                    },bg.border=NA, track.index = 2)

#define a color scale for all rna-seq log2(rpkm) values
rna_all = rbind(lemur_rna_human,lemur_rna_lemur)
rna_scale = colorRamp2(unique(sort(rna_all$value1)),viridis(length(unique(rna_all$value1))))

#draw the genomics tracks of rna-seq analysis
circos.genomicTrack(lemur_rna_human, stack = TRUE, 
                    panel.fun = function(region, value, ...) {
                        col = rna_scale(value[[1]])
                        circos.genomicRect(region, value, border = col, ...)
                    },bg.border=NA, track.index = 1)
circos.genomicTrack(lemur_rna_lemur, stack = TRUE, 
                    panel.fun = function(region, value, ...) {
                        col = rna_scale(value[[1]])
                        circos.genomicRect(region, value, border = col, ...)
                    },bg.border=NA, track.index = 2)



######################################################


#draw the legends of rna-seq and xr-seq log2(rpkm) values
legend_xr = Legend(at = c(min(xr_all$value1),median(xr_all$value1),max(xr_all$value1)), 
                     col_fun = xr_scale, 
                     title_position = "topleft", 
                     title = "XR - HeatMap",
                     title_gap = unit(4, "mm"),
                     title_gp = gpar(col=viridis(1),fontface="bold"),
                     grid_height = unit(20, "mm"),
                     grid_width = unit(5, "mm"),
                     legend_gp = gpar(fill),
                     direction = "vertical")

legend_rna = Legend(at = c(min(rna_all$value1),median(rna_all$value1),max(rna_all$value1)), 
                     col_fun = rna_scale, 
                     title_position = "topleft", 
                     title = "log2 RNA - RPKM",
                     title_gap = unit(4, "mm"),
                     title_gp = gpar(col=viridis(1),fontface="bold"),
                     grid_height = unit(20, "mm"),
                     grid_width = unit(5, "mm"),
                     legend_gp = gpar(fill),
                     direction = "vertical")

draw(legend_xr, x = unit(30,"mm"), y = unit(15,"mm"),just = "bottom")
draw(legend_rna, x = unit(376,"mm"), y = unit(15,"mm"),just = "bottom")

#draw the chr track and label by sector_label
#lemur chr track will be a constant color
colors_chr = viridis(56)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    index = CELL_META$sector.numeric.index
    chr = sector_label[index]
    chr_colors = colors_chr[index]
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    
    if (index <= 33) {
    circos.rect(xlim[1], 0, xlim[2], 0.45, col = viridis(1))
    } else {circos.rect(xlim[1], 0, xlim[2], 0.45, col = chr_colors)}

    circos.text(mean(xlim), mean(ylim)+0.25, chr, cex = 0.375, col = viridis(1),
                facing = "inside", niceFacing = TRUE)
}, track.height = 0.05, bg.border = NA)

# cut the percentages of alingment to breaks between %80, %90 and %100 in order to make transparent with constant ratios
# where the low percentage alignments (%80-%90) will be represented lightly
Alpha = as.vector(cut(perc$links.V9, breaks = c(80, 90, 100), labels = c(0.85,0.2), include.lowest = T))
# links will be colored depending on the color of chr track of human
hchr = chrLengths$chr[c(34:56)]
names(hchr) = viridis(56)[34:56]

# only for printing out to report how many links left to draw
i_perc = nrow(link1)/10

for (i in 1:nrow(link1)) {
    color = names(hchr[hchr==link2[i]$chr])
    color1 = transparent(color, as.numeric(Alpha[i]))
    sector.index1 = link2$chr[i]
    point1=c(link2$start[i],link2$end[i])
    sector.index2 =link1$chr[i]
    point2=c(link1$start[i],link1$end[i])
    if (i %% i_perc == 0){print(paste0("%",i/i_perc*10," completed"))}
    circos.link(
        sector.index1 = sector.index1,
        point1 = point1,
        sector.index2 =sector.index2, 
        point2 = point2,
        col = color1,
        h.ratio=0.85
    )
}

#clear the circos session and save the plot
circos.clear()
dev.off()
