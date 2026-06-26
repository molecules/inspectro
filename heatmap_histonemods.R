#!/bin/env -S Rscript --vanilla
#SBATCH --cpus-per-task 24
#SBATCH --mem-per-cpu 4G
#SBATCH --exclude cc011        # unreliable node
#SBATCH --mail-type END,FAIL   # get an email when this finishes or fails

library(rtracklayer)
library(tidyverse) #load dplyr, ggplot2, 

#Read in arguments from the command line, starting after this script's name
args <- commandArgs(trailingOnly = TRUE)

cluster_file <- args[1]
eigvecs_file <- args[2]
bigwig_file <- args [3]
mod_name <- args[4]

clusterbasename <- basename(cluster_file)
base_out_name <- strsplit(clusterbasename,"\\.") [[1]][1]

print("we made it to 21")

bw <- import(format = "BigWig", con = bigwig_file)

bins_50kb <- tileGenome(seqlengths(bw), tilewidth = 50000, cut.last.tile.in.chrom=TRUE)

bin_result <- binnedAverage(bins_50kb, numvar = coverage(bw, weight = "score"), varname = "avg_score")
print("we made it to 28!!")

data_ranges_df <- bin_result %>% as.data.frame()

df <- read.csv(eigvecs_file)

clust <- read.table(cluster_file, header=TRUE)
k9 <- clust[, c("chrom", "start", "end", "kmeans_sm9") ]
combo <- merge(df, k9)

small <- combo[,c("chrom","start","end","GC","arm","armlen","centel","centel_abs","is_bad","kmeans_sm9","E0","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10")]
longer <- small %>% pivot_longer(cols=c("E0","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10"), values_to = "value")
bin_result
head(data_ranges_df)
longer[1:20, ]
longer[1:40, ]
longer[30:40, ]
longer[3010:3040, ]
longer2 <- longer
longer2 <- longer %>% mutate(start = start +1 )
longer2
big_combo <- merge(longer2, data_ranges_df)
ls()
head(small)
small2 <- small %>% mutate(start = start +1)
big_combo <- merge(small2, data_ranges_df)
head(big_combo)
head(data_ranges_df)
data_ranges_df2 <- data_ranges_df %>% mutate(chrom = seqnames)
print("we made it to 57!!")
# Did we want all?
big_combo <- merge(small2, data_ranges_df2, all = TRUE)
big_combo <- merge(small2, data_ranges_df2, all = FALSE)

plot_Histone <- ggplot(big_combo, aes(x=centel, y=1, color = avg_score)) +
                geom_tile() +
                scale_color_viridis_c() +
                facet_grid(.~kmeans_sm9)
print("we made it to 66!!")
pdf_filename <-paste0(base_out_name,".", mod_name,".kmeanssm9.pdf")
pdf(pdf_filename, width =18, height =2)
print(plot_Histone)
dev.off()
print("we finished!!")