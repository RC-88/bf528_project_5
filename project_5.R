
library(BiocManager)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(ComplexHeatmap)
library(purrr)

##output results to group folder
file_path = "/Users/reinachau/Downloads/"

##programmer part
fpkm_tracking = read.table(file.path(file_path, "genes.fpkm_tracking"), header=T) %>% 
  dplyr::filter(FPKM > 0)

min(fpkm_tracking$FPKM)
max(fpkm_tracking$FPKM)

filtered_fpkm_tracking <- fpkm_tracking %>% 
  dplyr::filter(FPKM >= 1)

filtered_fpkm_tracking %>% 
  ggplot(aes(log(FPKM))) + 
  geom_histogram(fill='gray', color="white") +
  theme_classic()

###Analyst Part####

# part 6.1 - PO versus Adult"
gene_exp_diff = read.table(file.path(file_path, "gene_exp.diff"), header=T)
  
top10_diff_genes <- gene_exp_diff %>% 
  arrange(q_value) %>% 
  slice(1:10) %>% 
  select(gene, log2_foldchange=log2.fold_change., p_value, q_value)

write.csv(top10_diff_genes, file.path(file_path, "top_10_diff_exp_genes.csv"), row.names=F)

# part 6.2 - histogram of log2foldchange
fold_change_plot <- gene_exp_diff %>% 
  filter(log2.fold_change. >= -8 & log2.fold_change. <= 8) %>% 
  ggplot(aes(log2.fold_change.)) + 
  geom_histogram(bins=30, fill='gray', color="white") +
  scale_x_continuous(breaks=-8:8, label=-8:8) +
  xlab("Log2 Fold Change") + 
  theme_classic()

ggsave(filename="log2_fold_change.png", plot=fold_change_plot, path=file_path, width=4, height=5, units="in")

# part 6.3 - subset with significant genes only (2139 genes)
gene_significant <- gene_exp_diff %>% 
  filter(significant == "yes")

write.csv(gene_significant, file.path(file_path, "significant_genes.csv"), row.names=F)

# part 6.4 - histogram of log2foldchange of significant genes only
sig_fold_change_plot <- gene_significant %>% 
  filter(log2.fold_change. >= -8 & log2.fold_change. <= 8) %>% 
  ggplot(aes(log2.fold_change.)) + 
  geom_histogram(bins=30, fill='gray', color="white") +
  scale_x_continuous(breaks=-8:8, label=-8:8) +
  xlab("Log2 Fold Change") + 
  theme_classic()

ggsave(filename="sig_log2_fold_change.png", plot=sig_fold_change_plot, path=file_path, width=4, height=5, units="in")

# part 6.5 & 6.6 - up-regulated genes from significant genes only (n=1084)
up_regulated_genes <- gene_significant %>% 
  filter(log2.fold_change. > 0) %>% 
  select(gene)

write.table(up_regulated_genes, file.path(file_path, "up_regulated_genes.txt"), col.names=F, row.names=F, quote=F)

# part 6.5 & 6.6 - down-regulated genes from significant genes only (n=1055)
down_regulated_genes <- gene_significant %>% 
  filter(log2.fold_change. < 0) %>% 
  select(gene)

write.table(down_regulated_genes, file.path(file_path, "down_regulated_genes.txt"), col.names=F, row.names=F, quote=F)

# part 6.7 - DAVID Functional Annotation Clustering groups gene sets 


## Biologist

# part 7.1 
sample_path = "/project/bf528/project_2/data/samples/"

samples = c("Ad_1","Ad_2", "P0_1", "P0_2", "P4_1", "P4_2", "P7_1", "P7_2")

# create fpkm matrix for all samples
sample_fpkm_matrix <- NULL

for(s in seq_along(samples)){
  #s=1;
  if(samples[s] %in% "P0_1"){
    sample_genes_fpkm_tracking <- read.table("/projectnb/bf528/users/wheeler/project_2/P0_1_cufflinks/genes.fpkm_tracking", header=T) %>% 
      select(tracking_id, gene=gene_short_name, locus, FPKM) %>% 
      rename(!!paste0(samples[s]) := FPKM) %>% 
      arrange(tracking_id)
  }else{
    sample_genes_fpkm_tracking <- read.table(file.path(sample_path, samples[s], "genes.fpkm_tracking"), header=T) %>% 
      select(tracking_id, gene=gene_short_name, locus, FPKM) %>% 
      rename(!!paste0(samples[s]) := FPKM) %>% 
      arrange(tracking_id)
  }
  
  if(is.null(sample_fpkm_matrix)){
    sample_fpkm_matrix <- sample_genes_fpkm_tracking
  }else{
    sample_fpkm_matrix <- sample_fpkm_matrix %>% inner_join(sample_genes_fpkm_tracking, by=c("tracking_id", "gene", "locus"))
  }
  
}

## PLOT SAMPLES
plot_samples <- c("P0", "P4", "P7", "Ad")

## SARCOMERE PLOT ####
title <- "Sarcomere"
gene_samples <- c("Pdlim5", "Pygm", "Myoz2", "Des", "Csrp3", "Tcap", "Cryab") 
gene_colors <- brewer.pal(length(gene_samples), "Dark2")

fpkm_plot <- NULL;

for(s in seq_along(plot_samples)){
  #s=1;
  coln <- grep(plot_samples[s], colnames(sample_fpkm_matrix), value=T)
  fpkm <- sample_fpkm_matrix %>% 
    select(all_of(c("gene", coln))) %>% 
    rename(sample1 := !!coln[1], sample2 := !!coln[2])  %>% 
    filter(gene %in% gene_samples) %>% 
    mutate(sample_avg = (sample1 + sample2)/2, label=plot_samples[s], value=s)
  
  fpkm_plot <- rbind(fpkm_plot, fpkm)
}

# Order the gene names
fpkm_plot$gene <- factor(fpkm_plot$gene, levels=gene_samples, ordered=is.ordered(gene_samples))

# line plot
sarcomere_plot <- fpkm_plot %>% 
  ggplot() + 
  geom_line(aes(x=value, y=sample_avg, color=factor(gene))) +
  geom_point(aes(x=value, y=sample_avg, shape=factor(gene))) +
  scale_color_manual("Gene:", values=gene_colors) +
  scale_shape_manual("Gene:", values=1:length(gene_samples)) +
  scale_x_continuous(breaks=1:4, labels=plot_samples) +
  scale_y_continuous(limits=c(0, 1400), breaks=seq(0, 1400, by=200), labels=seq(0, 1400, by=200)) +
  theme_classic() + 
  ggtitle(title) +
  xlab("") +
  ylab("FPKM") +
  theme(plot.title=element_text(hjust=0.5))

ggsave(filename="sarcomere_plot.png", plot=sarcomere_plot, path=file_path, width=5, height=4, units="in")

## MITOCHONDRIA PLOT ####
title <- "Mitochrondria"
gene_samples <- c("Mpc1", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh") 
gene_colors <- brewer.pal(length(gene_samples), "Dark2")

fpkm_plot <- NULL

for(s in seq_along(plot_samples)){
  #s=1;
  coln <- grep(plot_samples[s], colnames(sample_fpkm_matrix), value=T)
  fpkm <- sample_fpkm_matrix %>% 
    select(all_of(c("gene", coln))) %>% 
    rename(sample1 := !!coln[1], sample2 := !!coln[2])  %>% 
    filter(gene %in% gene_samples) %>% 
    mutate(sample_avg = (sample1 + sample2)/2, label=plot_samples[s], value=s)
  
  fpkm_plot <- rbind(fpkm_plot, fpkm)
}

fpkm_plot$gene <- factor(fpkm_plot$gene, levels=gene_samples, ordered=is.ordered(gene_samples))

mitochrondria_plot <- fpkm_plot %>% 
  ggplot() + 
  geom_line(aes(x=value, y=sample_avg, color=factor(gene))) +
  geom_point(aes(x=value, y=sample_avg, shape=factor(gene))) +
  scale_color_manual("Gene:", values=gene_colors) +
  scale_shape_manual("Gene:", values=1:length(gene_samples)) +
  scale_x_continuous(breaks=1:4, labels=plot_samples) +
  scale_y_continuous(limits=c(0, 300), breaks=c(0, 100, 200, 300), labels=c(0, 100, 200, 300)) +
  theme_classic() + 
  ggtitle(title) +
  xlab("") +
  ylab("FPKM") +
  theme(plot.title=element_text(hjust=0.5))

ggsave(filename="mitochrondria_plot.png", plot=mitochrondria_plot, path=file_path, width=5, height=4, units="in")

## CELL CYCLE PLOT ####
title <- "Cell Cycle"
gene_samples <- c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "E2f1", "Cdc27", "Bora", "Cdc45", "Rad51", "Aurkb", "Cdc23") 
gene_colors <- colorRampPalette(brewer.pal(8, "Dark2"))

fpkm_plot <- NULL

for(s in seq_along(plot_samples)){
  #s=1;
  coln <- grep(plot_samples[s], colnames(sample_fpkm_matrix), value=T)
  fpkm <- sample_fpkm_matrix %>% 
    select(all_of(c("gene", coln))) %>% 
    rename(sample1 := !!coln[1], sample2 := !!coln[2])  %>% 
    filter(gene %in% gene_samples) %>% 
    mutate(sample_avg = (sample1 + sample2)/2, label=plot_samples[s], value=s)
  
  fpkm_plot <- rbind(fpkm_plot, fpkm)
}

fpkm_plot$gene <- factor(fpkm_plot$gene, levels=gene_samples, ordered=is.ordered(gene_samples))

cell_cycle_plot <- fpkm_plot %>% 
  ggplot() + 
  geom_line(aes(x=value, y=sample_avg, color=factor(gene))) +
  geom_point(aes(x=value, y=sample_avg, shape=factor(gene))) +
  scale_color_manual("Gene:", values=gene_colors(length(gene_samples))) +
  scale_shape_manual("Gene:", values=1:length(gene_samples)) +
  scale_x_continuous(breaks=1:4, labels=plot_samples) +
  scale_y_continuous(limits=c(0, 40), breaks=seq(0, 40, by=10), labels=seq(0, 40, by=10)) +
  theme_classic() + 
  ggtitle(title) +
  xlab("") +
  ylab("FPKM") +
  theme(plot.title=element_text(hjust=0.5))

ggsave(filename="cell_cycle_plot.png", plot=cell_cycle_plot, path=file_path, width=5, height=4, units="in")

## part 7.3 ####

set.seed(12345678)

# Top 1000 differential expressed genes
gene_labels <- c(
  "Pdlim5", "Pygm", "Myoz2", "Des", "Csrp3", "Tcap", "Cryab",
  "Mpc1", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh",
  "Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "E2f1", "Cdc27", "Bora", "Cdc45", "Rad51", "Aurkb", "Cdc23"
) 

# Top 1000 differential expressed genes
top_1000_sig_genes <-  gene_significant %>% 
  arrange(q_value) %>% 
  slice(1:1000) %>% 
  select(gene) %>% 
  unlist()

# filter fpkm matrix by top 1000 significant genes
sig_genes <- sample_fpkm_matrix %>% 
  filter(gene %in% top_1000_sig_genes) %>% 
  select(gene) %>% 
  unlist()

## fpkm matrix with the samples only
mat <- sample_fpkm_matrix %>% 
  filter(gene %in% sig_genes) %>% 
  mutate(row_id=paste0(gene, "_", locus)) %>% 
  column_to_rownames(var="row_id") %>% 
  select(-tracking_id, -gene, -locus) %>% 
  as.matrix()

## scale the data to Z-scores for heatmap
heat <- zFPKM(mat)

## Remove infinite values
heat <- heat[!is.infinite(rowSums(heat)),]

## set colour scheme and choose breaks
myCol <- colorRampPalette(c('blue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

## colum color
col.sample <- brewer.pal(length(colnames(heat)), 'Dark2')

column_df <- data.frame(
  Sample = colnames(heat),
  stringsAsFactors = FALSE
)

column_col <- list()

for(l in seq_along(column_df$sample)){
  column_col[[l]] <- column_df$col[l]
}

names(column_col) <- column_df$sample

ha_column = HeatmapAnnotation(df=column_df, col=column_col)

genelabels <- rowAnnotation(
  Genes = anno_mark(
    at = which(sig_genes %in% gene_labels),
    which = "row",
    labels = sig_genes[which(sig_genes %in% gene_labels)],
    labels_gp = gpar(fontsize = 10),
    padding = unit(1, "mm"))
)

ht1 = Heatmap(
  heat, 
  name = 'Z-score',
  col = colorRamp2(myBreaks, myCol),

  # parameters for the colour-bar that represents gradient of expression
  heatmap_legend_param = list(
    color_bar = 'continuous',
    legend_direction = 'horizontal',
    legend_width = unit(8, 'cm'),
    legend_height = unit(5.0, 'cm'),
    title_position = 'topcenter',
    title_gp=gpar(fontsize = 12, fontface = 'bold'),
    labels_gp=gpar(fontsize = 12)
  ),
  
  # column (sample) parameters
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  column_title = '',
  column_title_side = 'bottom',
  column_title_gp = gpar(fontsize = 12),
  column_title_rot = 0,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 10),
  column_names_max_height = unit(10, 'cm'),
  column_dend_height = unit(25,'mm'),

  # row (gene) parameters
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  #row_title = 'Statistically significant genes',
  row_title_side = 'left',
  row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
  row_title_rot = 90,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  row_names_side = 'left',
  row_dend_width = unit(25,'mm'),
  
  # cluster methods for rows and columns
  clustering_distance_columns = function(x) dist(x),
  clustering_method_columns = 'ward.D2',
  clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
  clustering_method_rows = 'ward.D2',
  top_annotation = ha_column)

ht_list = ht1 + genelabels

draw(
  ht_list, 
  heatmap_legend_side = 'bottom',
  annotation_legend_side = 'right',
  row_sub_title_side = 'left'
)


