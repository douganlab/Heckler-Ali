# CDK4Analysis.R
# Copyright (C) 2020 Lestat Ali
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(Seurat)
library(cowplot)

source("CDK4AnalysisHelpers.R")
source("CDK4PlotHelpers.R")

healthy_ids <- paste0("S", 1:4)
other_patient_ids <- c("P616", "P673")
treated_patient_ids <- c("P119", "P560", "P671", "P687", "P823")
all_patient_ids <- c(treated_patient_ids, other_patient_ids)
all_subject_ids <- c(healthy_ids,
                     paste0(all_patient_ids, "_Pre"),
                     paste0(other_patient_ids, "_Post"),
                     paste0(treated_patient_ids, "_Post"))

#### Individual Sample Loading ####

s1_seurat <- create.seurat.object.for.healthy.subject("S1")
s1_seurat <- annotate.with.tcr.data(s1_seurat, "S1")

s2_seurat <- create.seurat.object.for.healthy.subject("S2")
s2_seurat <- annotate.with.tcr.data(s2_seurat, "S2")

s3_seurat <- create.seurat.object.for.healthy.subject("S3")
s3_seurat <- annotate.with.tcr.data(s3_seurat, "S3")

s4_seurat <- create.seurat.object.for.healthy.subject("S4")
s4_seurat <- annotate.with.tcr.data(s4_seurat, "S4")

p119_seurat <- create.seurat.object.for.patient("P119")
p119_seurat <- annotate.with.tcr.data(p119_seurat, "P119")

p560_seurat <- create.seurat.object.for.patient("P560")
p560_seurat <- annotate.with.tcr.data(p560_seurat, "P560")

p616_seurat <- create.seurat.object.for.patient("P616")
p616_seurat <- annotate.with.tcr.data(p616_seurat, "P616")

p671_seurat <- create.seurat.object.for.patient("P671")
p671_seurat <- annotate.with.tcr.data(p671_seurat, "P671")

p673_seurat <- create.seurat.object.for.patient("P673")
p673_seurat <- annotate.with.tcr.data(p673_seurat, "P673")

p687_seurat <- create.seurat.object.for.patient("P687")
p687_seurat <- annotate.with.tcr.data(p687_seurat, "P687")

p823_seurat <- create.seurat.object.for.patient("P823")
p823_seurat <- annotate.with.tcr.data(p823_seurat, "P823")

s1_seurat <- RenameCells(s1_seurat,  add.cell.id = "S1")
s2_seurat <- RenameCells(s2_seurat,  add.cell.id = "S2")
s3_seurat <- RenameCells(s3_seurat,  add.cell.id = "S3")
s4_seurat <- RenameCells(s4_seurat,  add.cell.id = "S4")
p119_seurat <- RenameCells(p119_seurat,  add.cell.id = "P119")
p560_seurat <- RenameCells(p560_seurat,  add.cell.id = "P560")
p616_seurat <- RenameCells(p616_seurat,  add.cell.id = "P616")
p671_seurat <- RenameCells(p671_seurat,  add.cell.id = "P671")
p673_seurat <- RenameCells(p673_seurat,  add.cell.id = "P673")
p687_seurat <- RenameCells(p687_seurat,  add.cell.id = "P687")
p823_seurat <- RenameCells(p823_seurat,  add.cell.id = "P823")

#### Sample Integration ####

cd45dp_all <- list(s1_seurat, 
                   s2_seurat, 
                   s3_seurat, 
                   s4_seurat,
                   p119_seurat,
                   p560_seurat,
                   p616_seurat,
                   p671_seurat,
                   p673_seurat, 
                   p687_seurat, 
                   p823_seurat)

sample_treatment_id_map <- c(rep("Healthy", length(healthy_ids)),
                             rep("Cancer-NoCDK4/6i", length(all_patient_ids) + length(other_patient_ids)),
                             rep("Cancer-PostCDK4/6i", length(treated_patient_ids)))
names(sample_treatment_id_map) <- all_subject_ids

treatments <- sapply(Cells(cd45dp_full_int), function(cellname) {
  matching_id <- all_subject_ids[sapply(all_subject_ids, grepl, cellname)]
  return(sample_treatment_id_map[matching_id])
})
names(treatments) <- NULL

cd45dp_anchors <- FindIntegrationAnchors(object.list = cd45dp_all, dims = 1:15, scale = F, anchor.features = 2000)
cd45dp_full_int <- IntegrateData(anchorset = cd45dp_anchors, dims = 1:15)
cd45dp_full_int <- AddMetaData(cd45dp_full_int, treatments, col.name = "treatment")

DefaultAssay(cd45dp_full_int) <- "integrated"
cd45dp_full_int <- ScaleData(cd45dp_full_int)
cd45dp_full_int <- RunPCA(cd45dp_full_int)
cd45dp_full_int <- RunUMAP(cd45dp_full_int, dims = 1:14)
cd45dp_full_int <- FindNeighbors(cd45dp_full_int, dims = 1:14)
cd45dp_full_int <- FindClusters(cd45dp_full_int, resolution = 0.4)

# Exclude contaminating cluster 11 from further analysis
cd45dp_no11_int <- cd45dp_full_int[, Idents(cd45dp_full_int) != 11]

cluster_labels <- c("0" = "CD45RA+ naive",
                    "1" = "naive transitional",
                    "2" = "MCHII high",
                    "3" = "CTLs",
                    "4" = "CD45RO+ memory",
                    "5" = "transitional",
                    "6" = "GATA3+",
                    "7" = "NK-like effectors",
                    "8" = "naive",
                    "9" = "Trm",
                    "10" = "IFN-activated")

# Figure 3B
draw.cluster.umap(cd45dp_no11_int)

# Supplemental Figure
draw.clonal.umap(cd45dp_no11_int)

# Figure 3D
draw.monotone.umap(cd45dp_no11_int[, cd45dp_no11_int$treatment == "Healthy"], "grey20", T)
draw.monotone.umap(cd45dp_no11_int[, cd45dp_no11_int$treatment == "Cancer-NoCDK4/6i"], "#009E73", T)
draw.monotone.umap(cd45dp_no11_int[, cd45dp_no11_int$treatment == "Cancer-PostCDK4/6i"], "#E69F00", T)
all_cluster_ids <- c("1", "2", "3", "5E", "5M", "6", "7", "8", "9", "10")
pole_cluster_ids <- c("0", "4")

cluster_percents_by_patient <- data.frame("ClusterID" = rep(all_cluster_ids, length(all_subject_ids)),
                                          "PatientID" = unlist(lapply(all_subject_ids, rep, length(all_cluster_ids))),
                                          "Percentage" = rep(0, length(all_cluster_ids) * length(all_subject_ids)),
                                          stringsAsFactors = F)

cluster_percents_mat <- matrix(data = 0, nrow = length(all_subject_ids), ncol = length(all_cluster_ids))
rownames(cluster_percents_mat) <- all_subject_ids
colnames(cluster_percents_mat) <- all_cluster_ids
for (i in 1:nrow(cluster_percents_mat)) {
  for (j in 1:ncol(cluster_percents_mat)) {
    current_patient <- rownames(cluster_percents_mat)[i]
    current_cluster <- colnames(cluster_percents_mat)[j]
    current_cluster_count <- length(which(startsWith(Cells(cd45dp_no11_int), current_patient) &
                                            Idents(cd45dp_no11_int) == current_cluster))
    all_clusters_total <- length(which(startsWith(Cells(cd45dp_no11_int), current_patient) &
                                         !(Idents(cd45dp_no11_int) %in% pole_cluster_ids)))
    percentage <- (100 * current_cluster_count / all_clusters_total)
    cluster_percents_mat[i, j] <- round(percentage, 3)
  }
}
write.csv(cluster_percents_mat, file = "tables/cluster_percentages_by_patient_nonpole_updated.csv")

# Figure 3C
markers <- c("TCF7", 
             "JUN",
             "IL7R",
             "HLA-DRB1",
             "GZMH",
             "AQP3",
             "DUSP2",
             "CCL5",
             "GATA3",
             "NKG7",
             "GNLY",
             "S100A9",
             "ZNF683",
             "IRF7")
draw.violin.sandwich(cd45dp_no11_int, markers, T)

# Patient Contribution to Clusters (Figure 3E)
all_cluster_ids <- 0:10
cluster_percents_by_patient <- data.frame("ClusterID" = rep(all_cluster_ids, length(all_subject_ids)),
                                          "PatientID" = unlist(lapply(all_subject_ids, rep, length(all_cluster_ids))),
                                          "Percentage" = rep(0, length(all_cluster_ids) * length(all_subject_ids)),
                                          stringsAsFactors = F)

for (i in 1:nrow(cluster_percents_by_patient)) {
  current_patient <- cluster_percents_by_patient[i, "PatientID"]
  current_cluster <- cluster_percents_by_patient[i, "ClusterID"]
  current_cluster_count <- length(which(startsWith(Cells(cd45dp_no11_int), current_patient) &
                                          cd45dp_no11_int$seurat_clusters == current_cluster))
  all_clusters_total <- length(which(startsWith(Cells(cd45dp_no11_int), current_patient)))
  percentage <- (100 * current_cluster_count / all_clusters_total)
  cluster_percents_by_patient$Percentage[i] <- percentage
}

cluster_percents_by_patient$ClusterID <- factor(cluster_percents_by_patient$ClusterID)
cluster_percents_by_patient$PatientID <- factor(cluster_percents_by_patient$PatientID, 
                                                levels = c("P616_Pre", 
                                                           "P616_Post",
                                                           "P673_Pre",
                                                           "P673_Post",
                                                           "P119_Pre", 
                                                           "P119_Post",
                                                           "P560_Pre",
                                                           "P560_Post",
                                                           "P671_Pre",
                                                           "P671_Post",
                                                           "P687_Pre",
                                                           "P687_Post",
                                                           "P823_Pre",
                                                           "P823_Post",
                                                           "S1",
                                                           "S2",
                                                           "S3",
                                                           "S4"))

contribution_plot <- ggplot(data = cluster_percents_by_patient, aes(x = PatientID, y = Percentage, fill = ClusterID)) +
  geom_col() + 
  scale_fill_manual(values = intg_clusters_colors) +
  scale_y_continuous(expand = expand_scale(0)) + 
  theme.lra() + 
  theme(panel.grid = element_blank(), 
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 10, hjust = 1, angle = 45))
contribution_plot


#### Sub-Clustering of Cluster 5 (Transitional Cells) #### 

cd45dp_sub_int <- cd45dp_full_int[, Idents(cd45dp_full_int) == 5]
cd45dp_sub_int <- ScaleData(cd45dp_sub_int)
cd45dp_sub_int <- RunPCA(cd45dp_sub_int)
cd45dp_sub_int <- RunUMAP(cd45dp_sub_int, dims = 1:12, seed.use = 11)
cd45dp_sub_int <- FindNeighbors(cd45dp_sub_int, dims = 1:12, k.param = 300)
cd45dp_sub_int <- FindClusters(cd45dp_sub_int, resolution = 0.3)
cd45dp_sub_int <- RenameIdents(cd45dp_sub_int, "0" = "5E", "1" = "5M")

# Figure 4A
draw.cluster.umap(cd45dp_sub_int, cluster_colors = c("brown", "gold2"), save = F)

# Retrospectively label main Seurat object (Figure 4B)
cd45dp_no11_int <- SetIdent(cd45dp_no11_int, cells = WhichCells(cd45dp_sub_int, idents = "5E"), value = "5E")
cd45dp_no11_int <- SetIdent(cd45dp_no11_int, cells = WhichCells(cd45dp_sub_int, idents = "5M"), value = "5M")
draw.highlighted.umap(cd45dp_no11_int, highlight_cells = WhichCells(cd45dp_sub_int, idents = "5E"), highlight_color = "brown", T)
draw.highlighted.umap(cd45dp_no11_int, highlight_cells = WhichCells(cd45dp_sub_int, idents = "5M"), highlight_color = "gold2", T)

# Figure 4C
markers_5e_vs_5m <- FindMarkers(cd45dp_no11_int, ident.1 = "5E", ident.2 = "5M")
write.csv(markers_5e_vs_5m, file = "markers_5e_vs_5m.csv")

fate_markers <- c("LTB",
                  "FOS",
                  "NR4A2",
                  "DUSP2",
                  "ZFP36L2",
                  "JUNB",
                  "CSRNP1",
                  "IL7R",
                  "DUSP1",
                  "KLF6",
                  "HLA-DRB1",
                  "ANXA1",
                  "GZMB",
                  "FGFBP2",
                  "NKG7",
                  "LGALS1",
                  "GNLY",
                  "GZMH")

fate_exp <- AverageExpression(cd45dp_no11_int, assay = "RNA", features = fate_markers)
fate_exp_mat <- as.matrix(fate_exp$RNA[, c("4", "5M", "5E", "3", "7")])
fate_exp_mat <- cbind(fate_exp_mat, rowMeans(fate_exp_mat[, c("3", "7")]))
fate_exp_mat <- fate_exp_mat[, -c(4, 5)]
colnames(fate_exp_mat)[4] <- "Effector"
colnames(fate_exp_mat)[1] <- "Memory"
cluster5_exp_means <- rowMeans(fate_exp_mat[,c("5M", "5E")])
fate_exp_z_mat <- t(scale(t(fate_exp_mat), center = cluster5_exp_means))
fate_exp_df <- reshape2::melt(fate_exp_z_mat)
colnames(fate_exp_df) <- c("Gene", "ClusterID", "Expression")
fate_exp_df$ClusterID <- factor(fate_exp_df$ClusterID, levels = rev(c("Memory", "5M", "5E", "Effector")))

fate_heatmap_vertical <- ggplot(data = fate_exp_df, aes(x = ClusterID, y = Gene, fill = Expression)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "royalblue3", mid = "white", high = "red3", limits = c(-1.5, 1.5), oob = scales::squish, breaks = c(-1, 0, 1), name = "Expression\n(z-score)") +
  scale_x_discrete(expand = expand_scale(0)) +
  scale_y_discrete(expand = expand_scale(0)) +
  theme.lra() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 12),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "right",
        legend.direction = "vertical")
fate_heatmap_vertical

ggsave(fate_heatmap_vertical, file = "plots/fate_heatmap_updated_vertical.pdf", width = 5, height = 6.5, device = cairo_pdf)

#### MYC Score #### 

# Gene set obtained from MSigDB
hallmark_myc_genes <- read.csv("myc_target_genes_hallmark.csv", header = F)$V1
myc_genes <- intersect(rownames(cd45dp_no11_int$RNA), hallmark_myc_genes)

myc_gene_expr <- as.matrix(GetAssayData(cd45dp_no11_int$RNA)[myc_genes,])
myc_gene_expr_avg <- rowMeans(myc_gene_expr)
myc_gene_expr_sd <- rowSds(myc_gene_expr)
myc_gene_expr_z <- (myc_gene_expr - myc_gene_expr_avg) / myc_gene_expr_sd
myc_score <- colSums(myc_gene_expr_z)

sample_id_by_cell <- sapply(Cells(cd45dp_no11_int), function(cellname) {
  all_subject_ids[sapply(all_subject_ids, grepl, cellname)]
})

myc_score_df <- data.frame("SampleID" = sample_id_by_cell, 
                           "ClusterID" = Idents(cd45dp_no11_int),
                           "MycScore" = myc_score)
myc_score_avg_df <- aggregate(MycScore ~ ClusterID + SampleID, myc_score_df, mean)
myc_score_cluster5_df <- subset(myc_score_df, ClusterID %in% c("5E", "5M"))
myc_score_cluster5_df$ClusterID <- factor(myc_score_cluster5_df$ClusterID, levels = c("5E", "5M"))

cd45dp_sub_int_pre <- cd45dp_sub_int[, (grepl("Pre", Cells(cd45dp_sub_int)) |
                                          startsWith(Cells(cd45dp_sub_int), "P616_Post") |
                                          startsWith(Cells(cd45dp_sub_int), "P673_Post")) ]

cd45dp_sub_int_post <- cd45dp_sub_int[, (grepl("Post", Cells(cd45dp_sub_int)) &
                                           !(startsWith(Cells(cd45dp_sub_int), "P616_Post") |
                                               startsWith(Cells(cd45dp_sub_int), "P673_Post"))) ]

myc_score_cl5_pre_df <- myc_score_cluster5_df[Cells(cd45dp_sub_int_pre),]
myc_score_cl5_post_df <- myc_score_cluster5_df[Cells(cd45dp_sub_int_post),]

# Figure 4E
draw.gradient.umap(cd45dp_sub_int_pre, myc_score_cl5_pre_df$MycScore, save = T)
draw.gradient.umap(cd45dp_sub_int_post, myc_score_cl5_post_df$MycScore, save = T)

# Figure 4F
myc_score_mat <- matrix(data = 0, nrow = length(all_subject_ids), ncol = length(levels(myc_score_avg_df$ClusterID)))
rownames(myc_score_mat) <- all_subject_ids
colnames(myc_score_mat) <- levels(myc_score_avg_df$ClusterID)
for (i in 1:nrow(myc_score_mat)) {
  for (j in 1:ncol(myc_score_mat)) {
    myc_score_mat[i, j] <- mean(subset(myc_score_df, ClusterID == colnames(myc_score_mat)[j] &
                                         SampleID == rownames(myc_score_mat)[i])$MycScore)
  }
}
write.csv(myc_score_mat, file = "tables/myc_score_matrix_msigdb.csv")

#### Clonotype Sharing Grids #### 

p119_clonotypes <- get.unified.clonotyped.barcodes("P119-Pre", "P119-Post", p119_seurat)
p560_clonotypes <- get.unified.clonotyped.barcodes("P560-Pre", "P560-Post", p560_seurat)
p671_clonotypes <- get.unified.clonotyped.barcodes("P671-Pre", "P671-Post", p671_seurat)
p687_clonotypes <- get.unified.clonotyped.barcodes("P687-Pre", "P687-Post", p687_seurat)
p823_clonotypes <- get.unified.clonotyped.barcodes("P823-Pre", "P823-Post", p823_seurat)

treated_list <- list(p119_clonotypes,
                     p560_clonotypes,
                     p671_clonotypes,
                     p687_clonotypes,
                     p823_clonotypes)

# Normalization Matrix
pre_cluster_ids <- c("0", "1", "2", "3", "4", "5E", "5M", "6", "7", "8", "9", "10")
correction_factors <- rep(0, length(pre_cluster_ids))
names(correction_factors) <- pre_cluster_ids
for (pre_cluster_id in pre_cluster_ids) {
  current_subfactors <- lapply(treated_list, function(clonotypes) {
    unique_cluster_clonotypes <- unique(subset(clonotypes, ClusterID == pre_cluster_id)$Clonotype)
    length(which(startsWith(unique_cluster_clonotypes, "U_")))
  })
  current_factor <- sum(unlist(current_subfactors))
  correction_factors[pre_cluster_id] <- current_factor
}

mini_correction_mat <- matrix(correction_factors, nrow = length(pre_cluster_ids), ncol = length(pre_cluster_ids), byrow = F)

# Figure 4G
p119_pre_sharing_mat <- compute.pre.clonotypes(p119_clonotypes)
p560_pre_sharing_mat <- compute.pre.clonotypes(p560_clonotypes)
p671_pre_sharing_mat <- compute.pre.clonotypes(p671_clonotypes)
p687_pre_sharing_mat <- compute.pre.clonotypes(p687_clonotypes)
p823_pre_sharing_mat <- compute.pre.clonotypes(p823_clonotypes)
total_pre_sharing_mat <- (p119_pre_sharing_mat
                          + p560_pre_sharing_mat
                          + p671_pre_sharing_mat
                          + p687_pre_sharing_mat
                          + p823_pre_sharing_mat)

corrected_pre_sharing_mat <- round(100 * total_pre_sharing_mat / mini_correction_mat, 0)
plot.sharing.mat(corrected_pre_sharing_mat, prefix = "pre_norm_", high_color = "green4", limit_scale = 60, save = F)
plot.sharing.mat(p671_pre_sharing_mat, prefix = "pre_", high_color = "green4", limit_scale = 35, save = F)

# Supplmentary Figure 7B
p119_added_sharing_mat <- compute.added.clonotypes(p119_clonotypes)
p560_added_sharing_mat <- compute.added.clonotypes(p560_clonotypes)
p671_added_sharing_mat <- compute.added.clonotypes(p671_clonotypes)
p687_added_sharing_mat <- compute.added.clonotypes(p687_clonotypes)
p823_added_sharing_mat <- compute.added.clonotypes(p823_clonotypes)
total_added_sharing_mat <- (p119_added_sharing_mat
                            + p560_added_sharing_mat
                            + p671_added_sharing_mat
                            + p687_added_sharing_mat
                            + p823_added_sharing_mat)

plot.sharing.mat(total_added_sharing_mat, prefix = "added_", limit_scale = 49, save = F)

# Supplmentary Figure 7B
p119_removed_sharing_mat <- compute.removed.clonotypes(p119_clonotypes)
p560_removed_sharing_mat <- compute.removed.clonotypes(p560_clonotypes)
p671_removed_sharing_mat <- compute.removed.clonotypes(p671_clonotypes)
p687_removed_sharing_mat <- compute.removed.clonotypes(p687_clonotypes)
p823_removed_sharing_mat <- compute.removed.clonotypes(p823_clonotypes)
total_removed_sharing_mat <- (p119_removed_sharing_mat
                              + p560_removed_sharing_mat
                              + p671_removed_sharing_mat
                              + p687_removed_sharing_mat
                              + p823_removed_sharing_mat)

plot.sharing.mat(total_removed_sharing_mat, prefix = "removed_", limit_scale = -57, save = F)

# Figure 4I
p119_net_sharing <- compute.net.change.clonotypes(p119_clonotypes)
p560_net_sharing <- compute.net.change.clonotypes(p560_clonotypes)
p671_net_sharing <- compute.net.change.clonotypes(p671_clonotypes)
p687_net_sharing <- compute.net.change.clonotypes(p687_clonotypes)
p823_net_sharing <- compute.net.change.clonotypes(p823_clonotypes)
total_net_sharing <- (p119_net_sharing
                      + p560_net_sharing
                      + p671_net_sharing
                      + p687_net_sharing
                      + p823_net_sharing)

net_sharing_df <- data.frame("ClusterID" = names(total_net_sharing), "NetFlux" = total_net_sharing)
net_sharing_df$ClusterID <- factor(net_sharing_df$ClusterID, levels = pre_cluster_ids)
net_sharing_plot <- ggplot(net_sharing_df, aes(x = ClusterID, y = NetFlux, fill = ClusterID)) +
  scale_y_continuous(limits = c(-50, 50), expand = expand_scale(0)) +
  scale_x_discrete() +
  annotate("rect", xmin = 0, xmax = 12.5, ymin = 0, ymax = 50, alpha = 0.1, fill = "skyblue2") +
  annotate("rect", xmin = 0, xmax = 12.5, ymin = -50, ymax = 0, alpha = 0.1, fill = "red2") +
  annotate("text", x = 11, y = 40, label = "Clonotypes\nGained", colour = "navy") +
  annotate("text", x = 11, y = -40, label = "Clonotypes\nLost", colour = "red2") +
  geom_hline(yintercept = 0) +
  geom_col(aes(fill = ClusterID), colour = "black", size = 0.25) +
  scale_fill_manual(values = c(intg_clusters_colors, "5E" = "brown", "5M" = "gold2"), guide = "none") +  
  ylab("Change in Clonotypes") + xlab("") +
  theme.lra() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())
net_sharing_plot
ggsave(net_sharing_plot, file = "plots/net_sharing_plot.pdf", width = 5, height = 4, device = cairo_pdf)

# Newly Spread Clonotypes (Figure 4J)
clonotype_spread_df <- data.frame(ClusterID = "", ClonotypeCount = 0, SpreadCount = 0, stringsAsFactors = F)
for (cluster_id in pre_cluster_ids) {
  total_cluster_clonotypes <- 0
  total_cluster_spread <- 0
  for (clonotype_df in treated_list) {
    pre_clonotypes <- subset(clonotype_df, !endsWith(ClusterID, "'"))
    pre_agg <- aggregate(ClusterID ~ Clonotype, pre_clonotypes, c)
    cluster_clonotypes <- subset(pre_agg, ClusterID %in% cluster_id)$Clonotype
    cluster_clonotypes <- unique(cluster_clonotypes[startsWith(cluster_clonotypes, "U_")])
    pre_post_ids <- c(cluster_id, paste0(cluster_id, "'"))
    cluster_spread_clonotypes <- unique(subset(clonotype_df, Clonotype %in% cluster_clonotypes & !(ClusterID %in% pre_post_ids))$Clonotype)
    total_cluster_clonotypes <- total_cluster_clonotypes + length(cluster_clonotypes)
    total_cluster_spread <- total_cluster_spread + length(cluster_spread_clonotypes)
    
    post_clonotypes <- subset(clonotype_df, endsWith(ClusterID, "'"))  
    post_agg <- aggregate(ClusterID ~ Clonotype, post_clonotypes, c)
    post_subset <- subset(post_agg, Clonotype %in% cluster_spread_clonotypes)
  }
  
  current_spread_df <- data.frame(ClusterID = paste(cluster_id, collapse = "-"), ClonotypeCount = total_cluster_clonotypes, SpreadCount = total_cluster_spread, stringsAsFactors = F)
  clonotype_spread_df <- rbind(clonotype_spread_df, current_spread_df)
}

clonotype_spread_df <- clonotype_spread_df[-c(1), ]
clonotype_spread_df$PercentSpread <- round(100 * (clonotype_spread_df$SpreadCount / clonotype_spread_df$ClonotypeCount), 1)
clonotype_spread_df[10:12, "PercentSpread"] <- 0
clonotype_spread_df$ClusterID <- factor(spread_df$ClusterID, levels = pre_cluster_ids)

clonotype_spread_plot <- ggplot(data = clonotype_spread_df, aes(x = ClusterID, y = PercentSpread)) +
  geom_col(aes(fill = ClusterID), colour = "black", size = 0.2) +
  scale_fill_manual(values = c(intg_clusters_colors, "5E" = "brown", "5M" = "gold2"), guide = "none") +  
  scale_y_continuous(expand = expand_scale(0), limits = c(0, 62)) + 
  xlab("") +
  ylab("Newly Shared Clonotypes (%)") +
  theme.lra() + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank())
spreadclonotype_spread_plot_plot
ggsave(clonotype_spread_plot, file = "spread_plot_extended_treated.pdf", width = 5, height = 4, device = cairo_pdf)

