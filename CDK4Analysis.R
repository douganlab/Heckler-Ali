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
cd45dp_full_int <- RunPCA(cd45dp_full_int)=
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

# Figure 3C
draw.cluster.umap(cd45dp_no11_int)

# Figure 4F
draw.clonal.umap(cd45dp_no11_int)

# Supplementary Figure X
draw.monotone.umap(cd45dp_no11_int[, cd45dp_no11_int$treatment == "Healthy"], "grey20", T)
draw.monotone.umap(cd45dp_no11_int[, cd45dp_no11_int$treatment == "Cancer-NoCDK4/6i"], "#009E73", T)
draw.monotone.umap(cd45dp_no11_int[, cd45dp_no11_int$treatment == "Cancer-PostCDK4/6i"], "#E69F00", T)

# Figure 3F
markers <- c("TCF7", "JUN", "IL7R", "HLA-DRB1","GZMH", "AQP3", "DUSP2", "CCL5", "GATA3", "NKG7", "GNLY", "S100A9", "ZNF683", "IRF7")
draw.violin.sandwich(cd45dp_no11_int, markers, T)

# Patient Contribution to Clusters (Figure 3D)
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

# Figure 3H
trm_genes <- c("ZNF683", "CD7", "DUSP2",
               "KLRC3", "KLRC2", "NCR3",
               "FXYD2", "LEF1", "STK17A",
               "FXYD7", "IL7R", "LINC02446",
               "CD300A", "JUNB", "IFITM3", 
               "TCF7", "GNLY", "ITGAE", "CD69")

trm_exp <- AverageExpression(cd45dp_no11_int, assay = "RNA", features = trm_genes)
trm_exp_mat <- as.matrix(trm_exp$RNA[trm_genes,])
trm_exp_z_mat <- t(scale(t(trm_exp_mat)))

trm_exp_df <- reshape2::melt(trm_exp_z_mat)
colnames(trm_exp_df) <- c("Gene", "ClusterID", "Expression")
trm_exp_df$Gene <- factor(trm_exp_df$Gene, levels = rev(trm_genes))
trm_exp_df$ClusterID <- factor(trm_exp_df$ClusterID, levels =  as.character(0:10))

trm_heatmap <- ggplot(data = trm_exp_df, aes(x = ClusterID, y = Gene, fill = Expression)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "royalblue3", mid = "white", high = "red3", limits = c(-2, 3)) +
  scale_x_discrete(expand = expand_scale(0), position = "top", labels = 0:10) +
  scale_y_discrete(expand = expand_scale(0)) +
  theme.lra() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 12),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom")
trm_heatmap


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

# Retrospectively label main Seurat object
cd45dp_no11_int <- SetIdent(cd45dp_no11_int, cells = WhichCells(cd45dp_sub_int, idents = "5E"), value = "5E")
cd45dp_no11_int <- SetIdent(cd45dp_no11_int, cells = WhichCells(cd45dp_sub_int, idents = "5M"), value = "5M")
draw.highlighted.umap(cd45dp_no11_int, highlight_cells = WhichCells(cd45dp_sub_int, idents = "5E"), highlight_color = "brown", T)
draw.highlighted.umap(cd45dp_no11_int, highlight_cells = WhichCells(cd45dp_sub_int, idents = "5M"), highlight_color = "gold2", T)

# Figure 4D
fate_markers <- c("NKG7", "GZMH", "GNLY", 
                  "LGALS1", "ITGB1", "ITGB2",
                  "ZEB2", "HLA-DRB1", "CD74",
                  "KLRC3", "MT1F", "TIGIT",
                  "GZMK", "DUSP2", "LTB", 
                  "NR4A2", "FOS", "ZFP36L2", 
                  "IL7R", "AQP3")

fate_exp <- AverageExpression(cd45dp_no11_int, assay = "RNA", features = fate_markers)
fate_exp_mat <- as.matrix(fate_exp$RNA[, c("4", "5M", "5E", "3", "7")])
fate_exp_mat <- cbind(fate_exp_mat, rowMeans(fate_exp_mat[, c("3", "7")]))
fate_exp_mat <- fate_exp_mat[, -c(4,5)]
colnames(fate_exp_mat) <- c("Memory", "5M", "5E", "Effector")
fate_exp_z_mat <- t(scale(t(fate_exp_mat)))
fate_exp_df <- reshape2::melt(fate_exp_z_mat)
colnames(fate_exp_df) <- c("Gene", "ClusterID", "Expression")
fate_exp_df$ClusterID <- factor(fate_exp_df$ClusterID, levels = rev(c("Memory", "5M", "5E", "Effector")))

fate_heatmap <- ggplot(data = fate_exp_df, aes(x = Gene, y = ClusterID, fill = Expression)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "royalblue3", mid = "white", high = "red3", limits = c(-1.5, 1.5), breaks = c(-1, 0, 1), name = "Expression\n(z-score)") +
  scale_y_discrete(expand = expand_scale(0), position = "left") +
  scale_x_discrete(expand = expand_scale(0)) +
  theme.lra() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, hjust = 1, angle = 90),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")
fate_heatmap
ggsave(fate_heatmap, file = "plots/fate_heamtap.pdf", width = 6, height = 4, device = cairo_pdf)

#### MYC Score #### 

# List was generated with ARACNE
myc_genes <- read.csv("myc_target_genes.csv", stringsAsFactors = F, header = F)$V1
myc_gene_expr <- as.matrix(GetAssayData(cd45dp_no11_int$RNA)[myc_genes,])
myc_gene_expr_avg <- rowMeans(myc_gene_expr)
myc_gene_expr_sd <- rowSds(myc_gene_expr)
myc_gene_expr_z <- (myc_gene_expr - myc_gene_expr_avg) / myc_gene_expr_sd
myc_score <- colSums(myc_gene_expr_z)

sample_id_by_cell <- sapply(Cells(cd45dp_no11_int), function(cellname) {
  all_subject_ids[sapply(all_subject_ids, grepl, cellname)]
})

myc_score_df <- data.frame("SampleID" = sample_id_by_cell, "ClusterID" = Idents(cd45dp_no11_int),  "MycScore" = myc_score)
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

# Figure 4C
draw.gradient.umap(cd45dp_sub_int_pre, myc_score_cl5_pre_df$MycScore, save = T)
draw.gradient.umap(cd45dp_sub_int_post, myc_score_cl5_post_df$MycScore, save = T)

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

calculate.correction.matrix <- function() {
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
  return(matrix(correction_factors, nrow = length(pre_cluster_ids), ncol = length(pre_cluster_ids), byrow = F))
}

mini_correction_mat <- calculate.correction.matrix()

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
plot.sharing.mat(corrected_pre_sharing_mat, prefix = "pre_norm_", high_color = "green4", limit_scale = 60, save = T)

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

plot.sharing.mat(total_added_sharing_mat, prefix = "added_", save = T)

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

plot.sharing.mat(total_removed_sharing_mat, prefix = "removed_", save = T)


