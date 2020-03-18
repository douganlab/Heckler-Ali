# CDK4PlotHelpers.R
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

library(ggplot2)

intg_clusters_colors <- c("0" = "#CB449F", "1" = "palevioletred", "2" = "#F9989A", 
                          "3" = "#6CC6C3", "4" = "#68397A", "5" = "#F5CB1A" , 
                          "6" = "#B98E92", "7" = "#3D5EA9", "8" = "#CBC8BB",
                          "9" = "#DF0000", "10" = "#ADACCB", "11" = "black")


draw.cluster.umap <- function(seurat_object, save = FALSE) {
  umap_df <- as.data.frame(Reductions(seurat_object, slot = "umap")@cell.embeddings)
  umap_df$Cluster <- Idents(seurat_object)
  umap_cluster_intg <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(shape = 16, aes(colour = Cluster), alpha = 1, size = 0.75) +
    scale_colour_manual(values = intg_clusters_colors) +
    guides(colour = guide_legend(override.aes = list(size = 8))) + 
    xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + 
    theme.lra() + 
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.line = element_line(size = 0.5, arrow = arrow(length = unit(0.1, "inches"))))
  if (save) {
    ggsave("plots/umap_cluster_intg.pdf", umap_cluster_intg, width = 12, height = 9.5, device = cairo_pdf)
  } else {
    umap_cluster_intg
  }
}


draw.violin.sandwich <- function(seurat_object, genes, save = TRUE) {
  count_data <- GetAssayData(seurat_object$RNA)
  cluster_ids <- Idents(seurat_object)
  gene_expr <- count_data[genes,]
  plottable_expr <- reshape2::melt(t(as.matrix(gene_expr)))
  plottable_expr$Cluster <- rep(cluster_ids, length(genes))
  colnames(plottable_expr) <- c("CellID", "Gene", "Expression", "ClusterID")
  
  expr_violinp <- ggplot(data = plottable_expr, aes(x = ClusterID, y = Expression)) +
    geom_violin(aes(fill = ClusterID, colour = ClusterID), size = 0.8, alpha = 0.9, scale = "width") + 
    facet_grid(rows = vars(Gene)) + 
    scale_colour_manual(values = intg_clusters_colors) +
    scale_fill_manual(values = intg_clusters_colors) +
    xlab("") + ylab("") +
    theme.lra() + guides(fill = "none", colour = "none") +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "cm"),
          plot.title = element_blank(),
          strip.text = element_text(size = 8),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  if (save) {
    ggsave(paste0("plots/", gene, ".pdf"), plot = expr_violinp, width = 5, height = 5, device = cairo_pdf)
  } else {
    expr_violinp
  }
}


draw.monotone.umap <- function(seurat_object, colour = "black", save = FALSE) {
  umap_df <- as.data.frame(Reductions(seurat_object, slot = "umap")@cell.embeddings)
  umap_plot <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(alpha = 0.9, shape = 19, colour = colour, size = 0.75) +
    guides(colour = guide_legend(override.aes = list(size = 4))) + 
    xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + 
    theme.lra() + 
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 0.5, arrow = arrow(length = unit(0.1, "inches"))))
  if (save) {
    ggsave(paste0("plots/mono_umap.pdf"), umap_plot, width = 6, height = 5.5, device = cairo_pdf)
  } else {
    umap_plot
  }
}


draw.clonal.umap <- function(seurat_object, save = F) {
  umap_df <- as.data.frame(Reductions(seurat_object, slot = "umap")@cell.embeddings)
  umap_df$TCRFamily <- seurat_object$tcr_family
  umap_plot <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(shape = 16, aes(colour = TCRFamily), size = 0.75, alpha = 0.9) +
    scale_colour_manual(values = c("hotpink", "hotpink", "grey85", "grey85"), guide = "none") +
    xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + 
    theme.lra() + 
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 0.5, arrow = arrow(length = unit(0.1, "inches"))))
  if (save) {
    ggsave(paste0("plots/clonal_umap.pdf"), umap_plot, width = 6, height = 5.5, device = cairo_pdf)
  } else {
    umap_plot
  }
}


plot.sharing.mat <- function(sharing_mat, prefix = "", high_color = "red3", save = TRUE) {
  sharing_df <- reshape2::melt(sharing_mat)
  colnames(sharing_df) <- c("ClusterIDY", "ClusterIDX", "Count")
  grid_plot <- ggplot(data = sharing_df, aes(x = ClusterIDX, y = ClusterIDY, fill = Count)) +
    geom_tile(colour = "grey60") +
    geom_text(aes(label = Count), colour = ifelse(abs(sharing_df$Count) > 0, "black", "grey75"), family = "GillSans") + 
    scale_x_discrete(expand = expand_scale(add = 0)) + 
    scale_y_discrete(expand = expand_scale(add = 0)) + 
    geom_vline(xintercept = 13.5) +
    geom_hline(yintercept = 13.5) +
    geom_vline(xintercept = 0.5) +
    geom_hline(yintercept = 0.5) +
    geom_vline(xintercept = 26.5) +
    geom_hline(yintercept = 26.5) +
    xlab("") + ylab("") + theme.lra() + 
    theme(panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank())
  
  if (sum(sharing_mat, na.rm = T) > 0) {
    grid_plot <- grid_plot + scale_fill_gradient(limits = c(0, max(sharing_mat, na.rm = T) - 1), oob = scales::censor, low = "white", high = high_color, na.value = "grey90", guide = "none")
  } else {
    grid_plot <- grid_plot + scale_fill_gradient(limits = c(min(sharing_mat, na.rm = T) / 2, 0), low = "royalblue3", high = "white", na.value = "grey90", guide = "none")
  }
  
  if (save) {
    ggsave(paste0("plots/", prefix, "clonotypes_grid.pdf"), grid_plot, width = 6, height = 6, device = cairo_pdf)
  } else {
    grid_plot
  }
}

theme.lra <- function() {
  return(theme_bw() + theme(text = element_text(family = "Gill Sans", size = 14),
                            plot.title = element_text(family = "Gill Sans", face = "bold", size = 12),
                            legend.text = element_text(size = 12),
                            legend.title = element_text(size = 12),
                            axis.line = element_line(size = 0.35),
                            axis.title.x = element_text(size = 16, margin = margin(t = 10)),
                            axis.title.y = element_text(size = 16, margin = margin(r = 10)),
                            axis.text.x = element_text(size = 13, margin = margin(t = 5)),
                            axis.text.y = element_text(size = 13)))
}



