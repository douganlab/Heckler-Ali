#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3

# CDK4RNAVelocity.py
# Copyright (C) 2019 Lestat Ali
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

import scvelo as scv

# Loom files were generated from each sample's BAM using Velocyto.
# They were then concatenated into one large loom file using loompy.combine().
loom_data = scv.read("Combined.loom", cache = True)

# Clean-up
noncluster11_cells = pandas.read_csv("seurat_cells.csv", index_col = 0)
noncluster11_cells = noncluster11_cells.to_numpy().flatten().tolist()
loom_data = loom_data[noncluster11_cells,:]

scv.pp.filter_and_normalize(loom_data, min_shared_counts = 30, n_top_genes = 2000)
scv.pp.moments(loom_data, n_pcs = 50, n_neighbors = 50)

scv.tl.velocity(loom_data)
scv.tl.velocity_graph(loom_data)
scv.tl.rank_velocity_genes(loom_data, resolution = 0.4)
scv.tl.velocity_clusters(loom_data)
scv.tl.velocity_pseudotime(loom_data)

# Import UMAP from Seurat
umap_data = pandas.read_csv("seurat_umap_emb.csv", index_col = 0)
loom_data.obsm["X_umap"] = umap_data.to_numpy()

# Import clusters from Seurat
cluster_data = pandas.read_csv("seurat_clusters.csv", index_col = 0)
loom_data.obs["seurat_cluster"] = cluster_data
cluster_colors = pandas.read_csv("seurat_cluster_colors.csv", index_col = 0)
loom_data.uns["seurat_cluster_colors"] = cluster_colors.to_numpy().flatten().tolist()

# Visualize
scv.pl.velocity_embedding_stream(loom_data, basis = "umap", save = "stream.png", color = "seurat_cluster", alpha = 0.5, title = "", legend_fontsize = 0, dpi = 300)
scv.pl.velocity_embedding_grid(loom_data, basis = "umap", save = "grid.png", color = "seurat_cluster", arrow_length = 4, dpi = 300, arrow_color = "black", alpha = 0.5, title = "")
