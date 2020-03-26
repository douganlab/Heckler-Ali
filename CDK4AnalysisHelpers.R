# CDK4AnalysisHelpers.R
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


create.seurat.object.for.patient <- function(patient_id) {
  print("Loading pre-treatment data...")
  raw_data_pre <- Read10X(data.dir = paste0("10x/", patient_id, "-Pre/raw_feature_bc_matrix"))
  print("Loading post-treatment data...")
  raw_data_post <- Read10X(data.dir = paste0("10x/", patient_id, "-Post/raw_feature_bc_matrix"))
  
  print("Merging...")
  pre_seurat <- create.unnormalized.seurat.object(raw_data_pre, patient_id)
  post_seurat <- create.unnormalized.seurat.object(raw_data_post, patient_id)
  
  merged_seurat <- merge(pre_seurat, y = post_seurat, add.cell.ids = c("Pre", "Post"))
  
  merged_seurat$orig.ident <- rep(paste0(patient_id, "-Pre"), ncol(merged_seurat))
  merged_seurat$orig.ident[startsWith(Cells(merged_seurat), "Post")] <- paste0(patient_id, "-Post")
  
  merged_seurat <- NormalizeData(merged_seurat)
  merged_seurat <- FindVariableFeatures(merged_seurat, nfeatures = 2500)
  merged_seurat <- ScaleData(merged_seurat, features = rownames(merged_seurat))
  merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))
  return(merged_seurat)
}


create.unnormalized.seurat.object <- function(raw_data, patient_id) {
  seurat_object <- CreateSeuratObject(counts = raw_data, project = patient_id, min.cells = 3, min.features = 100)
  seurat_object[["percent_mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  feature_cutoff_top <- mean(seurat_object$nFeature_RNA) + 2 * sd(seurat_object$nFeature_RNA)
  feature_cutoff_btm <- mean(seurat_object$nFeature_RNA) - 2 * sd(seurat_object$nFeature_RNA)
  mt_cutoff <- median(seurat_object$percent_mt) + 2 * mad(seurat_object$percent_mt)
  filter <- (seurat_object[["percent_mt"]] < mt_cutoff
             & seurat_object[["nFeature_RNA"]] < feature_cutoff_top 
             & seurat_object[["nFeature_RNA"]] > feature_cutoff_btm)
  seurat_object <- seurat_object[, filter]
  if ("IGHM" %in% rownames(seurat_object)) {
    seurat_object <- subset(seurat_object, subset = (IGHM == 0))
  } 
  if ("IGHD" %in% rownames(seurat_object)) {
    seurat_object <- subset(seurat_object, subset = (IGHD == 0))
  } 
  if ("CSF3R" %in% rownames(seurat_object)) {
    seurat_object <- subset(seurat_object, subset = (CSF3R == 0))
  } 
  return(seurat_object)
}


create.seurat.object.for.healthy.subject <- function(subject_id) {
  raw_data <- Read10X(paste0(data.dir = "10x/", subject_id, "/raw_feature_bc_matrix"))
  seurat_object <- CreateSeuratObject(counts = raw_data, project = subject_id, min.cells = 3, min.features = 200)
  seurat_object[["percent_mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  feature_cutoff_top <- mean(seurat_object$nFeature_RNA) + 3 * sd(seurat_object$nFeature_RNA)
  feature_cutoff_btm <- mean(seurat_object$nFeature_RNA) - 3 * sd(seurat_object$nFeature_RNA)
  mt_cutoff <- mean(seurat_object$percent_mt) + 2 * sd(seurat_object$percent_mt)
  filter <- (seurat_object[["percent_mt"]] < mt_cutoff
             & seurat_object[["nFeature_RNA"]] < feature_cutoff_top 
             & seurat_object[["nFeature_RNA"]] > feature_cutoff_btm)
  seurat_object <- seurat_object[, filter]
  if ("IGHM" %in% rownames(seurat_object)) {
    seurat_object <- subset(seurat_object, subset = (IGHM == 0))
  } 
  if ("IGHD" %in% rownames(seurat_object)) {
    seurat_object <- subset(seurat_object, subset = (IGHD == 0))
  } 
  if ("CSF3R" %in% rownames(seurat_object)) {
    seurat_object <- subset(seurat_object, subset = (CSF3R == 0))
  } 
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, nfeatures = 2500)
  seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  return(seurat_object)
}


annotate.with.tcr.data <- function(merged_seurat_object, patient_id) {
  if (startsWith(patient_id, "S")) {
    tcr_annot_all <- get.tcr.family.annotations(patient_id, merged_seurat_object)
  } else {
    pre_id <- paste0(patient_id, "-Pre")
    seurat_pre <- merged_seurat_object[, merged_seurat_object$orig.ident == pre_id]
    seurat_pre <- RenameCells(seurat_pre, new.names = gsub("Pre_", "", Cells(seurat_pre)))
    tcr_annot_pre <- get.tcr.family.annotations(pre_id, seurat_pre)
    
    post_id <- paste0(patient_id, "-Post")
    seurat_post <- merged_seurat_object[, merged_seurat_object$orig.ident == post_id]
    seurat_post <- RenameCells(seurat_post, new.names = gsub("Post_", "", Cells(seurat_post)))
    tcr_annot_post <- get.tcr.family.annotations(post_id, seurat_post)
    
    tcr_annot_all <- c(tcr_annot_pre, tcr_annot_post)
  }
  
  merged_seurat_object <- AddMetaData(merged_seurat_object, tcr_annot_all, "tcr_family")
  return(merged_seurat_object)
}


get.tcr.family.annotations <- function(sample_id, seurat_object) {
  tcr_family <- rep("x", ncol(seurat_object))
  names(tcr_family) <- colnames(seurat_object)
  
  all_tcr_df <- read.csv(paste0("10x/VDJ-", sample_id, "/all_contig_annotations.csv"), stringsAsFactors = F)
  
  # Identify MAITs
  good_v_barcodes <- subset(all_tcr_df, v_gene == "TRAV1-2")$barcode
  good_j_barcodes <- subset(all_tcr_df, j_gene == "TRAJ33" | j_gene == "TRAJ12" | j_gene == "TRAJ20")$barcode
  mait_barcodes <- intersect(good_v_barcodes, good_j_barcodes)
  mait_barcodes <- sub("-1", "", mait_barcodes)
  valid_mait_barcodes <- intersect(mait_barcodes, colnames(seurat_object))
  tcr_family[valid_mait_barcodes] <- "mait"
  
  # Identity clonally expanded cells
  tcr_with_clonotype_df <- subset(all_tcr_df, raw_clonotype_id != "None" & is_cell == "True")
  tcr_barcodes <- unique(tcr_with_clonotype_df$barcode)
  tcr_clonotypes <-  unique(tcr_with_clonotype_df$raw_clonotype_id)
  cloned_barcodes <- c()
  for (current_barcode in tcr_barcodes) {
    matching_tcr_subset <- subset(tcr_with_clonotype_df, barcode == current_barcode)
    clonotypes <- matching_tcr_subset$raw_clonotype_id
    unique_clonotypes <- unique(clonotypes)
    if (length(unique_clonotypes) > 1) {
      print(paste("WARNING: this cell", current_barcode, "has more than 1 matching clonotype"))
    }
    clones <- subset(tcr_with_clonotype_df, raw_clonotype_id %in% unique_clonotypes & barcode != current_barcode)
    if (nrow(clones) > 0) {
      cloned_barcodes <- c(cloned_barcodes, current_barcode)
    }
  }
  cloned_barcodes <- sub("-1", "", cloned_barcodes)
  cloned_nonmait_barcodes <- setdiff(cloned_barcodes, mait_barcodes)
  valid_cloned_nonmait_barcodes <- intersect(cloned_nonmait_barcodes, colnames(seurat_object))
  tcr_family[valid_cloned_nonmait_barcodes] <- "clonal"
  
  return(tcr_family)
}


get.unified.clonotyped.barcodes <- function(pre_sample_id, post_sample_id, seurat_object) {
  pre_clonotypes <- get.clonotyped.barcodes(pre_sample_id, seurat_object)
  post_clonotypes <- get.clonotyped.barcodes(post_sample_id, seurat_object)
  
  pre_clonotype_metadata <- read.csv(paste0("10x/VDJ-", pre_sample_id, "/clonotypes.csv"), stringsAsFactors = F)
  post_clonotype_metadata <- read.csv(paste0("10x/VDJ-", post_sample_id, "/clonotypes.csv"), stringsAsFactors = F)
  
  pre_resolution_map <- rep("", nrow(pre_clonotype_metadata))
  names(pre_resolution_map) <- pre_clonotype_metadata$clonotype_id
  post_resolution_map <-  rep("", nrow(post_clonotype_metadata))
  names(post_resolution_map) <- post_clonotype_metadata$clonotype_id
  
  for (i in 1:nrow(pre_clonotype_metadata)) {
    current_cd3_nt <- pre_clonotype_metadata$cdr3s_nt[i]
    current_pre_clonotype_id <- pre_clonotype_metadata$clonotype_id[i]
    if (any(post_clonotype_metadata$cdr3s_nt == current_cd3_nt)) {
      matching_post_clonotype_id <- post_clonotype_metadata$clonotype_id[post_clonotype_metadata$cdr3s_nt == current_cd3_nt]
      pre_resolution_map[i]  <- paste0("U_", current_pre_clonotype_id, matching_post_clonotype_id)
    } else {
      pre_resolution_map[i] <- paste0("A_", current_pre_clonotype_id)
    }
  }
  
  for (j in 1:nrow(post_clonotype_metadata)) {
    current_cd3_nt <- post_clonotype_metadata$cdr3s_nt[j]
    current_post_clonotype_id <- post_clonotype_metadata$clonotype_id[j]
    if (any(pre_clonotype_metadata$cdr3s_nt == current_cd3_nt)) {
      matching_pre_clonotype_id <- pre_clonotype_metadata$clonotype_id[pre_clonotype_metadata$cdr3s_nt == current_cd3_nt]
      post_resolution_map[j]  <- paste0("U_", matching_pre_clonotype_id, current_post_clonotype_id)
    } else {
      post_resolution_map[j] <- paste0("B_", current_post_clonotype_id)
    }
  }
  
  pre_clonotypes_resolved <- pre_clonotypes
  pre_clonotypes_resolved$Clonotype <- pre_resolution_map[pre_clonotypes$Clonotype]
  
  post_clonotypes_resolved <- post_clonotypes
  post_clonotypes_resolved$Clonotype <- post_resolution_map[post_clonotypes$Clonotype]
  post_clonotypes_resolved$ClusterID <- paste0(post_clonotypes_resolved$ClusterID, "'")
  
  all_clonotypes_resolved <- rbind(pre_clonotypes_resolved, post_clonotypes_resolved)
  
  return(all_clonotypes_resolved)
}


# EXCLUDES MAITs
get.clonotyped.barcodes <- function(sample_id, seurat_object) {
  print(paste0("Processing ", sample_id, "..."))
  if (startsWith(sample_id, "P")) {
    filtered_seurat <- seurat_object[, seurat_object$orig.ident == sample_id]
  } else {
    filtered_seurat <- seurat_object
  }
  all_tcr_df <- read.csv(paste0("10x/VDJ-", sample_id, "/all_contig_annotations.csv"), stringsAsFactors = F)
  # These conditions are rather redundant but it doesn't hurt to be thorough
  filtered_tcr_df <- subset(all_tcr_df, is_cell == "True" & raw_clonotype_id != "None" & high_confidence == "True")
  filtered_tcr_df <- subset(filtered_tcr_df, !(v_gene == "TRAV1-2" & j_gene %in% c("TRAJ33", "TRAJ12", "TRAJ20")))
  filtered_tcr_df$barcode <- paste0(gsub("-", "_", sample_id), "_", filtered_tcr_df$barcode)
  filtered_tcr_df$barcode <- gsub("-1", "", filtered_tcr_df$barcode)
  unique_barcodes <- unique(filtered_tcr_df$barcode)
  unique_valid_barcodes <- unique_barcodes[unique_barcodes %in% Cells(filtered_seurat)]
  barcode_to_clonotype_df <- data.frame()
  for (barcode in unique_valid_barcodes) {
    matching_clonotypes <- filtered_tcr_df[filtered_tcr_df$barcode == barcode, "raw_clonotype_id"]
    unique_matching_clonotype <- unique(matching_clonotypes)
    if (length(unique_matching_clonotype) != 1) {
      print(paste("WARNING!", barcode))
    }
    int_cluster_id <- as.character(Idents(cd45dp_no11_int)[barcode])
    barcode_to_clonotype_df <- rbind(barcode_to_clonotype_df, c(barcode, unique_matching_clonotype, int_cluster_id), stringsAsFactors = F)
  }
  colnames(barcode_to_clonotype_df) <- c("Barcode", "Clonotype", "ClusterID")
  
  duplicated(barcode_to_clonotype_df$Clonotype)
  expanded_clonotype_df <- subset(barcode_to_clonotype_df, duplicated(Clonotype))
  print(paste0("Found ", nrow(expanded_clonotype_df), " clonally expanded T cells. Done!"))
  return(expanded_clonotype_df)
}


enumerate.shared.clonotypes <- function(patient_clonotypes, post.treatment = FALSE) {
  origin_cluster_ids <- c("0", "1", "2", "3", "4", "5E", "5M", "6", "7", "8", "9", "10")
  if (post.treatment) {
    target_cluster_ids <- paste0(origin_cluster_ids, "'")
  } else {
    target_cluster_ids <- origin_cluster_ids
  }
  shared_clonotype_df <- as.data.frame(matrix(data = 0,
                                              nrow = length(origin_cluster_ids), 
                                              ncol = length(target_cluster_ids)))
  rownames(shared_clonotype_df) <- origin_cluster_ids
  colnames(shared_clonotype_df) <- target_cluster_ids
  
  for (i in 1:nrow(shared_clonotype_df)) {
    for (j in 1:ncol(shared_clonotype_df)) {
      all_shared_clonotypes <- c()
      i_clonotypes <- subset(patient_clonotypes, ClusterID == origin_cluster_ids[i])$Clonotype
      j_clonotypes <- subset(patient_clonotypes, ClusterID == target_cluster_ids[j])$Clonotype
      all_shared_clonotypes <- c(all_shared_clonotypes, intersect(i_clonotypes, j_clonotypes))
      shared_clonotype_df[i, j] <- paste(all_shared_clonotypes, collapse = ",")
    }
  }
  
  return(shared_clonotype_df)
}


# method = {"new", "lost", "net"}
identify.delta.clonotypes <- function(pre_sharing_df, post_sharing_df, patient_clonotypes, method = "new") {
  delta_df <- post_sharing_df
  for (i in 1:nrow(delta_df)) {
    for (j in 1:ncol(delta_df)) {
      if (i == j) {
        delta_df[i,j] <- NA
      } else {
        pre_clonotypes <- unlist(strsplit(pre_sharing_df[i,j], ","))
        pre_clonotypes <- pre_clonotypes[startsWith(pre_clonotypes, "U_")]
        post_clonotypes <- unlist(strsplit(post_sharing_df[i,j], ","))
        post_clonotypes <- post_clonotypes[startsWith(post_clonotypes, "U_")]
        if (method == "new") {
          delta_df[i,j]  <- length(setdiff(post_clonotypes, pre_clonotypes))
        } else if (method == "lost") {
          delta_df[i,j]  <- -1 * length(setdiff(pre_clonotypes, post_clonotypes))
        } else {
          delta_df[i,j] <- length(post_clonotypes) - length(pre_clonotypes)
        }
      }
    }
  }
  return(data.matrix(delta_df))
}


compute.pre.clonotypes <- function(patient_clonotypes) {
  pre_sharing_df <- enumerate.shared.clonotypes(patient_clonotypes)
  pre_sharing_counts_df <- pre_sharing_df
  for (i in 1:nrow(pre_sharing_counts_df)) {
    for (j in 1:ncol(pre_sharing_counts_df)) {
      if (i == j) {
        pre_sharing_counts_df[i,j] <- NA
      } else {
        pre_clonotypes <- unlist(strsplit(pre_sharing_df[i,j], ","))
        pre_clonotypes <- pre_clonotypes[startsWith(pre_clonotypes, "U_")]
        pre_sharing_counts_df[i,j] <- length(pre_clonotypes)
      }
    }
  }
  return(data.matrix(pre_sharing_counts_df))
}


compute.added.clonotypes <- function(patient_clonotypes) {
  pre_sharing_df <- enumerate.shared.clonotypes(patient_clonotypes)
  post_sharing_df <- enumerate.shared.clonotypes(patient_clonotypes, T)
  added_sharing_df <- identify.delta.clonotypes(pre_sharing_df, post_sharing_df, patient_clonotypes, "new")
  return(added_sharing_df)
}


compute.removed.clonotypes <- function(patient_clonotypes) {
  pre_sharing_df <- enumerate.shared.clonotypes(patient_clonotypes)
  post_sharing_df <- enumerate.shared.clonotypes(patient_clonotypes, T)
  removed_sharing_df <- identify.delta.clonotypes(pre_sharing_df, post_sharing_df, patient_clonotypes, "lost")
  return(removed_sharing_df)
}

compute.net.change.clonotypes <- function(patient_clonotypes) {
  pre_cluster_ids <- c("0", "1", "2", "3", "4", "5E", "5M", "6", "7", "8", "9", "10")
  post_cluster_ids <- paste0(origin_cluster_ids, "'")
  
  pre_clonotype_count <- unlist(lapply(pre_cluster_ids, function(cluster_id) {
    unique_cluster_clonotypes <- unique(subset(patient_clonotypes, ClusterID == cluster_id)$Clonotype)
    length(which(startsWith(unique_cluster_clonotypes, "U_")))
  }))
  
  post_clonotype_count <- unlist(lapply(post_cluster_ids, function(cluster_id) {
    unique_cluster_clonotypes <- unique(subset(patient_clonotypes, ClusterID == cluster_id)$Clonotype)
    length(which(startsWith(unique_cluster_clonotypes, "U_")))
  }))
  
  delta_clonotype_count <- post_clonotype_count - pre_clonotype_count
  names(delta_clonotype_count) <- pre_cluster_ids
  return(delta_clonotype_count)
}

