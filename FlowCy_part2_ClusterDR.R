rm(list = ls())

# Load packages
library(rstudioapi)
library(devtools)
library("flowCore")
library("flowWorkspace")
library(cytofCore)
library(FlowSOM)
# library(cytofkit)
library(cluster)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(flowViz)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(uwot)
library(CATALYST)
library(Rphenograph)
library(diffcyt)
library(SummarizedExperiment)
library(stringr)
library(ggcyto)
library(SingleCellExperiment)
library(scran)
library(scater)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory
# Define workingDirectory
wdName <- "Working_DirectoryFCS"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")

setwd(workingDirectory)

# Load workspace and SCEobject
load("workspaceSCE.rds")

sce
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "mean")
# Run FlowSOM and ConsensusClusterPlus
seed <- 123456
set.seed(seed)
sce <- cluster(sce, features = "type", xdim = 10, ydim = 10, maxK = 20, 
               verbose = TRUE, seed = seed)
delta_area(sce)
# Run dimensionality reduction
n_cells <- 1000
n_cells < n_events
exaggeration_factor <- 12.0
eta <- n_cells/exaggeration_factor
sce <- runDR(sce, dr = "TSNE", cells = n_cells, features = "type", theta = 0.5, max_iter = 1000, 
             distMethod = "euclidean",
             PCA = TRUE, eta = eta, exaggeration_factor = 12.0)
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")
sce <- runDR(sce, dr = "DiffusionMap", cells = n_cells, features = "type", assay = "exprs")

save(list = ls(), file = "workspaceSCEclusterDR.rds")

# Plots
display.brewer.all(colorblindFriendly = TRUE)
delta_area(sce)
cdx <- type_markers(sce)
plotMultiHeatmap(sce, k = "meta10",
                 hm1 = cdx, hm2 = "abundances", 
                 bars = TRUE, perc = TRUE, row_anno = FALSE)
plotMultiHeatmap(sce, k = "meta10",
                 hm1 = cdx, hm2 = "abundances", 
                 bars = TRUE, perc = TRUE, row_anno = FALSE)

plotClusterHeatmap(sce, hm2 = NULL, k = "meta8", m = NULL, cluster_anno = TRUE, draw_freqs = TRUE,
                   draw_dend = TRUE)
plotClusterHeatmap(sce, hm2 = NULL, k = "meta12", m = NULL, cluster_anno = TRUE, draw_freqs = TRUE,
                   draw_dend = TRUE, scale = "last")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta8", by = "cluster_id",  fun = "mean", scale = "last")

plotDR(sce, dr = "UMAP", color_by = "condition") + scale_color_manual(values = c("blue", "red"))
plotDR(sce, dr = "UMAP", color_by = "condition")
# UMAP plot color_by = "meta8", facet_by = "condition"
plotDR(sce, dr = "DiffusionMap", color_by = "meta10", facet_by = "condition")

plotDR(sce, dr = "UMAP", color_by = "condition", facet_by = "condition")
plotDR(sce, dr = "UMAP", color_by = "meta8", facet_by = "condition")
plotDR(sce, dr = "TSNE", color_by = "condition", facet_by = "condition")
plotDR(sce, dr = "DiffusionMap", color_by = "condition") + scale_color_manual(values = c("blue", "red"))
markers <- c("CD27", "CD69", "PD1", "CD57")
plotDR(sce, dr = "DiffusionMap", color_by = markers)
plotDR(sce, dr = "UMAP", color_by = markers)
plotDR(sce, dr = "TSNE", color_by = markers)
# UMAP plot color_by = "meta8", facet_by = "sample_id"
plot <- plotDR(sce, dr = "UMAP", color_by = "meta8", facet_by = "sample_id")
plot$facet$params$ncol <- 3
plot

plot <- plotDR(sce, dr = "DiffusionMap", color_by = "meta8", facet_by = "sample_id")
plot$facet$params$ncol <- 3
plot
# UMAP color_by CD27 DNAM1
plotDR(sce, dr = "UMAP", color_by = c("CD27", "DNAM1"), facet_by = "condition")

plotDR(sce, dr = "UMAP", color_by = "meta8", facet_by = "sample_id")
plotAbundances(sce, k = "meta12", by = "cluster_id", group_by = "condition")
plotAbundances(sce, k = "meta12", by = "sample_id", group_by = "condition")
