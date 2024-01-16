#### Small simulated dataset to compare TreePPL and RevBayes ####

#install_github("treeppl/treepplr")

library(tidyverse)
library(ape)
library(ggtree)
library(evolnets)
library(treepplr)

tree <- rcoal(4, rooted = TRUE)
tree$tip.label <- paste0("S", 1:4)
plot(tree)
axisPhylo()

height <- node.depth.edgelength(tree)[1]
scaling_factor <- 2.0/height
tree$edge.length <- tree$edge.length*scaling_factor
plot(tree)
axisPhylo()

is.binary(tree)
is.ultrametric(tree)
is.rooted(tree)

host_tree <- rcoal(3, rooted = TRUE)
host_tree$tip.label <- paste0("H", 1:3)
plot(host_tree)
axisPhylo()

matrix <- matrix(data = c(2,0,2, 2,0,0, 2,2,0, 0,2,0), nrow = 4, ncol = 3, byrow = TRUE)
rownames(matrix) <- tree$tip.label
colnames(matrix) <- host_tree$tip.label

# write.csv(matrix, "data/matrix.csv", row.names = TRUE)
# write.nexus.data(matrix, "data/matrix.nex", format = "standard")
# write.tree(host_tree, "data/host_tree.tre")
#
# write.tree(tree, "data/tree.tre")
# Add subroot branch to tree
tree_string <- readLines("data/tree.tre")
tree_tiny_stem_string <- sub(");$", "):0.01;", tree_string)
# writeLines(tree_tiny_stem_string, "data/tree_tiny_stem.tre")
tree_long_stem_string <- sub(");$", "):2.0;", tree_string)
# writeLines(tree_long_stem_string, "data/tree_long_stem.tre")



###### Prepare input for treepplr ######
#
# symbiont_tree: TreeLabeled, ntips: Int, nhosts: Int,
# interactions: Int[], host_distances: Real[],
# dMean: Real, tune: Real

symbiont_tree <- read_tree_from_revbayes("data/tree_tiny_stem_Rev.tre")
plot_data <- plot_matrix_phylo(matrix, at_nodes = NULL, symbiont_tree, host_tree, find_modules = FALSE)
plot_data[[1]] <- plot_data[[1]] + geom_nodelab(size = 3, hjust = -0.1)
plot_data

ntips <- Ntip(symbiont_tree)
nhosts <- Ntip(host_tree)
interactions <- matrix
host_distances <- cophenetic.phylo(host_tree)
dMean <- sum(host_distances) / factorial(nhosts)
tune <- 0.9

###



