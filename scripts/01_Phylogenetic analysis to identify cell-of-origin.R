## 01_Phylogenetic analysis to identify cell-of-origin.R
## 
## Purpose: To identify cell-of-origing (COO) on B cell developmental axis
## 
## 
##
##
########################################################################################################################



#########################################################
# libraries
#
library(ape)
library(openxlsx)
library(phylobase)
library(gdata)




###########################################################
# directories
#
analysis.dir <-  "/Users/justynaanna/Documents/CLL/Github/MethylCOOM/analysis/" # output directory for the analysis
data.dir <-  "/Users/justynaanna/Documents/CLL/Github/MethylCOOM/data/" # directory containing linear B cell-specific CpGs






############################################################
# load methylation data for linear B-cell specific CpGs and normal B cells
meth.data.CLL <- read.xlsx(file.path(data.dir, "linear.B.cell.specificCpGs.Oakes.NG.cohort.xlsx"))
meth.data.full <- readRDS(file.path(data.dir, "450k.meth.table.bcells.RDS"))

#check for the presence of all B cells IDs
normal.Bcells <-  c("cg", "NBCs", "GCFs", "loMBCs", "intMBCs","MGZs","hiMBCs")
meth.data.Bcell <- meth.data.full[, normal.Bcells]
  




############################################################
# methylation data - merged B cells and CLLs
merged.methylation.data <- merge(meth.data.Bcell, meth.data.CLL,  by.x="cg", by.y="cg")
rownames(merged.methylation.data) <- merged.methylation.data$cg



########################################################
## data processing (transpose the matrix)
meth.data.without.rownames <- merged.methylation.data[,-1]
transposed.meth.matrix <- t(meth.data.without.rownames)




########################################################
## Calculate Euclidean and Manhattan distances. 
## Euclidean distance is included in case someone would like to use ot
## In Methyl-COOM method we are using Manhattan distances onlz

## Specify distance metrics
DISTANCE.METRICS <- c("euclidean" = "euclidean", "manhattan" = "manhattan")

# Step1: Calculate distances and infer phylogenetic tree using fastme.bal()
edge2node <- list()
best.tree <- list()
dist.mm <- list()

for(metric in 1:length(DISTANCE.METRICS)){
  names.metric <- DISTANCE.METRICS[metric]
  dist.mm[[metric]] <- stats::dist(transposed.meth.matrix, method=DISTANCE.METRICS[metric]) # computes distances between rows of the data matrix
  attr(dist.mm, "Labels") <- rownames(transposed.meth.matrix)
  #attr(dist.mm[[metric]], "Labels") <- as.character(1:40)
  best.tree[[metric]] <- fastme.bal(dist.mm[[metric]])
  edge2node[[metric]] <- apply(best.tree[[metric]]$edge, 1, min)
}

names(dist.mm) <- DISTANCE.METRICS #add names based on distance metrics
names(best.tree) <- DISTANCE.METRICS # add names based on distance metrics



# Step2: Test plot - unrooted tree for Manhattan distance
pdf(paste0(analysis.dir, "/01_phylogeny/test.Manhattan.distances.phylogeny.NG.cohort", Sys.Date(), ".pdf"))
plot.phylo(best.tree[["manhattan"]],
           type= "unrooted",
           cex = 0.5,
           use.edge.length=T,
           #lab4ut = "axial",
           edge.lty=1,
           no.margin=TRUE,
           direction="rightwards",
           rotate.tree=-10,
           align.tip.label=TRUE,
           label.offset = 0.6,
           edge.color = "black",
           edge.width =1,
           font=2,
           root.edge=T)
#nodelabels(frame="n")
dev.off()





# Step3: Assign the colors to the edges of the phylogenetic tree
# Green color (forestgreen) for normal B cells. 
# Palette of orange color for CLLs - colorRampPalette(c("gold", "goldenrod"))
tree <- best.tree[["manhattan"]]
edge.tip.labels=tree$tip.label[tree$edge[,2]]

bcells.edge <- which(edge.tip.labels  %in% c("NBCs", "GCFs", "loMBCs", "intMBCs", "MGZs", "hiMBCs"))
Oakes.egde <- which(edge.tip.labels %in% colnames(meth.data.without.rownames)[grep("CLL", colnames(meth.data.without.rownames) )])

e.colors <- vector(length=length(edge.tip.labels))
e.colors[bcells.edge] <- "forestgreen"
colfunc <- colorRampPalette(c("gold", "goldenrod"))
colors2 <- colfunc(length(Oakes.egde))

e.colors[Oakes.egde] <- colors2
e.colors[which(e.colors=="FALSE")] <- "black"



# Step4:  plot the unrooted tree with the color annotation for the egdes & tip labels
pdf(paste0(analysis.dir, "/01_phylogeny/Manhattan.Oakes.NG.cohort", Sys.Date(), ".pdf"))
plot.phylo(best.tree[["manhattan"]],
           type= "unrooted",
           cex = 0.5,
           use.edge.length=T,
           #lab4ut = "axial",
           edge.lty=1,
           no.margin=TRUE,
           direction="rightwards",
           rotate.tree=-30,
           align.tip.label=TRUE,
           label.offset = 0.6,
           edge.color = e.colors,
           edge.width =1,
           font=2,
           root.edge=T,
           show.tip.label = T)
#nodelabels(frame="n")
dev.off()


# Step5:  plot the unrooted tree with the color annotation for the egdes -but no- tip labels
pdf(paste0(analysis.dir, "/01_phylogeny/Manhattan.Oakes.NG.cohort.no.tip.labels", Sys.Date(), ".pdf"))
plot.phylo(best.tree[["manhattan"]],
           type= "unrooted",
           cex = 0.5,
           use.edge.length=T,
           #lab4ut = "axial",
           edge.lty=1,
           no.margin=TRUE,
           direction="rightwards",
           rotate.tree=-30,
           align.tip.label=TRUE,
           label.offset = 0.6,
           edge.color = e.colors,
           edge.width =1,
           font=2,
           root.edge=T,
           show.tip.label = F)
#nodelabels(frame="n")
dev.off()




# Step6:  Determine the full structure of MANHATTAN DISTANCE-BASED TREE 
selected_tree = best.tree[["manhattan"]]
full_structure_of_the_tree = phylo4(best.tree[["manhattan"]])
full_structure_of_the_tree_df = as(full_structure_of_the_tree, "data.frame")






# Step7:  Extract normal B cell differentiation axis (between NBC and hiMBC) 

## extract the position of nodes corresponding to normal B cell differentiation axis 
nbc_node = full_structure_of_the_tree_df$node[full_structure_of_the_tree_df$label %in% "NBCs"] 
himbc_node = full_structure_of_the_tree_df$node[full_structure_of_the_tree_df$label %in% "hiMBCs"] 


## extract the normal B cell differentiation axis (full path on the phylogenetic tree between NBCs and hiMBCs)
differentiation_axis_initial = nodepath(selected_tree, from=nbc_node, to=himbc_node)
differentiation_axis = differentiation_axis_initial[c(-1, -length(differentiation_axis_initial))]






# Step8:  Determine B cell ancestor nodes (cell-of-origin) for CLL samples
# This step will assign the ancestor present on normal B cell differentiation axis
# for CLL samples


## determine cell-of-origin (COO) for each CLL sample.
## COO means the closest normal B cell ancestor for CLL case
samples_ancestor_annotation_df = data.frame(matrix(nrow=length(1:nrow(transposed.meth.matrix)), ncol=2)) 
colnames(samples_ancestor_annotation_df) = c("sample", "node")



for(rows in 1:nrow(samples_ancestor_annotation_df)) {
  
  print(rows)
  sample.label <-  full_structure_of_the_tree_df$label[rows]
  ancestor = full_structure_of_the_tree_df$ancestor[full_structure_of_the_tree_df$label %in% sample.label]
  
  while (!(ancestor %in% differentiation_axis)) {
    new.node <-  ancestor
    ancestor = full_structure_of_the_tree_df$ancestor[full_structure_of_the_tree_df$node %in% new.node]
    
  }
  samples_ancestor_annotation_df$sample[rows] = sample.label
  samples_ancestor_annotation_df$node[rows] = ancestor
}




## determine the location of COO node for each CLL (location of the node for each CLL sample)
bcell.maturation.axis.nodes <- full_structure_of_the_tree_df[full_structure_of_the_tree_df$node.type=="internal", ]
bcell.maturation.axis.nodes <- bcell.maturation.axis.nodes[bcell.maturation.axis.nodes$node %in%  differentiation_axis, ]

bcell.maturation.axis.nodes$edge.length[1] = full_structure_of_the_tree_df$edge.length[full_structure_of_the_tree_df$label %in% "NBCs"]

samples.to.exclude <- c("NBCs", "GCFs", "loMBCs", "intMBCs", "MGZs", "hiMBCs")
samples.and.distances <- samples_ancestor_annotation_df[!(samples_ancestor_annotation_df$sample %in% samples.to.exclude), ]


for (samples in 1:nrow(samples.and.distances)) {
  
  print(samples)
  patient.id = samples.and.distances$sample[samples]
  
  node.location = samples.and.distances$node[samples.and.distances$sample %in% patient.id]
  start = full_structure_of_the_tree_df$ancestor[full_structure_of_the_tree_df$label %in% "NBCs"]
  start.node <- which(bcell.maturation.axis.nodes$node %in% start)
  stop = node.location
  stop.node <- which(bcell.maturation.axis.nodes$node == stop)
  
  sum.of.nodes = sum(bcell.maturation.axis.nodes$edge.length[start.node:stop.node])
  samples.and.distances$location[samples] <- sum.of.nodes
  
}

Distance.object <- samples.and.distances
DIstance.object.CLL <- Distance.object




####################################################
## save the object with the cell-of-origin assignment
saveRDS(Distance.object, paste0(analysis.dir, "/01_phylogeny//COO.assignment.Chris.NGcohort.on.linear.Bcell.specifics.CpGs", Sys.Date(), ".RDS"))







# Step9:  Determine the location on the phylogenetci tree for normal B cells
# calculate distances. intMBCs and MGZs will be treated as one subtype
# NBC will be treated as root. 


Bcells <- c("NBCs", "GCFs", "loMBCs", "intMBCs", "MGZs", "hiMBCs")
samples.and.distances <- samples_ancestor_annotation_df[samples_ancestor_annotation_df$sample %in% Bcells, ]


for (samples in 1:nrow(samples.and.distances)) {
  
  print(samples)
  patient.id = samples.and.distances$sample[samples]
  
  node.location = samples.and.distances$node[samples.and.distances$sample %in% patient.id]
  start = full_structure_of_the_tree_df$ancestor[full_structure_of_the_tree_df$label %in% "NBCs"]
  start.node <- which(bcell.maturation.axis.nodes$node %in% start)
  stop = node.location
  stop.node <- which(bcell.maturation.axis.nodes$node == stop)
  
  sum.of.nodes = sum(bcell.maturation.axis.nodes$edge.length[start.node:stop.node])
  samples.and.distances$location[samples] <- sum.of.nodes
  
}

Distance.object <- samples.and.distances





## path to hiMBC
edges.himbc <- Distance.object$location[Distance.object$sample=="hiMBCs"] + full_structure_of_the_tree_df$edge.length[full_structure_of_the_tree_df$label %in% "hiMBCs"]
Distance.object$location[Distance.object$sample=="hiMBCs"] <- edges.himbc

Distance.object <- Distance.object[-which(Distance.object$sample=="MGZs"),]

location.normal.B.cells <- data.frame(Bcells=c("NBC", "GCF", "loMBC", "intMBC_MGZ", "hiMBCs"),
                                      Location=Distance.object$location)

# save RDS object for B cells
saveRDS(location.normal.B.cells, paste0(analysis.dir, "/01_phylogeny//Manhattan.phylogeny.normal.B.cells.on.linear.Bcell.specifics.CpGs", Sys.Date(), ".RDS"))



# create a data frame storing phylogenetic distances
location.normal.B.cells <- data.frame(sample=c("NBC", "GCF", "loMBC", "intMBC_MGZ", "hiMBCs"),
                                      location=Distance.object$location)

location.normal.B.cells$location[location.normal.B.cells$sample=="NBC"] <- 0





# Step10: Combine two distance objects - for CLL and B cells
Distance.full <- rbind.data.frame(DIstance.object.CLL[,c(1,3)], location.normal.B.cells)
Distance.full$percent.programming <- Distance.full$location/Distance.full$location[112] *100
Distance.full <- Distance.full[order(Distance.full$percent.programming, decreasing = F), ]
write.xlsx(Distance.full, paste0(analysis.dir, "/01_phylogeny//Manhattan.phylogeny.normal.B.cells.OakesCLLs.on.linear.Bcell.specifics.CpGs", Sys.Date(), ".xlsx"))
