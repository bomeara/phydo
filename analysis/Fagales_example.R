# mergingTreesWithMRP analysis

# devtools::install_github("phylotastic/datelife")
# devtools::install_github("dwbapst/paleotree", ref="developmentBranch")

# necessary to load ape though
library(ape)

# to test with otol Fagales tree
dq <- datelife::make_datelife_query("fagales",
	get_spp_from_taxon = TRUE)
tree_backbone<- datelife::get_otol_synthetic_tree(dq$cleaned_names,
	ott_ids = dq$ott_ids, resolve = FALSE)
is.binary(tree_backbone)
# and with Fagales PBDB taxon tree
faData <- paleotree::getCladeTaxaPBDB("Fagales")
# make the taxon tree
tree_secondary <- paleotree::makePBDBtaxonTree(
     taxaDataPBDB = faData,
     rankTaxon = "species",
     method = "parentChild"
     )
# look at stats of the trees
ape::Ntip(tree_backbone)
ape::Nnode(tree_backbone)
ape::Ntip(tree_secondary)
ape::Nnode(tree_secondary)


source("~//phydo//R//merging_trees_with_MRP.R")
mergedTree <- merging_trees_with_MRP(
	tree_backbone, tree_secondary, 
	reduce_collapse = TRUE,
	backbone_reweighting=10,
	trace=1)
mergedTree


pdf_path <- paste0("~//phydo//analysis/",
		"merged_Fagales_tree_",
		format(Sys.time(), "%m-%d-%y"),
		".pdf")

pdf(file = pdf_path,
	height = 50)
plot(mergedTree,
	cex=0.2,
	show.tip.label = TRUE,
	no.margin = TRUE)
dev.off()
# need to use path expand to convert ~
shell.exec(path.expand(pdf_path))

# fan

pdf_path <- paste0("~//phydo//analysis/",
		"merged_Fagales_tree_fan_",
		format(Sys.time(), "%m-%d-%y"),
		".pdf")
pdf(file = pdf_path,
	 height = 30, width = 30)
mergedTreeBrlen <- ape::compute.brlen(mergedTree)
plot(mergedTreeBrlen, cex = 0.2, type = "fan")
dev.off()
# need to use path expand to convert ~
shell.exec(path.expand(pdf_path))


########################################################

# then get the strict consensus?
# strictMerged <- consensus(mergedTrees, p = 1, check.labels = TRUE)
# plot(strictMerged)


