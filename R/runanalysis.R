#' Do phydo analysis
#'
#' @param taxon Clade of interest
#' @return A list with all output
#' @export
run_phydo <- function(taxon) {
	wikipedia_summary <- wikipedia_pics <- datelife_biggest <- all_species <- misse_results <- pubmed_all <- genbank_count_by_gene <- genbank_count <- otol <- eol <- eol_tbl <- location_realm_biome <- NULL # very stupid way to initialize
  try(wikipedia_summary <- get_wikipedia_summary(taxon))
  try(wikipedia_pics <- get_wikipedia_pics(taxon))
  try(datelife_biggest <- get_datelife_biggest(taxon))
  try(all_species <- get_descendant_species(taxon))
  try(misse_results <- run_misse(datelife_biggest, all_species))
  try(pubmed_all <- get_pubmed(taxon, search.string=""))
  try(genbank_count_by_gene <- get_genbank_count_by_gene(taxon))
  try(genbank_count <- get_genbank_count(taxon))
  try(otol <- get_otol(taxon))
  try(eol <- get_eol(taxon))
  #eol_tbl <- eol_traits2(eol)
  try(location_realm_biome <- get_location_realm_biome(taxon))
  #return(list(wikipedia_summary=wikipedia_summary, datelife_biggest=datelife_biggest, pubmed=pubmed, genbank_count_by_gene=genbank_count_by_gene, genbank_count=genbank_count , otol=otol, location_realm_biome=location_realm_biome, eol=eol, eol_tbl=eol_tbl))
  return(list(wikipedia_summary=wikipedia_summary, wikipedia_pics=wikipedia_pics, datelife_biggest=datelife_biggest, all_species=all_species, misse_results=misse_results, genbank_count_by_gene=genbank_count_by_gene, genbank_count=genbank_count, otol=otol, eol=eol))

}

#' Create a file of results
#'
#' @param taxon Clade of interest
#' @param format Format: pdf or html
#' @param output_dir Where to put the output; by default, current working directory
#' @return Nothing, though a file is created in the current working directory
#' @export
#' @examples
#' render_phydo("Tyto")
render_phydo <- function(taxon, format="pdf", output_dir=getwd()) {
  rmarkdown::render(system.file("rmd", "summary.Rmd", package="phydo"), params=list(taxon=taxon),output_file=paste0("Report_",gsub(" ", "_",taxon), ".", format), output_dir=output_dir, encoding="UTF-8")
}


#' Create a file of results
#'
#' @param phydo_result The result of run_phydo
#' @param taxon Clade of interest
#' @return Nothing, though a file is created in the current working directory
#' @export
render_quarto <- function(phydo_result, taxon) {
  #quarto::quarto_render(system.file( "summary.qmd", package="phydo"), 
  quarto::quarto_render("test.qmd", 

	output_file=paste0("Report_",gsub(" ", "_",taxon), ".html"),
	execute_dir = tempdir(),
  	execute_params = list(
        taxon = taxon,
        all_results = jsonlite::serializeJSON(phydo_result)
    ))
}

#' Run biogeobears analyses
#'
#' Run a set of biogeobears analyses to estimate parameters of biogeographic models
#'
#' @param phy The phylo object
#' @param locations The states for biogeographic regions
#' @return A list with fits for each model
#' @export
run_biogeobears <- function(phy, locations) {

}

#' Run misse analysis
#' 
#' Run a misse analysis to estimate parameters of diversification models.
#' By passing in a vector of species names, it uses those to estimate the
#' sampling fraction based on how many are in the tree.
#' 
#' @param phy The phylo object
#' @param species A vector of species names in the clade, NOT just on the tree.
#' @return A list with class misse.fits objects
#' @export
run_misse <- function(phy, species=NULL) {
	if(!inherits(phy, "phylo")) {
		return(NULL)
	}
	misse_results <- NULL
	sampling_fraction <- 1
	if(!is.null(species)) {
		sampling_fraction <- ape::Ntip(phy)/length(species)
	}
	try({
		misse_results <- hisse::MiSSEGreedy(phy, f=sampling_fraction, generateMiSSEGreedyCombinations(max.param=max(3, min(52,ceiling(ape::Ntip(phy)/10))), fixed.eps=c(NA, 0, 0.9)))
	}, silent=TRUE)
	return(misse_results)
}

#' All pairwise correlations
#'
#' Use independent contrasts on each pair of traits. If a trait is missing for a taxon, drop that taxon for all analyses with that trait but not others. It runs cor.test on the positivized contrasts and returns the correlation, lower, and upper 95% values for the contrasts
#'
#' @param phy A phylo object
#' @param traits A data.frame of traits, with taxon names as rownames
#' @return A list with three objects: the correlations, lower, and upper.
#' @export
contrasts_correlations <- function(phy, traits) {
  correlations <- matrix(NA, ncol=ncol(traits), nrow=ncol(traits))
  rownames(correlations) <- colnames(traits)
  colnames(correlations) <- colnames(traits)
  correlations.lower <- correlations
  correlations.upper <- correlations
  for (i in sequence(nrow(correlations))) {
    for (j in sequence(ncol(correlations))) {
      if (j>i) {
        traits.local <- traits[,c(i,j)]
        traits.local <- traits.local[!is.na(traits.local[,1]),]
        traits.local <- traits.local[!is.na(traits.local[,2]),]

        pruned <- geiger::treedata(phy, traits.local, sort=TRUE, warnings=FALSE)
        pic.x <- ape::pic(pruned$data[,1], pruned$phy)
        pic.y <- ape::pic(pruned$data[,2], pruned$phy)
        # positivize the values
        pic.y[which(pic.x<0)] <- -pic.y[which(pic.x<0)]
        pic.x <- abs(pic.x)
        result <- cor.test(pic.x, pic.y)
        correlations[i,j] <- result$estimate
        correlations.lower[i,j] <- result$conf.int[1]
        correlations.upper[i,j] <- result$conf.int[2]
      }
    }
  }
  return(list(correlations=correlations, lower=correlations.lower, upper=correlations.upper))
}

#' Reconstruct an apomorphy down a tree using parsimony
#' 
#' Given a phylogeny, tip states, and which state is the apomorphy, reconstruct the apomorphy down the tree.
#' @param phy A phylo object
#' @param data A phyDat class object containing the states
#' @param apomorphy The state that is the apomorphy
#' @return A tree with node.labels with the reconstructed states
#' @export
reconstruct_apomorphy <- function(phy, data, apomorphy) {
	phy <- ape::multi2di(phy) # needs recon on fully resolved tree
	if(!ape::is.rooted(phy)) {
		stop("Tree must be rooted")
	}
	recon <- phangorn::ancestral.pars(phy, data, "MPR")
}