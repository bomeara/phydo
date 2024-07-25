#' Do phydo analysis
#'
#' @param taxon Clade of interest
#' @return A list with all output
#' @export
run_phydo <- function(taxon) {
	phydo_data <- get_phydo_data(taxon)
	misse_results <- NULL
  	try(misse_results <- run_misse(phydo_data$datelife_biggest, phydo_data$all_species))
	phydo_results <- phydo_data
	phydo_results$misse_results <- misse_results
	corhmm_results <- list()
	try({corhmm_results <- run_corHMM(phydo_results$datelife_biggest, phydo_results$all_traits_qualitative)})
	phydo_results$corhmm_results <- corhmm_results
	phydo_results$best_corhmm <- list()
	try({phydo_results$best_corhmm <- get_best_corhmm(corhmm_results)})
	phydo_results$contrasts_correlations <- NULL
	try({phydo_results$contrasts_correlations <- contrasts_correlations(phydo_results$datelife_biggest, phydo_results$all_traits_numeric)})
  	return(phydo_results)
}

#' Get phydo raw data
#' 
#' @param taxon Clade of interest
#' @return A list with all output
#' @export
get_phydo_data <- function(taxon) {
	wikipedia_summary <- wikipedia_pics <- datelife_biggest <- all_species <- pubmed_all <- genbank_count_by_gene <- genbank_count <- eol <- eol_tbl <- location_realm_biome <- gbif <- eol_tbl <- NULL # very stupid way to initialize
  try(wikipedia_summary <- get_wikipedia_summary(taxon))
  try(wikipedia_pics <- get_wikipedia_pics(taxon))
  try(datelife_biggest <- get_datelife_biggest(taxon))
  try(all_species <- get_descendant_species(taxon))
  #try(misse_results <- run_misse(datelife_biggest, all_species))
  try(pubmed_all <- get_pubmed(taxon, search.string=""))
  try(genbank_count_by_gene <- get_genbank_count_by_gene(taxon))
  try(genbank_count <- get_genbank_count(taxon))
  try(eol <- get_eol(taxon))
  try(eol_tbl <- eol_traits2(eol))
  try(gbif <- gbif_taxon_query(taxon))
  try(location_realm_biome <- get_location_realm_biome(gbif))
  phydo_data <- list(wikipedia_summary=wikipedia_summary, wikipedia_pics=wikipedia_pics, datelife_biggest=datelife_biggest, all_species=all_species,  location_realm_biome=location_realm_biome, genbank_count_by_gene=genbank_count_by_gene, genbank_count=genbank_count, pubmed_all=pubmed_all, eol=eol,eol_tbl=eol_tbl, gbif=gbif)
  phydo_data <- aggregate_phydo_traits(phydo_data)
  return(phydo_data)
}

remove_nchar_zero <- function(x) {
  x[!nchar(x)==0]
}

aggregate_phydo_traits <- function(phydo_data) {
	location_aggregated <- phydo_data$location_realm_biome |> dplyr::group_by(species) |> 
		dplyr::summarize(
			realm=paste(sort(unique(remove_nchar_zero(realm))), collapse="; "), 
			biome=paste(sort(unique(remove_nchar_zero(biome))), collapse="; "),
			continent=paste(sort(unique(remove_nchar_zero(continent))), collapse="; "),
			countryCode=paste(sort(unique(remove_nchar_zero(countryCode))), collapse="; "),
			mean_lat=mean(decimalLatitude, na.rm=TRUE),
			mean_lon=mean(decimalLongitude, na.rm=TRUE),
			mean_elevation=mean(elevation, na.rm=TRUE),
			median_lat=median(decimalLatitude, na.rm=TRUE),
			median_lon=median(decimalLongitude, na.rm=TRUE),
			median_elevation=median(elevation, na.rm=TRUE),
			max_975_lat=quantile(decimalLatitude, 0.975, na.rm=TRUE),
			max_975_lon=quantile(decimalLongitude, 0.975, na.rm=TRUE),
			max_975_elevation=quantile(elevation, 0.975, na.rm=TRUE),
			min_025_lat=quantile(decimalLatitude, 0.025, na.rm=TRUE),
			min_025_lon=quantile(decimalLongitude, 0.025, na.rm=TRUE),
			min_025_elevation=quantile(elevation, 0.025, na.rm=TRUE)
		)
	all_traits <- dplyr::full_join(phydo_data$eol_tbl, location_aggregated, by="species")
	all_traits <- as.data.frame(subset(all_traits, nchar(species)>0))
	all_traits_continuous <- as.data.frame(apply(all_traits, 2, as.numeric))
	all_traits_continuous$species <- all_traits$species
	numeric_traits <- c()
	qualitative_traits <- c()
	for (i in 2:ncol(all_traits_continuous)) {
		if(all(is.na(all_traits_continuous[,i]))) {
			qualitative_traits <- c(qualitative_traits, i)
		} else {
			numeric_traits <- c(numeric_traits, i)
		}
	}
	all_traits_qualitative <- all_traits[,c(1,qualitative_traits)]
	all_traits_numeric <- all_traits_continuous[,c(1,numeric_traits)]
	phydo_data$all_traits_qualitative <- all_traits_qualitative
	phydo_data$all_traits_numeric <- all_traits_numeric
	
	return(phydo_data)
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
		#misse_results <- hisse::MiSSEGreedy(phy, f=sampling_fraction, hisse::generateMiSSEGreedyCombinations(max.param=max(3, min(52,ceiling(ape::Ntip(phy)/10))), fixed.eps=c(NA, 0, 0.9)))
		misse_results <- hisse::MiSSEGreedy(phy, f=sampling_fraction, hisse::generateMiSSEGreedyCombinations(max.param=52, turnover.tries=sequence(5), eps.tries=sequence(5), fixed.eps=c(NA, 0, 0.9)))
	}, silent=TRUE)
	return(misse_results)
}

#' All pairwise correlations
#'
#' Use independent contrasts on each pair of traits. If a trait is missing for a taxon, drop that taxon for all analyses with that trait but not others. It runs cor.test on the positivized contrasts and returns the correlation, lower, and upper 95% values for the contrasts
#'
#' @param phy A phylo object
#' @param traits A data.frame of traits with first column as species
#' @return A list with three objects: the correlations, lower, and upper.
#' @export
contrasts_correlations <- function(phy, traits) {
	rownames(traits) <- traits[,1]
	traits <- traits[,-1]
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
        correlations.upper[j,i] <- result$conf.int[2]
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

#' Run a corHMM analysis
#' 
#' @param phy A phylo object
#' @param trait_data A data.frame with species in the first column and traits in the rest
#' @return A list with the results of the corHMM analysis
#' @export
run_corHMM <- function(phy, trait_data) {
	nchar <- ncol(trait_data)-1
	all_results <- list()
	if(!inherits(phy, "phylo")) {
		return(list())
	}
	for (i in 2:nchar) {
		print(paste0("Running character ", i, ": ", colnames(trait_data)[i+1]))
		focal_data <- trait_data[,c(1,i)]
		focal_data <- focal_data[!is.na(focal_data[,2]),]
		transformed_data <- focal_data
		if(length(unique(focal_data[,2]))>1) {
			global_states <- sort(unique(unlist(strsplit(as.character(focal_data[,2]), "; "))))
			print(paste0("There are ", length(global_states), " states in ", length(unique(focal_data[,2])), " combinations"))
			for (row_index in sequence(nrow(transformed_data))) {
				observed_states <- sort(unique(unlist(strsplit(as.character(focal_data[row_index,2]), "; "))))
				transformed_data[row_index,2] <- paste0(which(global_states %in% observed_states), collapse="&")
			}
			allowed_species <- intersect(phy$tip.label, focal_data$species)
			phy_focal <- ape::drop.tip(phy, setdiff(phy$tip.label, allowed_species))
			phy_focal <- ape::multi2di(phy_focal)
			transformed_data <- transformed_data[transformed_data$species %in% phy_focal$tip.label,]
			print(paste0("There are ", nrow(transformed_data), " species in the tree"))
			
			rate.cat_vector <- sequence(max(2, floor(ape::Ntip(phy_focal)/10)))
			model_vector <- c("ER", "ARD", "SYM")
			run_combinations <- expand.grid(rate.cat=rate.cat_vector, model=model_vector)
			if(length(global_states)>0.5*ape::Ntip(phy_focal)) {
				run_combinations <- subset(run_combinations, model=="ER") # too complex otherwise
			}
			
			if(ape::Ntip(phy_focal)>3 & length(global_states)<ape::Ntip(phy_focal)) {

				for (model_index in sequence(nrow(run_combinations))) {
					try({
						chosen_rate.cat <- run_combinations[model_index,1]
						chosen_model <- run_combinations[model_index,2]
						model_result <- corHMM::corHMM(phy_focal, transformed_data, model=chosen_model, rate.cat=chosen_rate.cat, get.tip.states=TRUE)
						model_result$chosen_model <- chosen_model
						model_result$global_states <- global_states
						model_result$character_index <- i
						model_result$character_name <- colnames(trait_data)[i+1]
						all_results[[length(all_results)+1]]<-  model_result
					}, silent=TRUE)
				}
			}
		}
	}
	return(all_results)
}

get_best_corhmm <- function(corhmm_result) {
	character_names <- sapply(corhmm_result, "[[", "character_name")
	all_best_models <- list()
	for (focal_name in sequence(unique(character_names))) {
		result_indices <- which(character_names==focal_name)
		local_results <- corhmm_result[result_indices]
		AICc_values <- sapply(local_results, "[[", "AICc")
		best_model_index <- which(AICc_values==min(AICc_values))
		best_model <- local_results[[best_model_index]]
		all_best_models[[length(all_best_models)+1]] <- best_model
	}
	return(all_best_models)
}

