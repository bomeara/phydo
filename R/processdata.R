#' Continuous to discrete
#'
#' Convert a set of latitudes and longitudes to biogeographic regions
#'
#' @param latlon The data.frame of points, with species, latitude, and longitude columns
#' @return A data.frame as above but with a region column added
#' @export
convert_continuous_to_discrete <- function(latlon) {
}


#' Get a tree and data, resolve to same set of taxa, and make a shared object
#'
#' @param tree a phylo object
#' @param data A data.frame. Rownames are species names.
#' @param warnings Complain about mismatches if TRUE.
#' @return The function returns a list of two elements (phy and data) that are manipulated to include only those species found in both the tree and data supplied by the user. The class of the returned object is "phydo".
#' @export
match_data <- function(tree, data, warnings = FALSE) {
  dm = length(dim(data))
  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  if (is.factor(data)) {
    data <- as.matrix(data)
  }
  if (is.array(data) & length(dim(data)) == 1) {
    data <- as.matrix(data)
  }
  if (is.null(rownames(data))) {
    stop("names for 'data' must be supplied")
  } else {
    data.names <- rownames(data)
  }
  nc <- geiger::name.check(tree, data)
  if (is.na(nc[[1]][1]) | nc[[1]][1] != "OK") {
    if (length(nc[[1]] != 0)) {
      tree <- ape::drop.tip(tree, as.character(nc[[1]]))
      if (warnings) {
        warning(paste("The following tips were not found in 'data' and were dropped from 'phy':\n\t",
                      paste(nc[[1]], collapse = "\n\t"), sep = ""))
      }
    }
    if (length(nc[[2]] != 0)) {
      m <- match(data.names, nc[[2]])
      data = as.data.frame(data[is.na(m), ],stringsAsFactors=FALSE)
      data.names <- data.names[is.na(m)]
      if (warnings) {
        warning(paste("The following tips were not found in 'phy' and were dropped from 'data':\n\t",
                      paste(nc[[2]], collapse = "\n\t"), sep = ""))
      }
    }
  }
  order <- match(data.names, tree$tip.label)
  rownames(data) <- tree$tip.label[order]

  index <- match(tree$tip.label, rownames(data))
  index <- index[!is.na(index)]
  data <- as.data.frame(data[index, ], stringsAsFactors=FALSE)
  if (dm == 2) {
    data <- as.data.frame(data, stringsAsFactors=FALSE)
  }
  MatchReturn<-list(phy = tree, data = data)
  class(MatchReturn)<-c("list","phydo")
  return(MatchReturn)
}


#' Convert the phydo format to one usable for Beaulieu-style
#'
#' Many of the packages authored by Jeremy Beaulieu with collaborators have the first column used to hold taxon names rather than store them as row names. Bless his heart. This converts a phydo object, which has its data with taxon names as row names, into having a data object internally with a column for taxon names.
#'
#' @param phydo The phydo object
#' @export
#' @return a list of two objects, phy and data, with the data in the format preferred by Beaulieu
phydo_convert_to_Beaulieu_data <- function(phydo) {
  return(list(phy=phydo$phy, data=cbind.data.frame(taxa=rownames(match_data$data), phydo$data[,1:ncol(phydo$data)])))
}

#' Check a vector to see if the elements are numeric
#'
#' @param x Vector of class whatever
#' @return A boolean on whether the elements, other than NAs, are numeric
check_numeric_vector <- function(x) {
  x <- x[!is.na(x)]
  return(all(grepl("^[0-9\\.]{1,}$",x)))
}

#' Check a data.frame to see which columns are continuous
#'
#' Continuous columns are not integers or have letters. They may have NAs
#'
#' @param x data.frame of traits
#' @return A vector of booleans
#' @export
check_continuous <- function(x) {
  continuous <- sapply(x, check_numeric_vector)
  numeric_cols <- which(continuous)
  for (i in seq_along(numeric_cols)) {
    local_vector <- x[,numeric_cols[i]]
    local_vector <- local_vector[!is.na(local_vector)]
    local_vector <- as.numeric(local_vector)
    if(all(local_vector==round(local_vector))) {
      continuous[numeric_cols[i]] <- FALSE
    }
  }
  return(continuous)
}

#' Return a phydo object with only the selected kind of traits
#'
#' @param phydo a phydo class object, with tree and data
#' @param keep Whether you want to keep continuous or discrete data
#' @return A phydo object with the tree and only the chosen kind of data
#' @export
phydo_drop_type <- function(phydo, keep=c("continuous", "discrete")) {
  keep_continuous <- TRUE
  if(keep=="discrete") {
    keep_continuous <- FALSE
  }

  column_check <- check_continuous(phydo$data)
  if(!keep_continuous) {
    column_check <- !column_check
  }
  phydo$data  <- phydo$data[, column_check, drop=FALSE]
  return(phydo)
}

#' Try geiger's fitContinuous or fitDiscrete on all the continuous datasets using all geiger models
#'
#' Note that this ignores SE at the moment. If a character is missing information for a taxon, that taxon is deleted for just that character only.
#'
#' @param phydo a phydo class object, with tree and data
#' @param models which models to use. If NULL, uses all available
#' @param keep continuous or discrete data
#' @param ncores how many cores to use; if NULL, detects automatically
#' @return A two dimensional list. The first dimension is model, the second is character
#' @export
#' @examples
#' data(geospiza,package="geiger")
#' phydo <- match_data(geospiza$phy, geospiza$dat)
#' results <- phydo_fitContinuous(phydo)
#' # Look at model for OU for character 1
#' print(results[["OU"]][["wingL"]])
phydo_fitGeiger <- function(phydo, models=NULL, keep=c("continuous", "discrete"), ncores=NULL){

  if(keep=="continuous") {
    if(is.null(models)) {
      models=c("BM","OU","EB","trend","lambda","kappa","delta","drift","white")
    }
    geiger_dat <- phydo_drop_type(phydo, keep="continuous")
    fitGeiger <- geiger::fitContinuous
  } else {
    if(is.null(models)) {
      models=c("ER","SYM","ARD","meristic")
    }
    geiger_dat <- phydo_drop_type(phydo, keep="discrete")
    fitGeiger <- geiger::fitDiscrete
  }

  if(is.null(colnames(geiger_dat$data))) {
    colnames(geiger_dat$data) <- paste0(keep, "_", sequence(ncol(geiger_dat$data)))
  }
  fitGeigerResList <- list()
  for(model_index in seq_along(models)){
    for (char_index in sequence(ncol(geiger_dat$data))) {
      sliced_data <- geiger_dat
      sliced_data$data <- sliced_data$data[,char_index, drop=TRUE]
      names(sliced_data$data) <- rownames(geiger_dat$data)
      sliced_data$data <- sliced_data$data[!is.na(sliced_data$data)]
      sliced_data_cleaned <- match_data(sliced_data$phy, sliced_data$data)
      fitGeigerResList[[models[model_index]]][[colnames(geiger_dat$data)[char_index]]] <- fitGeiger(phy = sliced_data_cleaned$phy, dat = sliced_data_cleaned$data, model = models[model_index], ncores=ncores)
    }
  }
  return(fitGeigerResList)
}

#
# #' Try geiger's fitContinuous on their discrete datasets using all geiger models
# #'
# #' @param phydo a phydo class object, with tree and datasets
# #' @param models which models to used
# #' @param ncores how many cores to use; if NULL, detects automatically
# #' @return A two dimensio list. The first dimension is model, the second is character
# #' @export
#
# phydo_fitDiscrete <- function (phydo, models = c("ER","SYM","ARD","meristic"), ncores=NULL){
#   DiscDat<-phydo_drop_type(phydo, keep="discrete")
#   fitDiscreteResList <- list()
#   for (model_index in seq_along(models)){
#     for (char_index in sequence(ncol(DiscDat$data))) {
#       sliced_data <- DiscDat
#       sliced_data$data <- sliced_data$data[,char_index, drop=TRUE]
#       names(sliced_data$data) <- rownames(DiscDat$data)
#       sliced_data$data <- sliced_data$data[!is.na(sliced_data$data)]
#       sliced_data_cleaned <- match_data(sliced_data$phy, sliced_data$data)
#       fitDiscreteResList[[models[model_index]]][[colnames(DiscDat$data)[char_index]]]<-geiger::fitDiscrete(phy = sliced_data_cleaned$phy, dat = sliced_data_cleaned$data, model = models[model_index], ncores = ncores)
#     }
#   }
#   return(fitDiscreteResList)
# }

#' Root a tree
#' 
#' From Jeremy Beaulieu, slight mods from Brian O'Meara
#' @param tree a phylo object
#' @return A rooted tree
#' @export
reroot.tree <- function(tree){
	tree <- phytools::bind.tip(tree, tip.label="outgroup", where=ape::Ntip(tree)+1)
	p<-root(tree,outgroup="outgroup",resolve.root=TRUE)
	bb <- which(p$tip.label=="outgroup")
	if(is.null(p$edge.length)) {
		p$edge.length <- rep(0, length(p$edge[,1]))
	}
	ee <- p$edge.length[which(p$edge[,2]==bb)]
	p$edge.length[which(p$edge[,2]==bb)] <- 0
	cc <- p$edge[which(p$edge[,2]==bb),1]
	dd <- setdiff(p$edge[which(p$edge[,1]==cc),2],bb)

	p$edge.length[which(p$edge[,2]==dd)]<-p$edge.length[which(p$edge[,2]==dd)]+(ee/2)
	p$edge.length[which(p$edge[,2]==bb)] <- (ee/2)
	#p <- ape::drop.tip(p, "outgroup", trim.internal=FALSE)
	#p <- ladderize(p)
	return(p)
}