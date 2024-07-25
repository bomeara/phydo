

#' Query gbif to return latitude and longitude for members of a clade of interest
#'
#' @param query name of the clade of interest
#' @return A data frame (tibble) of key, scientific name, decimal latitude, and decimal longitude
#' @export
gbif_taxon_query <- function (query){
  key <- rgbif::name_backbone(query)
  gbif_download <- rgbif::occ_download(
	rgbif::pred('taxonKey', unname(key[1,'usageKey'])),
	rgbif::pred("hasGeospatialIssue", FALSE),
	rgbif::pred("hasCoordinate", TRUE),
	rgbif::pred_or(
    	rgbif::pred_not(rgbif::pred_in("establishmentMeans",c("MANAGED","INTRODUCED"))),
    	rgbif::pred_isnull("establishmentMeans")
    )
  )
  rgbif::occ_download_wait(gbif_download)
  dat <- rgbif::occ_download_get(gbif_download, path=tempdir()) |> rgbif::occ_download_import()
  return(dat)
}

# e.g. gbif_clade_query("fagales", "class")


#' Query gbif to return latitude and longitude for a vector of species
#'
#' @param species a vector of species names
#' @param gbif_limit Maximum number of records to return (hard limit is 200000)
#' @return A data frame (tibble) of key, scientific name, decimal latitude, and decimal longitude
#' @export
gbif_species_query <- function (species, gbif_limit=200000){
  all.records <- data.frame()
  for (species_index in seq_along(species)) {
    all.records <- plyr::rbind.fill(all.records, rgbif::occ_search(scientificName = species[species_index], fields="minimal", limit = gbif_limit)$data)
  }
  return(all.records)
}

# e.g.  species <- c("Puma concolor", "Myrmecocystus mexicanus", "Ursus arctos")
# gbif_species_query(species=species)

#' Query many sources at once using spocc to get latitude and longitude
#'
#' @param taxon The taxon to get. If a vector, loops over it.
#' @param limit The limit of records per site to get (default is maximum of any site)
#' @param sources Vector of sources (see ?spocc::occ)
#' @param has_coords Boolean for whether to only return records with longitude and latitude data
#' @param by_species Boolean: if TRUE, separates the taxon into species first (and searches for the higher level taxon as well)
#' @param verbose Boolean: if TRUE, print out progress
#' @param ... Other arguments to pass to spocc::occ
#' @return data.frame of results
#' @export
#' @examples
#' locations <- spocc_taxon_query("Myrmecocystus", limit=50)
spocc_taxon_query <- function(taxon, limit=10000, sources=c("gbif", "inat", "idigbio"), has_coords=TRUE, by_species=TRUE, verbose=TRUE, ...) {
  all.records <- data.frame()
  all.taxa <- c(taxon)
  if(by_species) {
    for (taxon_index in seq_along(taxon)) {
      all.taxa <- c(all.taxa,get_descendant_species(taxon[taxon_index]))
    }
  } else {
    all.taxa <- taxon
  }
  for (taxon_index in seq_along(all.taxa)) {
    local.records <- spocc::occ2df(spocc::occ(query=all.taxa[taxon_index], from=sources, limit=limit, has_coords=has_coords))
    if(verbose) {
      print(paste0("Now finished with ", nrow(local.records), " records for ", all.taxa[taxon_index], " which is taxon ", taxon_index, " of ", length(all.taxa), " taxa"))
    }
    if(nrow(local.records)>0) {
      local.records$taxon <- all.taxa[taxon_index]
      all.records <- plyr::rbind.fill(all.records, local.records)
    }
  }
  all.records$longitude <- as.numeric(all.records$longitude)
  all.records$latitude <- as.numeric(all.records$latitude)
  return(all.records)
}

get_all_species_only_from_datelife <- function(taxon) {
  all_species <- unname(datelife::get_all_species(taxon)$tnrs_names)
  all_species <- all_species[which(!grepl("_.+_", all_species))]
  return(all_species)
}

#' Clean locality information
#'
#' This uses the ropensci CoordinateCleaner package to clean up points.
#' @param locations Data.frame containing points (latitude and longitude, perhaps other data as columns)
#' @return Cleaned data.frame
#' @export
locality_clean <- function(locations) {
  locations <- CoordinateCleaner::cc_val(locations, lon="decimalLongitude", lat="decimalLatitude", value="clean")
  locations <- CoordinateCleaner::clean_coordinates(locations, lon="decimalLongitude", lat="decimalLatitude", species=NULL, tests=c( "centroids", "equal", "gbif", "institutions","zeros"), value="clean")
  locations$longitude <- locations$decimalLongitude
  locations$latitude <- locations$decimalLatitude
  return(locations)
}

#' Plot a map of points
#' @param gbif_points Data.frame containing points (latitude and longitude, perhaps other data as columns)
#' @return A ggplot2 map
#' 
#' Note that we will need to set zoom automatically, also add points, also make the ocean full of water.
#' @examples
#' gis_points <- spocc_taxon_query("Tyto", limit=100, sources="gbif")
#' plot_gis(gbif_points)
#' @export
plot_gis <- function(gbif_points) {
	min_lat <- min(gbif_points$latitude)
	max_lat <- max(gbif_points$latitude)
	min_lon <- min(gbif_points$longitude)
	max_lon <- max(gbif_points$longitude)
	
	med_bbox <- sf::st_bbox(c(xmin = min_lon, xmax = max_lon, ymin = min_lat, ymax = max_lat),
						crs = 4326)

	med_bbox_df <- data.frame(x = c(min_lon, max_lon),
							y = c(min_lat, max_lat))


	extent_zoomed <- raster::extent(med_bbox)
	prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

	elev_med <- elevatr::get_elev_raster(med_bbox_df, prj =  prj_dd, z = 3, clip = "bbox")
	elev_med_mat <- rayshader::raster_to_matrix(elev_med)

	base_map <- elev_med_mat |> 
	height_shade() |> 
	add_overlay(
		sphere_shade(elev_med_mat,
					texture = rayshader::create_texture(
					lightcolor = "#b8ff78",
					shadowcolor = "#193600",
					leftcolor = "#80d453",
					rightcolor = "#80d453",
					centercolor = "#568a27"),
					sunangle = 0,
					colorintensity = 5)
	)

	base_map |> plot_map()
	return(base_map)
}

#' Use azizka/speciesgeocodeR/ and WWF data to encode locations for habitat and biome
#' @param x Data.frame containing points (latitude and longitude, perhaps other data as columns)
#' @return data.frame with columns for habitat and biome.
#' @export
WWFload <- function(x = NULL) {
    if (missing(x)) {
        x <- getwd()
    }
    download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip", 
        destfile = file.path(x, "wwf_ecoregions.zip"))
    unzip(file.path(x, "wwf_ecoregions.zip"), exdir = file.path(x, "WWF_ecoregions"))
    file.remove(file.path(x, "wwf_ecoregions.zip"))
    wwf <- sf::st_read(file.path(x, "WWF_ecoregions", "official", 
        "wwf_terr_ecos.shp"))
    return(wwf)
}





#'
#' Uses info from http://omap.africanmarineatlas.org/BIOSPHERE/data/note_areas_sp/Ecoregions_Ecosystems/WWF_Ecoregions/WWFecoregions.htm to convert codes to more readable text
#'
#' @param locations Data.frame containing points (latitude and longitude, perhaps other data as columns)
#' @return data.frame with columns for habitat and biome.
#' @export
#' @examples
#' locations <- spocc_taxon_query("Myrmecocystus", limit=50)
#' locations <- locality_clean(locations)
#' locations <- locality_add_habitat_biome(locations)
#' print(head(locations))
locality_add_habitat_biome <- function(locations) {
  locations_spatial <- sf::st_as_sf(locations, coords=c("decimalLongitude", "decimalLatitude"), crs=4326)

  wwf <- WWFload(tempdir())
  sf::sf_use_s2(FALSE)
  wwf_flat<- sf::st_transform(wwf, 4326)
  mappedregions <- sapply(sf::st_intersects(locations_spatial,wwf_flat), function(z) if (length(z)==0) NA_integer_ else z[1])

  realms <- data.frame(code=c("AA", "AN", "AT", "IM", "NA", "NT", "OC", "PA"), realm=c("Australasia", "Antarctic", "Afrotropics", "IndoMalay", "Nearctic", "Neotropics", "Oceania", "Palearctic"), stringsAsFactors=FALSE)
  biomes <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga", "Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra", "Mediterranean Forests, Woodlands & Scrub", "Deserts & Xeric Shrublands", "Mangroves")
  locations$eco_name <- wwf$ECO_NAM[mappedregions]
  locations$biome <- biomes[wwf$BIOME[mappedregions]]
  locations$realm <- wwf$REALM[mappedregions]
  for(realm_index in seq_along(realms$code)) {
	locations$realm <- gsub(realms$code[realm_index], realms$realm[realm_index], locations$realm)
  }
  return(locations)
}

#' Aggregate count by category
#'
#' Get count (or frequency) of number of entries per taxon for each category. The return will be a matrix with rows equal to your taxa and columns all possible categories for the focal column, with entries being the number / frequency of records for that taxon for that category
#'
#' @param locations Data.frame of locations (i.e., from locality_add_habitat_biome)
#' @param focal Column name to aggregate data over.
#' @param group_by What column name to use for grouping
#' @param return_frequency Boolean; if TRUE, give frequency, not counts
#' @return A matrix of counts or frequencies
#' @export
#' @examples
#' locations <- spocc_taxon_query("Bubo", limit=500)
#' locations <- locality_clean(locations)
#' locations <- locality_add_habitat_biome(locations)
#' biome_counts <- aggregate_category(locations, focal="biome", group_by="taxon")
#' print(head(biome_counts))
#' realm_frequencies <- aggregate_category(locations, focal="realm", group_by="taxon", return_frequency=TRUE)
#' print(head(realm_frequencies))
aggregate_category <- function(locations, focal='realm', group_by = "taxon", return_frequency=FALSE) {
  categories <- sort(unique(locations[,focal]))
  taxa <- sort(unique(locations[,group_by]))
  result <- matrix(0, nrow=length(taxa), ncol=length(categories))
  rownames(result) <- taxa
  colnames(result) <- categories
  for (taxon_index in seq_along(taxa)) {
    for (category_index in seq_along(categories)) {
      result[taxon_index, category_index] <- nrow(subset(locations, locations[,focal]==categories[category_index] & locations[,group_by]==taxa[taxon_index]))
    }
  }
  if(return_frequency) {
    for (taxon_index in seq_along(taxa)) {
      result[taxon_index,] <- result[taxon_index,] / sum(result[taxon_index,])
    }
  }
  return(result)
}

#' Get all descendant species of the taxon
#'
#' Uses taxize, datelife, and sources of GBIF and OpenTree
#'
#' @param taxon Clade of interest
#' @return vector of species names
#' @export
get_descendant_species <- function(taxon, strict=TRUE) {
  species <- c()
  #col_id <- taxize::gnr_resolve(taxon, data_source_ids=1, ask=FALSE, fields="all", best_match_only=TRUE)

  #col_id <- taxize::get_colid_(taxon)[[1]]$id[1]
  #species <- taxize::downstream(col_id, downto = "species", db = "col")[[1]]$childtaxa_name

  try(gbif_id <- taxize::get_gbifid_(taxon)[[1]]$usagekey[1], silent=TRUE)
  try(species <- taxize::downstream(gbif_id, downto = "species", db = "gbif", limit=1000)[[1]]$name, silent=TRUE)
  try(species <- c(species, unname(datelife::get_all_descendant_species(taxon)$tnrs_names)), silent=TRUE)
  #try(col_id <- taxize::get_colid_(taxon)[[1]]$id[1])
  #try(species <- c(species,taxize::downstream(col_id, downto = "species", db = "col")[[1]]$childtaxa_name))
  try(species <- gsub("_", " ", species))
  try(species <- unique(species))
  return(species)
}


#' Get information on specimens from the paleobiology data base (largely for extinct taxa)
#'
#' @param taxon taxon name
#' @return a data frame of name, paleo longitude, paleo latitude, earliest record in ma (minma), latest record in ma (maxma), earliest interval, and latest interval
#' @export
pbdb_taxon_query <- function(taxon){
  pbdb_data <- read.csv(paste0("http://paleobiodb.org/",
    "data1.2/occs/list.txt?base_name=",utils::URLencode(taxon),"&level=3&show=paleoloc"),
    stringsAsFactors = FALSE)
  lat_long <- data.frame(pbdb_data$accepted_name, pbdb_data$paleolng, pbdb_data$paleolat, pbdb_data$max_ma, pbdb_data$min_ma, pbdb_data$early_interval, pbdb_data$late_interval)
  lat_long$searched_taxon <- taxon
  return(lat_long)
}
#
# remotes::install_github("ropensci/rfishbase")
# library("rfishbase")
# library("dplyr")
#
# dat <- distribution(species_list(Genus='Labroides'))
# dat <- dat[c("SpecCode", "Species", "NorthernLatitude", "NorthernLatitudeNS", "SouthernLatitude", "SouthernLatitudeNS", "WesternLongitude", "WesternLongitudeEW", "EasternLongitude", "EasternLongitudeEW")]


#' Get information on fish localities rfishbase
#'
#' @param genus fish genus name
#' @return a data frame of SpecCode, Species, NorthernLatitude, NorthernLatitudeNS, SouthernLatitude, SouthernLatitudeNS, WesternLongitude, WesternLongitudeEW, EasternLongitude, EasternLongitudeEW
#' @export

fishbase_genus_query <- function (genus) {
  dat <- rfishbase::distribution(species_list(Genus=genus))
  dat <- dat[c("SpecCode", "Species", "NorthernLatitude", "NorthernLatitudeNS", "SouthernLatitude", "SouthernLatitudeNS", "WesternLongitude", "WesternLongitudeEW", "EasternLongitude", "EasternLongitudeEW")]

  return(dat)
}

#' Get information on fish localities rfishbase
#'
#' @param species vector of fish species names
#' @return a data frame of SpecCode, Species, NorthernLatitude, NorthernLatitudeNS, SouthernLatitude, SouthernLatitudeNS, WesternLongitude, WesternLongitudeEW, EasternLongitude, EasternLongitudeEW
#' @export

fishbase_species_query <- function (species) {
  all.records <- data.frame()
  for (species_index in seq_along(species)) {
    all.records <- plyr::rbind.fill(all.records, rfishbase::distribution(species_list=species[species_index]))
  }
  return(all.records)
}

#' Get info from Open Tree of Life
#'
#' @param taxon The clade to investigate
#' @return list containing studies (a data.frame of info about that taxon. Columns include year of the study, number of taxa in the tree from that study, the study citation, and the study DOI) and ntaxa (the number of taxa in OpenTree for that taxon).
#' @export
#' @examples
#' info <- get_otol("Gallus")
#' histogram(info$studies$year) # Years in which chicken papers in OpenTree were published
get_otol <- function(taxon) {
  clade.info <- rotl::tnrs_match_names(taxon)
  clade.name <- clade.info$unique_name[1]
  id <- clade.info$ott_id[1]
  node.info <- rotl::tol_node_info(id)
  relevant.studies <- rotl::studies_find_trees(property="ot:ottTaxonName", value=clade.name)
  tree.info <- data.frame()
  all.trees <- rotl::list_trees(relevant.studies)
  for (study.index in sequence(nrow(relevant.studies))) {
    study.info <- rotl::get_publication(rotl::get_study_meta(relevant.studies$study_ids[study.index]))
    for (tree.index in sequence(length(all.trees[study.index]))) {
      phy <- NULL
      try(phy <- rotl::get_study_tree(study_id=relevant.studies$study_ids[study.index], tree_id = gsub('tree/', '',all.trees[[study.index]][tree.index])))
      local.result <- NULL
      if(!is.null(phy)) {
        try(local.result <- data.frame(Year=relevant.studies$study_year[study.index], Ntax=ape::Ntip(phy), Pub=study.info[1], DOI=attr(study.info, "DOI")))
      } else {
        try(local.result <- data.frame(Year=relevant.studies$study_year[study.index], Ntax=NA, Pub=study.info[1], DOI=attr(study.info, "DOI")))
      }
      if(!is.null(local.result)) {
        if(nrow(tree.info)==0) {
          tree.info <- local.result
        } else {
          tree.info <- rbind(tree.info, local.result)
        }
      }

    }
  }
  return(list(studies=tree.info, ntaxa=node.info$num_tips))
}

#' Information about the number of species in the clade in genbank
#'
#' @param taxon The clade of interest
#' @return The count of species in genbank
#' @export
#' @examples
#' taxon <- "Myrmecocystus"
#' print(paste("There are", get_genbank_count(taxon), "species of", taxon, "in GenBank"))
get_genbank_count <- function(taxon) {
  clade.name<- rotl::tnrs_match_names(taxon)$unique_name[1]
  genbank.species.query <- paste0(clade.name, '[subtree] AND species[rank] AND specified[prop]')
  genbank.species.count <-  rentrez::entrez_search(db="taxonomy", genbank.species.query, use_history=TRUE)$count
}

#' Get the count of number of sequences for different genes for this taxon
#'
#' This will use a few popular genes by default, but you can pass your own instead.
#'
#' @param taxon The clade of interest
#' @param focal.genes Vector of gene names
#' @return vector of counts by gene, with label being the given gene name
#' @export
#' @examples
#' get_genbank_count_by_gene("Myrmecocystus")
get_genbank_count_by_gene <- function(taxon, focal.genes=c("COI", "18S", "28S", "matk", "rbcl")) {
  clade.name<- rotl::tnrs_match_names(taxon)$unique_name[1]
  GetNucCount <- function(gene, taxon=clade.name) {
    gene.query <- paste0(taxon, '[organism] AND ',gene)
    Sys.sleep(3) #just to make sure we stay nice
    return(rentrez::entrez_search(db="nuccore", gene.query, use_history=TRUE)$count)
  }
  gene.count <- sapply(focal.genes, GetNucCount, taxon=clade.name) #make sure not to do mclapply or similar lest you
  return(gene.count)
}

#' Get information on the clade from pubmed
#'
#' This will give the total number of pubmed articles mentioning the taxon name and other search string and information on the most recent such papers.
#'
#' @param taxon The clade of interest
#' @param search.string Include spaces and AND and similar search elements
#' @param retmax how many papers to return
#' @return List with a count element (integer) and a recent.papers element (data.frame)
#' @export
#' @examples
#' taxon <- "Formicidae"
#' results <- get_pubmed(taxon)
#' print(paste("There are", results$count, "papers on", taxon, "and phylogeny"))
#' print(results$recent.papers[,c("Date", "FirstAuthor")])

get_pubmed <- function(taxon, search.string=' AND phylogeny',retmax=50) {
  clade.name<- rotl::tnrs_match_names(taxon)$unique_name[1]
  pubmed.query <- paste0(clade.name, "[TIAB] ", search.string)
  pubmed.result <- rentrez::entrez_search(db="pubmed", pubmed.query, use_history=TRUE, retmax=retmax)
  if(length(pubmed.result$id)>0 & !is.na(clade.name)) {
    pubmed.summaries <- rentrez::entrez_summary(db="pubmed", id=pubmed.result$id)
    pubmed.df <- data.frame(Date=rentrez::extract_from_esummary(pubmed.summaries, elements=c("sortpubdate")), FirstAuthor=rentrez::extract_from_esummary(pubmed.summaries, elements=c("sortfirstauthor")), Journal=rentrez::extract_from_esummary(pubmed.summaries, elements=c("fulljournalname")), Title=rentrez::extract_from_esummary(pubmed.summaries, elements=c("title")), row.names=NULL)
  } else {
    pubmed.result$count <- 0
    pubmed.df <- data.frame(Date=NA, FirstAuthor=NA, Journal=NA, Title=NA, row.names=NULL)
  }
  return(list(count=pubmed.result$count,   recent.papers =   pubmed.df ))
}

#' Get information on the count of papers from pubmed
#'
#' This will give the total number of pubmed articles mentioning the taxon name and other search string and information on the most recent such papers.
#'
#' @param taxon The clade of interest
#' @param search.string Include spaces and AND and similar search elements
#' @return List with a count element (integer) and a recent.papers element (data.frame)
#' @export
#' @examples
#' taxon <- "Formicidae"
#' results <- get_pubmed_count_only(taxon)
#' print(paste("There are", results, "papers on", taxon, "and phylogeny"))

get_pubmed_count_only <- function(taxon, search.string=' AND phylogeny') {
  clade.name<- rotl::tnrs_match_names(taxon)$unique_name[1]
  if(is.na(clade.name)) {
	return(0)
  }
  pubmed.query <- paste0(clade.name, "[TIAB] ", search.string)
  pubmed.result <- rentrez::entrez_search(db="pubmed", pubmed.query, use_history=TRUE)
  return(pubmed.result$count)
}

#' Get location, realm, and biome
#'
#' @param gbif_record Data from GBIF
#' @return list containing a data.frame of species and locations, a table of realms (biogeographic regions), and a table of biomes
#' @export
get_location_realm_biome <- function(gbif_record) {
  locations <- locality_clean(gbif_record)
  locations <- locality_add_habitat_biome(locations)
  #biome <- aggregate_category(locations, focal="biome", group_by="taxon")
  #realm <- aggregate_category(locations, focal="realm", group_by="taxon")
  return(locations)
}

#' Get biggest tree from datelife
#'
#' @param taxon Clade of interest
#' @return phylo object
#' @export
get_datelife_biggest <- function(taxon) {
  clade.name<- rotl::tnrs_match_names(taxon)$unique_name[1]
  all_species <- datelife::get_all_descendant_species(clade.name)
  datelife_biggest <- NULL
  try(datelife_biggest <- datelife::datelife_search(input=unname(all_species$tnrs_names), get_spp_from_taxon=FALSE, summary_format="phylo_biggest"))
  
  return(datelife_biggest)
}

#' Get all tree from datelife
#'
#' @param taxon Clade of interest
#' @return phylo object
#' @export
get_datelife_all <- function(taxon) {
  clade.name<- rotl::tnrs_match_names(taxon)$unique_name[1]
  all_species <- datelife::get_all_descendant_species(clade.name)
  datelife_all <- NULL
  try(datelife_biggest <- datelife::datelife_search(input=unname(all_species$tnrs_names), get_spp_from_taxon=FALSE, summary_format="phylo_all"))
  return(datelife_all)
}

#' Get Wikipedia summary
#'
#' @param taxon Clade of interest
#' @return text string of summary of page
#' @export
get_wikipedia_summary <- function(taxon) {
  URL <- paste0('https://en.wikipedia.org/w/api.php?format=json&action=query&prop=extracts&exintro&explaintext&redirects=1&titles=', utils::URLencode(taxon))
  return(jsonlite::fromJSON(URL)$query$pages[[1]]$extract)
}

#' Get Wikipedia images
#'
#' @param taxon 
#' @return images from Wikipedia
#' @export
get_wikipedia_pics <- function(taxon) {
  URL <- paste0('https://en.wikipedia.org/api/rest_v1/page/media-list/', utils::URLencode(taxon))
  image_URL <- (jsonlite::fromJSON(URL)$items$srcset)[[1]]$src[1]
  return(gsub('//', 'http://', image_URL))
}


#' Get Wikipedia images
#'
#' @param taxon 
#' @return images from Wikipedia
#' @export
get_combie_image_summary <- function(taxon){

  w1<-get_wikipedia_summary(taxon)
    
  #get the text from the body
#   html <- XML::htmlTreeParse(w1, useInternal = TRUE)
#   txt <- XML::xpathApply(html, "//body//text()[not(ancestor::script)][not(ancestor::style)][not(ancestor::noscript)]", xmlValue)
#   txt<-toString(txt)
  
  ww1 <-get_wikipedia_pics(taxon)
  
  wikiimage <- magick::image_read_svg(ww1, width = 350)

  #read file
  img<-readJPG("wikiimage")
  
  #get size
  h<-dim(img)[1]
  w<-dim(img)[2]
  
  #open new file for output
  jpg("out.jpg", width=w, height=h)
  par(mar=c(0,0,0,0), xpd=NA, mgp=c(0,0,0), oma=c(0,0,0,0), ann=F)
  plot.new()
  plot.window(0:1, 0:1)
  
  #fill plot with image
  usr<-par("usr")    
  rasterImage(img, usr[1], usr[3], usr[2], usr[4])
  
  #close image
  dev.off()
    
 return(total)
}



#' Get info from Encyclopedia of Life
#'
#' @param taxon The clade to investigate
#' @return data.frame with species name, trait category, trait value, and other information
#' @export
#' @examples
#' info <- get_eol("Gallus")
get_eol <- function(taxon) {
  descendants <- get_descendant_species(taxon)
  all_df <- data.frame()
  for (taxon_index in seq_along(descendants)) {
	  local_df <- data.frame()
	  try(local_df <- eol_data(descendants[taxon_index]), silent=TRUE)
	  Sys.sleep(1) #just to make sure we stay nice
	  if(nrow(local_df)>0) {
		  print(paste0("Found ", nrow(local_df), " EOL datapoints for ", length(unique(local_df$trait)), " traits for ", descendants[taxon_index]))
		  all_df <- rbind(all_df, local_df)
	  } else {
		 print(paste0("Found 0 EOL datapoints for ", descendants[taxon_index]))  
	  }
  }
  for (col_index in sequence(ncol(all_df))) {
	  all_df[,col_index] <- stringr::str_trim(all_df[,col_index])
  }
  return(all_df)
}

#' Aggregate pubmed data down OToLtree
#' 
#' @param taxon The clade to investigate
#' @param ignore_subspecies Should subspecies be ignored?
#' @param just_genera Should only genera be used?
#' @return List with tree and information about nodes and tips
#' @export
#' @examples
#' taxon <- "Pinaceae"
#' info <- get_pubmed_info_on_OTOL(taxon)
get_pubmed_info_on_OTOL <- function(taxon, ignore_subspecies=TRUE, just_genera=FALSE) {
	resolved_name <- rotl::tnrs_match_names(taxon)
	focal_tree <- rotl::tol_subtree(ott_id = resolved_name$ott_id, label_format="name")
	focal_tree <- reroot.tree(focal_tree)
	all_species <- focal_tree$tip.label
	if(ignore_subspecies) {
		bad_taxa <- all_species[which(grepl("_.+_", all_species))]
		if(length(bad_taxa)>0) {
			focal_tree <- ape::drop.tip(focal_tree, bad_taxa)
		}

	}
	if (just_genera) {
		focal_tree <- pruneToGenera(focal_tree)	
	}
	focal_tree <- reorder(focal_tree, order="postorder")
	focal_df <- as(as(focal_tree, "phylo4d"), "data.frame")
	focal_df$PubmedThisLevel <- NA
	for (i in sequence(nrow(focal_df))) {
		if (!is.na(focal_df$label[i])) {
			focal_df$PubmedThisLevel[i] <- get_pubmed_count_only(focal_df$label[i], search.string="")
			print(paste0("Found ", focal_df$PubmedThisLevel[i], " papers for ", focal_df$label[i]))
		}
	}
	
	focal_df$PubmedDescendantTipsOnlyLowerBound <- 0
	focal_df$PubmedDescendantTipsOrInternalLowerBound <- 0
	for (i in sequence(nrow(focal_df))) {
		
		if (focal_df$node.type[i]=="internal") {
			descendant_tips <- focal_df$node[which(focal_df$ancestor==focal_df$node[i])]
			focal_df$PubmedDescendantTipsOnlyLowerBound[i] <- max(focal_df$PubmedThisLevel[which(focal_df$node %in% descendant_tips)], na.rm=TRUE)
			if(!is.finite(focal_df$PubmedDescendantTipsOnlyLowerBound[i])) {
				focal_df$PubmedDescendantTipsOnlyLowerBound[i] <- 0
			}
			focal_df$PubmedDescendantTipsOrInternalLowerBound[i] <- max(c(focal_df$PubmedDescendantTipsOnlyLowerBound[i], focal_df$PubmedThisLevel[i]), na.rm=TRUE)
			
		}
	}
	return(list(tree=focal_tree, info=focal_df))
}


#' Plot pubmed data down OToLtree
#' 
#' This will use reconstruction of number of papers down a tree
#' @param pubmed_info_on_OTOL List with tree and information about nodes and tips, from get_pubmed_info_on_OTOL
#' @return Plot
#' @export
plot_pubmed_info_on_OTOL <- function(pubmed_info_on_OTOL) {
	p4d <- phylo4d(pubmed_info_on_OTOL$tree, all.data=pubmed_info_on_OTOL$info)	
	traits <- p4d@data
	traits <- subset(traits, node.type=="tip")
	rownames(traits) <- traits$label
	combined <- geiger::treedata(pubmed_info_on_OTOL$tree, traits)
	combined_data <- as.data.frame(combined$data)
	combined_data$PubmedThisLevel <- as.numeric(combined_data$PubmedThisLevel)
	combined_data$Log1pPubmedThisLevel <- log1p(combined_data$PubmedThisLevel)
	Log1pPubmedThisLevel <- combined_data$Log1pPubmedThisLevel
	names(Log1pPubmedThisLevel) <- combined_data$label
	PubmedThisLevel <- combined_data$PubmedThisLevel
	names(PubmedThisLevel) <- combined_data$label
	combined_tree <- ape::multi2di(combined$phy)
	combined_tree <- ape::compute.brlen(combined_tree)
	contmap_obj <- phytools::contMap(combined_tree, PubmedThisLevel, plot=FALSE)
	contmap_obj_viridis <- setMap(contmap_obj, c("gray", rev(viridis::viridis(100))))
	plot(
		contmap_obj_viridis, 
		type="fan", 
		legend = 0.7*max(nodeHeights(combined_tree)),
		fsize = c(0.5, 0.7)
	)

}