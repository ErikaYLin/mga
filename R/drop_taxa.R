#' @export
# drop_taxa1 <- function(x, # phyloseq object   ## MGA-class object ##
#                              pathogens = NULL # .csv object listing taxa to search for
#                              ) {
#
#   # Convert phyloseq tax_table to a data frame
#   pathogen <- phyloseq::tax_table(x)
#   pathogen <- as.data.frame(pathogen)
#
#   # Add genus to species name in Species column
#   spec.names <- pathogen$Species
#   pathogen$Species <- paste(pathogen$Genus, pathogen$Species, sep = " ")
#
#   # Add column to indicate possible pathogenicity (yes = 1, no = 0)
#   pathogen$pathogenic = 0
#   pathogen$pathogenic_taxon = NA
#
#   # Search for pathogens in taxonomy table
#   for (j in 1:nrow(pathogens)) {
#     # Search in Family
#     if ("Family" %in% pathogens$Taxon[j]){
#       for (i in 1:nrow(pathogen)) {
#         if (pathogens$Pathogen[j] %in% pathogen$Family[i]){
#           pathogen$pathogenic[i] = 1
#           pathogen$pathogenic_taxon[i] = pathogens$Pathogen[j]
#         }}
#       # Search in Genus
#     } else if ("Genus" %in% pathogens$Taxon[j]){
#       for (i in 1:nrow(pathogen)) {
#         if (pathogens$Pathogen[j] %in% pathogen$Genus[i]){
#           pathogen$pathogenic[i] = 1
#           if (is.na(pathogen$pathogenic_taxon[i])){
#             pathogen$pathogenic_taxon[i] = pathogens$Pathogen[j]
#           }}}
#       # search in Species
#     } else if ("Species" %in% pathogens$Taxon[j]){
#       for (i in 1:nrow(pathogen)) {
#         if (pathogens$Pathogen[j] %in% pathogen$Species[i]){
#           pathogen$pathogenic[i] = 1
#           if (is.na(pathogen$pathogenic_taxon[i])){
#             pathogen$pathogenic_taxon[i] = pathogens$Pathogen[j]
#           }}}
#     } else {
#       pathogen$pathogenic[i] = 0
#       pathogen$pathogenic_taxon[i] = NA}
#   }
#
#   # Remove genus from species names
#   pathogen$Species <- spec.names
#
#   # Update the tax_table
#   tax = phyloseq::tax_table(as.matrix(pathogen))
#   phyloseq::tax_table(ps.species) <- tax
#
#   # Extract pathogens only
#   pathogen.only <- pathogen[pathogen$pathogenic == 1,]
#   path.names <- row.names(pathogen.only)
#
#   # Filter pathogens from ASV table
#   pathogen.ps <- phyloseq::prune_taxa(path.names, ps.species)
#
#   # Return the filtered phyloseq object
#   return(pathogen.ps)
# }
#
#
#
# drop_taxa2 <- function(mga, # mga-class object ##
#                       taxa = NULL, # .csv object listing taxa to search for with columns named "taxon" and "group"
#                       network = TRUE,
#                       network.args = list(type = "taxa", # "samples"
#                                           distance = "jaccard", # "jaccard", "unifrac", "wunifrac", DPCoA", "jsd"
#                                           max.dist = 0.45, # default distance set by `make_network` function from `phyloseq` package is 0.4
#                                           keep.isolates = TRUE)) {
#
#
#   # Convert mga tax_table to a data frame
#   taxTab <- phyloseq::tax_table(mga$ps)
#   taxTab <- as.data.frame(taxTab)
#
#   # Add genus to species name in Species column
#   spec.names <- taxTab$Species
#   taxTab$Species <- paste(taxTab$Genus, taxTab$Species, sep = " ")
#
#   # Add column to indicate taxa to keep (yes = 1, no = 0)
#   taxTab$keep = 0
#   taxTab$keep_taxon = NA
#
#   # Search for taxa in taxonomy table
#   for (j in 1:nrow(taxa)) {
#     # Search in Family
#     if ("Family" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Family[i]){
#           taxTab$keep[i] = 1
#           taxTab$keep_taxon[i] = taxa$taxon[j]
#         }}
#       # Search in Genus
#     } else if ("Genus" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Genus[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # search in Species
#     } else if ("Species" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Species[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#     } else {
#       taxTab$keep[i] = 0
#       taxTab$keep_taxon[i] = NA}
#   }
#
#   # Remove genus from species names
#   taxTab$Species <- spec.names
#
#   # Update the tax_table
#   tax = phyloseq::tax_table(as.matrix(taxTab))
#   phyloseq::tax_table(mga$ps) <- tax
#
#   # Extract selected taxa only
#   keep.only <- taxTab[taxTab$keep == 1,]
#   keep.names <- row.names(keep.only)
#
#   # Filter taxa from ASV table
#   prune.ps <- phyloseq::prune_taxa(keep.names, mga$ps)
#
#   if (network) {
#
#     # Construct the network using `prune.ps`
#     prune.net <- phyloseq::make_network(prune.ps,
#                                   type = network.args$type,
#                                   distance = network.args$distance,
#                                   max.dist = network.args$max.dist,
#                                   keep.isolates = network.args$keep.isolates)
#     sampledata <- as.data.frame(phyloseq::sample_data(prune.ps)) # subset sample_info as a data frame
#
#     # Count vertices and edges in site network
#     v <- igraph::vcount(prune.net)
#     e <- igraph::ecount(prune.net)
#
#     connectivity <- e/v   # avg number associations between ASVs
#     connectance <- e/v^2  # fraction of edges out of all possible edges
#     # Degree of each node
#     degrees <- igraph::degree(prune.net, v = igraph::V(prune.net), mode = "all")
#
#   } else if (network == FALSE) {
#     sampledata <- as.data.frame(phyloseq::sample_data(prune.ps)) # subset sample_info as a data frame
#   }
#
#   # Aggregate duplicate species
#   species <- phyloseq::tax_glom(prune.ps, taxrank = 'Species', NArm = FALSE)
#
#   # Sum the presences in each sample for species richness
#   rich <- nrow(phyloseq::tax_table(species)) # species count in sample
#   # total individuals from all species in each sample
#   n <- apply(phyloseq::otu_table(prune.ps), 1, function(x) sum(x))
#   ASVs <- ncol(phyloseq::otu_table(prune.ps)) # total number of ASVs in site
#
#   # Alpha species diversity measures
#   # Shannon index
#   ps_alpha_div <- phyloseq::estimate_richness(prune.ps, split = TRUE, measure = "Shannon")
#   # Simpson index
#   ps_alpha_div2 <- phyloseq::estimate_richness(prune.ps, split = TRUE, measure = "Simpson")
#
#   # Phylogenetic diversity
#   # Sum all branch lengths in the phylogenetic tree
#   PD <- sum(phyloseq::tree_layout(phyloseq::phy_tree(prune.ps))$edgeDT[["edge.length"]])
#
#   if (network) {
#
#     # Display sample measures in a data frame
#     results.prune <- data.frame(ps_alpha_div,
#                                   ps_alpha_div2,
#                                   rich,
#                                   PD,
#                                   n.indivs = n,
#                                   ASVs,
#                                   vertices = v,
#                                   edges = e,
#                                   connectivity,
#                                   connectance)
#
#     results.prune <- cbind(results.prune, sampledata))
#     # results.samples <- dplyr::relocate(sample.ID, .before = Shannon) # move Sample.ID column to leftmost
#
#     degree.prune <- data.frame(phyloseq::otu_table(prune.ps)[1,], degrees, row.names = colnames(phyloseq::otu_table(prune.ps)))
#     colnames(degree.prune)[1] <- sampledata$sample.ID
#
#   } else if (network == FALSE) {
#
#     # Display sample measures in a data frame
#     results.prune <- data.frame(ps_alpha_div,
#                                   ps_alpha_div2,
#                                   rich,
#                                   PD,
#                                   n.indivs = n,
#                                   ASVs)
#
#     results.prune <- cbind(results.prune, sampledata)
#     # results.samples <- dplyr::relocate(sample.ID, .before = Shannon) # move Sample.ID column to leftmost
#
#   }
#
#   # list of all objects to return
#   ps_network <- list()
#   ps_network$seqtabAll <- mga$seqtabAll
#   ps_network$ASVtab <- mga$ASVtab
#   ps_network$taxTab <- mga$taxTab
#   ps_network$phylo_tree <- mga$phylo_tree
#   ps_network$sampledata <- mga$sampledata
#   ps_network$ps <- mga$ps
#   ps_network$results.samples <- mga$results.samples
#
#   if ("net" %in% mga) {
#     ps.network$net <- mga$net
#     ps.network$degree.samp <- mga$degree.samp
#   }
#
#   ps_network$prune.ps <- prune.ps
#   ps.network$results.prune <- results.prune
#
#   if (network) {
#     ps_network$prune.net <- prune.net
#     ps_network$degree.prune <- degree.prune
#   }
#
#   # Convert to mga class
#   ps_network <- new.mga(ps_network)
#
#   return(ps_network)
# }
#


drop_taxa <- function(mga, # mga-class object
                      taxa = NULL, # .csv object listing taxa to search for with columns named "taxon" and "group"
                      group.species = TRUE,
                      network = TRUE,
                      network.args = list(type = "taxa", # "samples"
                                          distance = "jaccard", # "jaccard", "unifrac", "wunifrac", DPCoA", "jsd"
                                          max.dist = 0.45, # default distance set by `make_network` function from `phyloseq` package is 0.4
                                          keep.isolates = TRUE)) {


  # Convert mga tax_table to a data frame
  taxTab <- phyloseq::tax_table(mga$ps)
  taxTab <- as.data.frame(taxTab)

  # Add genus to species name in Species column
  spec.names <- taxTab$Species
  taxTab$Species <- paste(taxTab$Genus, taxTab$Species, sep = " ")

  # Add column to indicate taxa to keep (yes = 1, no = 0)
  taxTab$keep = 0
  taxTab$keep_taxon = NA

  # Search for taxa in taxonomy table
  for (j in 1:nrow(taxa)) {
    # Search in Kingdom
    if ("Kingdom" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab)) {
        if (taxa$taxon[j] %in% taxTab$Kingdom[i]){
          taxTab$keep[i] = 1
          taxTab$keep_taxon[i] = taxa$taxon[j]
        }}
      # Search in Phylum
    } else if ("Phylum" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab)) {
        if (taxa$taxon[j] %in% taxTab$Phylum[i]){
          taxTab$keep[i] = 1
          if (is.na(taxTab$keep_taxon[i])){
            taxTab$keep_taxon[i] = taxa$taxon[j]
          }}}
      # Search in Class
    } else if ("Class" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab)) {
        if (taxa$taxon[j] %in% taxTab$Class[i]){
          taxTab$keep[i] = 1
          if (is.na(taxTab$keep_taxon[i])){
            taxTab$keep_taxon[i] = taxa$taxon[j]
          }}}
      # Search in Order
    } else if ("Order" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab)) {
        if (taxa$taxon[j] %in% taxTab$Order[i]){
          taxTab$keep[i] = 1
          if (is.na(taxTab$keep_taxon[i])){
            taxTab$keep_taxon[i] = taxa$taxon[j]
          }}}
      # Search in Family
    } else if ("Family" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab)) {
        if (taxa$taxon[j] %in% taxTab$Family[i]){
          taxTab$keep[i] = 1
          if (is.na(taxTab$keep_taxon[i])){
            taxTab$keep_taxon[i] = taxa$taxon[j]
          }}}
      # Search in Genus
    } else if ("Genus" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab)) {
        if (taxa$taxon[j] %in% taxTab$Genus[i]){
          taxTab$keep[i] = 1
          if (is.na(taxTab$keep_taxon[i])){
            taxTab$keep_taxon[i] = taxa$taxon[j]
          }}}
      # Search in Species
    } else if ("Species" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab)) {
        if (taxa$taxon[j] %in% taxTab$Species[i]){
          taxTab$keep[i] = 1
          if (is.na(taxTab$keep_taxon[i])){
            taxTab$keep_taxon[i] = taxa$taxon[j]
          }}}
    } else {
      taxTab$keep[i] = 0
      taxTab$keep_taxon[i] = NA}
  }

  # Remove genus from species names
  taxTab$Species <- spec.names

  # Update the tax_table
  tax = phyloseq::tax_table(as.matrix(taxTab))
  phyloseq::tax_table(mga$ps) <- tax

  # Extract selected taxa only
  keep.only <- taxTab[taxTab$keep == 1,]
  keep.names <- row.names(keep.only)

  # Filter taxa from ASV table
  mga$ps <- phyloseq::prune_taxa(keep.names, mga$ps)

  # Aggregate duplicate species
  ps.species <- phyloseq::tax_glom(mga$ps, taxrank = 'Species', NArm = FALSE)

  if (group.species) {
  mga$ps.species <- ps.species
  }

  if (network) {

    # Construct the network using `prune.ps`
    net <- phyloseq::make_network(if (group.species){ps.species} else{mga$ps},
                                  type = network.args$type,
                                  distance = network.args$distance,
                                  max.dist = network.args$max.dist,
                                  keep.isolates = network.args$keep.isolates)
    sampledata <- as.data.frame(phyloseq::sample_data(mga$ps)) # subset sample_info as a data frame

    # Count vertices and edges in site network
    v <- igraph::vcount(net)
    e <- igraph::ecount(net)

    connectivity <- e/v   # avg number associations between ASVs
    connectance <- e/v^2  # fraction of edges out of all possible edges
    # Degree of each node
    degrees <- igraph::degree(net, v = igraph::V(net), mode = "all")

  } else if (network == FALSE) {
    sampledata <- as.data.frame(phyloseq::sample_data(mga$ps)) # subset sample_info as a data frame
  }

  # Sum the presences in each sample for species richness
  rich <- nrow(phyloseq::tax_table(ps.species)) # species count in sample
  # total individuals from all species in each sample
  n <- apply(phyloseq::otu_table(mga$ps), 1, function(x) sum(x))
  ASVs <- ncol(phyloseq::otu_table(mga$ps)) # total number of ASVs in site

  # Alpha species diversity measures
  # Shannon index
  ps_alpha_div <- phyloseq::estimate_richness(if (group.species){ps.species} else{mga$ps}, split = TRUE, measure = "Shannon")
  # Simpson index
  ps_alpha_div2 <- phyloseq::estimate_richness(if (group.species){ps.species} else{mga$ps}, split = TRUE, measure = "Simpson")

  # Phylogenetic diversity
  # Sum all branch lengths in the phylogenetic tree
  PD <- sum(phyloseq::tree_layout(phyloseq::phy_tree(if (group.species){ps.species} else{mga$ps}))$edgeDT[["edge.length"]])

  if (network) {

    # Display sample measures in a data frame
    results.samples <- data.frame(ps_alpha_div,
                                  ps_alpha_div2,
                                  rich,
                                  PD,
                                  n.indivs = n,
                                  ASVs,
                                  vertices = v,
                                  edges = e,
                                  connectivity,
                                  connectance)

    results.samples <- cbind(results.samples, sampledata)
    # results.samples <- dplyr::relocate(sample.ID, .before = Shannon) # move Sample.ID column to leftmost

    degree.samp <- data.frame(if (group.species){otu_table(ps.species)[1,]} else{phyloseq::otu_table(mga$ps)[1,]}, degrees,
                              if (group.species){row.names = colnames(otu_table(ps.species))} else{row.names = colnames(phyloseq::otu_table(mga$ps))})
    colnames(degree.samp)[1] <- sampledata$sample.ID

  } else if (network == FALSE) {

    # Display sample measures in a data frame
    results.samples <- data.frame(ps_alpha_div,
                                  ps_alpha_div2,
                                  rich,
                                  PD,
                                  n.indivs = n,
                                  ASVs)

    results.samples <- cbind(results.samples, sampledata)
    # results.samples <- dplyr::relocate(sample.ID, .before = Shannon) # move Sample.ID column to leftmost

  }

  # Return filtered mga-class object
  mga$ASVtab <- as.integer(phyloseq::otu_table(mga$ps))
  mga$taxTab <- keep.only
  mga$results.samples <- results.samples

  if (network) {
    mga$net <- net
    mga$degree.samp <- degree.samp
  }

  # Convert to mga-class object
  mga <- new.mga(mga)

  return(mga)
}
