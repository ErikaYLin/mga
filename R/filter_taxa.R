#' @export
# Taxonomic filtering and pruning
filter_taxa <- function(mga, # mga-class object
                        keep = TRUE, # keep only selected taxa in output
                        drop = TRUE, # drop only selected taxa in output
                        taxa = NULL, # .csv object listing taxa to search for with columns named "taxon" and "group"
                        group.species = TRUE,
                        network = TRUE,
                        network.args = list(type = "taxa", # "samples"
                                            distance = "jaccard", # "jaccard", "unifrac", "wunifrac", DPCoA", "jsd"
                                            max.dist = 0.45, # default distance set by `make_network` function from `phyloseq` package is 0.4
                                            keep.isolates = TRUE)) {

  if (keep) {
    # New mga object for taxa kept
    mga_keep <- mga
  }
  if (drop) {
    mga_drop <- mga
  }

  # Convert mga ps,species tax_table to a data frame
  taxTab <- phyloseq::tax_table(mga$ps.species)
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
  if (keep) {
    phyloseq::tax_table(mga_keep$ps.species) <- tax
  }
  if (drop) {
    phyloseq::tax_table(mga_drop$ps.species) <- tax
  }

  # Extract selected taxa only
  if (keep) {
    keep.only <- taxTab[taxTab$keep == 1,]
    keep.names <- row.names(keep.only)
  }
  if (drop) {
    drop.only <- taxTab[taxTab$keep == 0,]
    drop.names <- row.names(drop.only)
  }

  if (1 %in% taxTab$keep) {

    message("Selected taxa have been identified in mga object. Now filtering.")

  # Filter taxa from ASV table
  if (keep) {
    mga_keep$ps.species <- phyloseq::prune_taxa(keep.names, mga_keep$ps.species)
  }
  if (drop) {
    mga_drop$ps.species <- phyloseq::prune_taxa(drop.names, mga_drop$ps.species)
  }


  # Convert mga ps tax_table to a data frame
  taxTab2 <- phyloseq::tax_table(mga$ps)
  taxTab2 <- as.data.frame(taxTab2)

  # Add genus to species name in Species column
  spec.names2 <- taxTab2$Species
  taxTab2$Species <- paste(taxTab2$Genus, taxTab2$Species, sep = " ")

  # Add column to indicate taxa to keep (yes = 1, no = 0)
  taxTab2$keep = 0
  taxTab2$keep_taxon = NA

  # Search for taxa in taxonomy table
  for (j in 1:nrow(taxa)) {
    # Search in Kingdom
    if ("Kingdom" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab2)) {
        if (taxa$taxon[j] %in% taxTab2$Kingdom[i]){
          taxTab2$keep[i] = 1
          taxTab2$keep_taxon[i] = taxa$taxon[j]
        }}
      # Search in Phylum
    } else if ("Phylum" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab2)) {
        if (taxa$taxon[j] %in% taxTab2$Phylum[i]){
          taxTab2$keep[i] = 1
          if (is.na(taxTab2$keep_taxon[i])){
            taxTab2$keep_taxon[i] = taxa$taxon[j]
          }}}
      # Search in Class
    } else if ("Class" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab2)) {
        if (taxa$taxon[j] %in% taxTab2$Class[i]){
          taxTab2$keep[i] = 1
          if (is.na(taxTab2$keep_taxon[i])){
            taxTab2$keep_taxon[i] = taxa$taxon[j]
          }}}
      # Search in Order
    } else if ("Order" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab2)) {
        if (taxa$taxon[j] %in% taxTab2$Order[i]){
          taxTab2$keep[i] = 1
          if (is.na(taxTab2$keep_taxon[i])){
            taxTab2$keep_taxon[i] = taxa$taxon[j]
          }}}
      # Search in Family
    } else if ("Family" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab2)) {
        if (taxa$taxon[j] %in% taxTab2$Family[i]){
          taxTab2$keep[i] = 1
          if (is.na(taxTab2$keep_taxon[i])){
            taxTab2$keep_taxon[i] = taxa$taxon[j]
          }}}
      # Search in Genus
    } else if ("Genus" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab2)) {
        if (taxa$taxon[j] %in% taxTab2$Genus[i]){
          taxTab2$keep[i] = 1
          if (is.na(taxTab2$keep_taxon[i])){
            taxTab2$keep_taxon[i] = taxa$taxon[j]
          }}}
      # Search in Species
    } else if ("Species" %in% taxa$group[j]){
      for (i in 1:nrow(taxTab2)) {
        if (taxa$taxon[j] %in% taxTab2$Species[i]){
          taxTab2$keep[i] = 1
          if (is.na(taxTab2$keep_taxon[i])){
            taxTab2$keep_taxon[i] = taxa$taxon[j]
          }}}
    } else {
      taxTab2$keep[i] = 0
      taxTab2$keep_taxon[i] = NA}
  }

  # Remove genus from species names
  taxTab2$Species <- spec.names2

  # Update the tax_table
  tax2 = phyloseq::tax_table(as.matrix(taxTab2))
  if (keep) {
    phyloseq::tax_table(mga_keep$ps) <- tax2
  }
  if (drop) {
    phyloseq::tax_table(mga_drop$ps) <- tax2
  }

  # Extract selected taxa only
  if (keep) {
    keep.only2 <- taxTab2[taxTab2$keep == 1,]
    keep.names2 <- row.names(keep.only2)
  }
  if (drop) {
    drop.only2 <- taxTab2[taxTab2$keep == 0,]
    drop.names2 <- row.names(drop.only2)
  }

  # Filter taxa from ASV table
  if (keep) {
    mga_keep$ps <- phyloseq::prune_taxa(keep.names2, mga_keep$ps)
  }
  if (drop) {
    mga_drop$ps <- phyloseq::prune_taxa(drop.names2, mga_drop$ps)
  }


  if (keep) {
    if (network) {

      # Construct the network
      mga_keep$net <- phyloseq::make_network(if (group.species){mga_keep$ps.species} else{mga_keep$ps},
                                    type = network.args$type,
                                    distance = network.args$distance,
                                    max.dist = network.args$max.dist,
                                    keep.isolates = network.args$keep.isolates)
      sampledata1 <- as.data.frame(phyloseq::sample_data(mga_keep$ps)) # subset sample_info as a data frame

      # Count vertices and edges in site network
      v1 <- igraph::vcount(mga_keep$net)
      e1 <- igraph::ecount(mga_keep$net)

      connectivity1 <- e1/v1   # avg number associations between ASVs
      connectance1 <- e1/v1^2  # fraction of edges out of all possible edges
      # Degree of each node
      degrees1 <- igraph::degree(mga_keep$net, v = igraph::V(mga_keep$net), mode = "all")

    } else if (network == FALSE) {
      sampledata1 <- as.data.frame(phyloseq::sample_data(mga_keep$ps)) # subset sample_info as a data frame
    }

    # Sum the presences in each sample for species richness
    rich1 <- nrow(phyloseq::tax_table(mga_keep$ps.species)) # species count in sample
    # total individuals from all species in each sample
    n1 <- apply(phyloseq::otu_table(mga_keep$ps), 1, function(x) sum(x))
    ASVs1 <- ncol(phyloseq::otu_table(mga_keep$ps)) # total number of ASVs in site

    # Alpha species diversity measures
    # Shannon index
    ps_alpha_div_1 <- phyloseq::estimate_richness(if (group.species){mga_keep$ps.species} else{mga_keep$ps}, split = TRUE, measure = "Shannon")
    # Simpson index
    ps_alpha_div2_1 <- phyloseq::estimate_richness(if (group.species){mga_keep$ps.species} else{mga_keep$ps}, split = TRUE, measure = "Simpson")

    # Phylogenetic diversity
    # Test if phy_tree is present
    if (nrow(phyloseq::tax_table(if (group.species){mga_keep$ps.species} else{mga_keep$ps})) <= 1) {
      warning("filter_taxa attempted to reduce tree to 1 or fewer tips. tree replaced with NULL. PD will have NA value.")
      PD1 <- NA
    } else {
      # Sum all branch lengths in the phylogenetic tree
      PD1 <- sum(phyloseq::tree_layout(phyloseq::phy_tree(if (group.species){mga_keep$ps.species} else{mga_keep$ps}))$edgeDT[["edge.length"]])
    }

    if (network) {

      # Display sample measures in a data frame
      results.samples1 <- data.frame(ps_alpha_div_1,
                                    ps_alpha_div2_1,
                                    rich = rich1,
                                    PD = PD1,
                                    n.indivs = n1,
                                    ASVs = ASVs1,
                                    vertices = v1,
                                    edges = e1,
                                    connectivity = connectivity1,
                                    connectance = connectance1)

      results.samples1 <- cbind(results.samples1, sampledata1)

      if (group.species) {
        degree.samp1 <- as.data.frame(t(phyloseq::otu_table(mga_keep$ps.species)))
        degree.samp1$degrees <- degrees1
      } else {
        degree.samp1 <- as.data.frame(t(phyloseq::otu_table(mga_keep$ps)))
        degree.samp1$degrees <- degrees1}

    } else if (network == FALSE) {

      # Display sample measures in a data frame
      results.samples1 <- data.frame(ps_alpha_div_1,
                                    ps_alpha_div2_1,
                                    rich = rich1,
                                    PD = PD1,
                                    n.indivs = n1,
                                    ASVs = ASVs1)

      results.samples1 <- cbind(results.samples1, sampledata1)

    }

    # Return filtered mga-class object
    mga_keep$ASVtab <- as.integer(phyloseq::otu_table(mga_keep$ps))
    mga_keep$taxTab <- if (group.species) {keep.only} else {keep.only2}
    mga_keep$results.samples <- results.samples1

    if (network) {
      mga_keep$degree.samp <- degree.samp1
    }

    # Convert to mga-class object
    mga_keep <- new.mga(mga_keep)
  }

  if (drop) {
    if (network) {

      # Construct the network
      mga_drop$net <- phyloseq::make_network(if (group.species){mga_drop$ps.species} else{mga_drop$ps},
                                    type = network.args$type,
                                    distance = network.args$distance,
                                    max.dist = network.args$max.dist,
                                    keep.isolates = network.args$keep.isolates)
      sampledata2 <- as.data.frame(phyloseq::sample_data(mga_drop$ps)) # subset sample_info as a data frame

      # Count vertices and edges in site network
      v2 <- igraph::vcount(mga_drop$net)
      e2 <- igraph::ecount(mga_drop$net)

      connectivity2 <- e2/v2   # avg number associations between ASVs
      connectance2 <- e2/v2^2  # fraction of edges out of all possible edges
      # Degree of each node
      degrees2 <- igraph::degree(mga_drop$net, v = igraph::V(mga_drop$net), mode = "all")

    } else if (network == FALSE) {
      sampledata2 <- as.data.frame(phyloseq::sample_data(mga_drop$ps)) # subset sample_info as a data frame
    }

    # Sum the presences in each sample for species richness
    rich2 <- nrow(phyloseq::tax_table(mga_drop$ps.species)) # species count in sample
    # total individuals from all species in each sample
    n2 <- apply(phyloseq::otu_table(mga_drop$ps), 1, function(x) sum(x))
    ASVs2 <- ncol(phyloseq::otu_table(mga_drop$ps)) # total number of ASVs in site

    # Alpha species diversity measures
    # Shannon index
    ps_alpha_div_2 <- phyloseq::estimate_richness(if (group.species){mga_drop$ps.species} else{mga_drop$ps}, split = TRUE, measure = "Shannon")
    # Simpson index
    ps_alpha_div2_2 <- phyloseq::estimate_richness(if (group.species){mga_drop$ps.species} else{mga_drop$ps}, split = TRUE, measure = "Simpson")

    # Phylogenetic diversity
    if (nrow(phyloseq::tax_table(if (group.species){mga_keep$ps.species} else{mga_keep$ps})) <= 1) {
      warning("filter_taxa attempted to reduce tree to 1 or fewer tips. tree replaced with NULL. PD will have NA value.")
      PD2 <- NA
    } else {
      # Sum all branch lengths in the phylogenetic tree
      PD2 <- sum(phyloseq::tree_layout(phyloseq::phy_tree(if (group.species){mga_keep$ps.species} else{mga_keep$ps}))$edgeDT[["edge.length"]])
    }

    if (network) {

      # Display sample measures in a data frame
      results.samples2 <- data.frame(ps_alpha_div_2,
                                    ps_alpha_div2_2,
                                    rich = rich2,
                                    PD = PD2,
                                    n.indivs = n2,
                                    ASVs = ASVs2,
                                    vertices = v2,
                                    edges = e2,
                                    connectivity = connectivity2,
                                    connectance = connectance2)

      results.samples2 <- cbind(results.samples2, sampledata2)

      if (group.species) {
        degree.samp2 <- as.data.frame(t(phyloseq::otu_table(mga_drop$ps.species)))
        degree.samp2$degrees <- degrees2
      } else {
        degree.samp2 <- as.data.frame(t(phyloseq::otu_table(mga_drop$ps)))
        degree.samp2$degrees <- degrees2}

    } else if (network == FALSE) {

      # Display sample measures in a data frame
      results.samples <- data.frame(ps_alpha_div_2,
                                    ps_alpha_div2_2,
                                    rich = rich2,
                                    PD = PD2,
                                    n.indivs = n2,
                                    ASVs = ASVs2)

      results.samples2 <- cbind(results.samples2, sampledata2)

    }

    # Return filtered mga-class object
    mga_drop$ASVtab <- as.integer(phyloseq::otu_table(mga_drop$ps))
    mga_drop$taxTab <- if (group.species) {drop.only} else {drop.only2}
    mga_drop$results.samples <- results.samples2

    if (network) {
      mga_drop$degree.samp <- degree.samp2
    }

    # Convert to mga-class object
    mga_drop <- new.mga(mga_drop)
  }

  if (keep == TRUE & drop == FALSE) {
    return(mga_keep)
  } else if (keep == FALSE & drop == TRUE) {
    return(mga_drop)
  } else if (keep & drop) {
    return(list(mga_keep, mga_drop))
  }

  } else { # Did NOT identify any selected taxa

    warning("No taxa selected were identified in the mga object. Please check if taxa are present in the taxonomy table.")
    return(mga)
  }
}







# # Taxonomic filtering and pruning
# drop_taxa <- function(mga, # mga-class object
#                       taxa = NULL, # .csv object listing taxa to search for with columns named "taxon" and "group"
#                       group.species = TRUE,
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
#     # Search in Kingdom
#     if ("Kingdom" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Kingdom[i]){
#           taxTab$keep[i] = 1
#           taxTab$keep_taxon[i] = taxa$taxon[j]
#         }}
#       # Search in Phylum
#     } else if ("Phylum" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Phylum[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Class
#     } else if ("Class" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Class[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Order
#     } else if ("Order" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Order[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Family
#     } else if ("Family" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Family[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Genus
#     } else if ("Genus" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Genus[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Species
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
#   mga$ps <- phyloseq::prune_taxa(keep.names, mga$ps)
#
#   # Aggregate duplicate species
#   ps.species <- phyloseq::tax_glom(mga$ps, taxrank = 'Species', NArm = FALSE)
#
#   if (group.species) {
#   mga$ps.species <- ps.species
#   }
#
#   if (network) {
#
#     # Construct the network using `prune.ps`
#     net <- phyloseq::make_network(if (group.species){mga$ps.species} else{mga$ps},
#                                   type = network.args$type,
#                                   distance = network.args$distance,
#                                   max.dist = network.args$max.dist,
#                                   keep.isolates = network.args$keep.isolates)
#     sampledata <- as.data.frame(phyloseq::sample_data(mga$ps)) # subset sample_info as a data frame
#
#     # Count vertices and edges in site network
#     v <- igraph::vcount(net)
#     e <- igraph::ecount(net)
#
#     connectivity <- e/v   # avg number associations between ASVs
#     connectance <- e/v^2  # fraction of edges out of all possible edges
#     # Degree of each node
#     degrees <- igraph::degree(net, v = igraph::V(net), mode = "all")
#
#   } else if (network == FALSE) {
#     sampledata <- as.data.frame(phyloseq::sample_data(mga$ps)) # subset sample_info as a data frame
#   }
#
#   # Sum the presences in each sample for species richness
#   rich <- nrow(phyloseq::tax_table(mga$ps.species)) # species count in sample
#   # total individuals from all species in each sample
#   n <- apply(phyloseq::otu_table(mga$ps), 1, function(x) sum(x))
#   ASVs <- ncol(phyloseq::otu_table(mga$ps)) # total number of ASVs in site
#
#   # Alpha species diversity measures
#   # Shannon index
#   ps_alpha_div <- phyloseq::estimate_richness(if (group.species){mga$ps.species} else{mga$ps}, split = TRUE, measure = "Shannon")
#   # Simpson index
#   ps_alpha_div2 <- phyloseq::estimate_richness(if (group.species){mga$ps.species} else{mga$ps}, split = TRUE, measure = "Simpson")
#
#   # Phylogenetic diversity
#   # Sum all branch lengths in the phylogenetic tree
#   PD <- sum(phyloseq::tree_layout(phyloseq::phy_tree(if (group.species){mga$ps.species} else{mga$ps}))$edgeDT[["edge.length"]])
#
#   if (network) {
#
#     # Display sample measures in a data frame
#     results.samples <- data.frame(ps_alpha_div,
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
#     results.samples <- cbind(results.samples, sampledata)
#
#     if (group.species) {
#       degree.samp <- as.data.frame(t(phyloseq::otu_table(mga$ps.species)))
#       degree.samp$degrees <- degrees
#     } else {
#       degree.samp <- as.data.frame(t(phyloseq::otu_table(mga$ps)))
#       degree.samp$degrees <- degrees}
#
#   } else if (network == FALSE) {
#
#     # Display sample measures in a data frame
#     results.samples <- data.frame(ps_alpha_div,
#                                   ps_alpha_div2,
#                                   rich,
#                                   PD,
#                                   n.indivs = n,
#                                   ASVs)
#
#     results.samples <- cbind(results.samples, sampledata)
#
#   }
#
#   # Return filtered mga-class object
#   mga$ASVtab <- as.integer(phyloseq::otu_table(mga$ps))
#   mga$taxTab <- keep.only
#   mga$results.samples <- results.samples
#
#   if (network) {
#     mga$net <- net
#     mga$degree.samp <- degree.samp
#   }
#
#   # Convert to mga-class object
#   mga <- new.mga(mga)
#
#   return(mga)
# }



# Taxonomic filtering and pruning
# drop_taxa <- function(mga, # mga-class object
#                       taxa = NULL, # .csv object listing taxa to search for with columns named "taxon" and "group"
#                       group.species = TRUE,
#                       network = TRUE,
#                       network.args = list(type = "taxa", # "samples"
#                                           distance = "jaccard", # "jaccard", "unifrac", "wunifrac", DPCoA", "jsd"
#                                           max.dist = 0.45, # default distance set by `make_network` function from `phyloseq` package is 0.4
#                                           keep.isolates = TRUE)) {
#
#   # Convert mga ps,species tax_table to a data frame
#   taxTab <- phyloseq::tax_table(mga$ps.species)
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
#     # Search in Kingdom
#     if ("Kingdom" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Kingdom[i]){
#           taxTab$keep[i] = 1
#           taxTab$keep_taxon[i] = taxa$taxon[j]
#         }}
#       # Search in Phylum
#     } else if ("Phylum" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Phylum[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Class
#     } else if ("Class" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Class[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Order
#     } else if ("Order" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Order[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Family
#     } else if ("Family" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Family[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Genus
#     } else if ("Genus" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab)) {
#         if (taxa$taxon[j] %in% taxTab$Genus[i]){
#           taxTab$keep[i] = 1
#           if (is.na(taxTab$keep_taxon[i])){
#             taxTab$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Species
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
#   phyloseq::tax_table(mga$ps.species) <- tax
#
#   # Extract selected taxa only
#   keep.only <- taxTab[taxTab$keep == 1,]
#   keep.names <- row.names(keep.only)
#
#   # Filter taxa from ASV table
#   mga$ps.species <- phyloseq::prune_taxa(keep.names, mga$ps.species)
#
#
#   # Convert mga ps tax_table to a data frame
#   taxTab2 <- phyloseq::tax_table(mga$ps)
#   taxTab2 <- as.data.frame(taxTab2)
#
#   # Add genus to species name in Species column
#   spec.names2 <- taxTab2$Species
#   taxTab2$Species <- paste(taxTab2$Genus, taxTab2$Species, sep = " ")
#
#   # Add column to indicate taxa to keep (yes = 1, no = 0)
#   taxTab2$keep = 0
#   taxTab2$keep_taxon = NA
#
#   # Search for taxa in taxonomy table
#   for (j in 1:nrow(taxa)) {
#     # Search in Kingdom
#     if ("Kingdom" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab2)) {
#         if (taxa$taxon[j] %in% taxTab2$Kingdom[i]){
#           taxTab2$keep[i] = 1
#           taxTab2$keep_taxon[i] = taxa$taxon[j]
#         }}
#       # Search in Phylum
#     } else if ("Phylum" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab2)) {
#         if (taxa$taxon[j] %in% taxTab2$Phylum[i]){
#           taxTab2$keep[i] = 1
#           if (is.na(taxTab2$keep_taxon[i])){
#             taxTab2$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Class
#     } else if ("Class" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab2)) {
#         if (taxa$taxon[j] %in% taxTab2$Class[i]){
#           taxTab2$keep[i] = 1
#           if (is.na(taxTab2$keep_taxon[i])){
#             taxTab2$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Order
#     } else if ("Order" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab2)) {
#         if (taxa$taxon[j] %in% taxTab2$Order[i]){
#           taxTab2$keep[i] = 1
#           if (is.na(taxTab2$keep_taxon[i])){
#             taxTab2$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Family
#     } else if ("Family" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab2)) {
#         if (taxa$taxon[j] %in% taxTab2$Family[i]){
#           taxTab2$keep[i] = 1
#           if (is.na(taxTab2$keep_taxon[i])){
#             taxTab2$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Genus
#     } else if ("Genus" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab2)) {
#         if (taxa$taxon[j] %in% taxTab2$Genus[i]){
#           taxTab2$keep[i] = 1
#           if (is.na(taxTab2$keep_taxon[i])){
#             taxTab2$keep_taxon[i] = taxa$taxon[j]
#           }}}
#       # Search in Species
#     } else if ("Species" %in% taxa$group[j]){
#       for (i in 1:nrow(taxTab2)) {
#         if (taxa$taxon[j] %in% taxTab2$Species[i]){
#           taxTab2$keep[i] = 1
#           if (is.na(taxTab2$keep_taxon[i])){
#             taxTab2$keep_taxon[i] = taxa$taxon[j]
#           }}}
#     } else {
#       taxTab2$keep[i] = 0
#       taxTab2$keep_taxon[i] = NA}
#   }
#
#   # Remove genus from species names
#   taxTab2$Species <- spec.names2
#
#   # Update the tax_table
#   tax = phyloseq::tax_table(as.matrix(taxTab2))
#   phyloseq::tax_table(mga$ps) <- tax
#
#   # Extract selected taxa only
#   keep.only <- taxTab2[taxTab2$keep == 1,]
#   keep.names <- row.names(keep.only)
#
#   # Filter taxa from ASV table
#   mga$ps <- phyloseq::prune_taxa(keep.names, mga$ps)
#
#
#   if (network) {
#
#     # Construct the network
#     net <- phyloseq::make_network(if (group.species){mga$ps.species} else{mga$ps},
#                                   type = network.args$type,
#                                   distance = network.args$distance,
#                                   max.dist = network.args$max.dist,
#                                   keep.isolates = network.args$keep.isolates)
#     sampledata <- as.data.frame(phyloseq::sample_data(mga$ps)) # subset sample_info as a data frame
#
#     # Count vertices and edges in site network
#     v <- igraph::vcount(net)
#     e <- igraph::ecount(net)
#
#     connectivity <- e/v   # avg number associations between ASVs
#     connectance <- e/v^2  # fraction of edges out of all possible edges
#     # Degree of each node
#     degrees <- igraph::degree(net, v = igraph::V(net), mode = "all")
#
#   } else if (network == FALSE) {
#     sampledata <- as.data.frame(phyloseq::sample_data(mga$ps)) # subset sample_info as a data frame
#   }
#
#   # Sum the presences in each sample for species richness
#   rich <- nrow(phyloseq::tax_table(mga$ps.species)) # species count in sample
#   # total individuals from all species in each sample
#   n <- apply(phyloseq::otu_table(mga$ps), 1, function(x) sum(x))
#   ASVs <- ncol(phyloseq::otu_table(mga$ps)) # total number of ASVs in site
#
#   # Alpha species diversity measures
#   # Shannon index
#   ps_alpha_div <- phyloseq::estimate_richness(if (group.species){mga$ps.species} else{mga$ps}, split = TRUE, measure = "Shannon")
#   # Simpson index
#   ps_alpha_div2 <- phyloseq::estimate_richness(if (group.species){mga$ps.species} else{mga$ps}, split = TRUE, measure = "Simpson")
#
#   # Phylogenetic diversity
#   # Sum all branch lengths in the phylogenetic tree
#   PD <- sum(phyloseq::tree_layout(phyloseq::phy_tree(if (group.species){mga$ps.species} else{mga$ps}))$edgeDT[["edge.length"]])
#
#   if (network) {
#
#     # Display sample measures in a data frame
#     results.samples <- data.frame(ps_alpha_div,
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
#     results.samples <- cbind(results.samples, sampledata)
#
#     if (group.species) {
#       degree.samp <- as.data.frame(t(phyloseq::otu_table(mga$ps.species)))
#       degree.samp$degrees <- degrees
#     } else {
#       degree.samp <- as.data.frame(t(phyloseq::otu_table(mga$ps)))
#       degree.samp$degrees <- degrees}
#
#   } else if (network == FALSE) {
#
#     # Display sample measures in a data frame
#     results.samples <- data.frame(ps_alpha_div,
#                                   ps_alpha_div2,
#                                   rich,
#                                   PD,
#                                   n.indivs = n,
#                                   ASVs)
#
#     results.samples <- cbind(results.samples, sampledata)
#
#   }
#
#   # Return filtered mga-class object
#   mga$ASVtab <- as.integer(phyloseq::otu_table(mga$ps))
#   mga$taxTab <- keep.only
#   mga$results.samples <- results.samples
#
#   if (network) {
#     mga$net <- net
#     mga$degree.samp <- degree.samp
#   }
#
#   # Convert to mga-class object
#   mga <- new.mga(mga)
#
#   return(mga)
# }
