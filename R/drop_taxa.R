#' @export
# Taxonomic filtering and pruning
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

    if (group.species) {
      degree.samp <- as.data.frame(t(phyloseq::otu_table(ps.species)))
      degree.samp$degrees <- degrees
    } else {
      degree.samp <- as.data.frame(t(phyloseq::otu_table(mga$ps)))
      degree.samp$degrees <- degrees}

  } else if (network == FALSE) {

    # Display sample measures in a data frame
    results.samples <- data.frame(ps_alpha_div,
                                  ps_alpha_div2,
                                  rich,
                                  PD,
                                  n.indivs = n,
                                  ASVs)

    results.samples <- cbind(results.samples, sampledata)

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
