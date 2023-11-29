#' @export
filter_pathogens <- function(x, # phyloseq object
                             pathogens = NULL # .csv object listing taxa to search for
                             ) {

  # Convert phyloseq tax_table to a data frame
  pathogen <- phyloseq::tax_table(ps.species)
  pathogen <- as.data.frame(pathogen)

  # Add genus to species name in Species column
  spec.names <- pathogen$Species
  pathogen$Species <- paste(pathogen$Genus, pathogen$Species, sep = " ")

  # Add column to indicate possible pathogenicity (yes = 1, no = 0)
  pathogen$pathogenic = 0
  pathogen$pathogenic_taxon = NA

  # Search for pathogens in taxonomy table
  for (j in 1:nrow(pathogens)) {
    # Search in Family
    if ("Family" %in% pathogens$Taxon[j]){
      for (i in 1:nrow(pathogen)) {
        if (pathogens$Pathogen[j] %in% pathogen$Family[i]){
          pathogen$pathogenic[i] = 1
          pathogen$pathogenic_taxon[i] = pathogens$Pathogen[j]
        }}
      # Search in Genus
    } else if ("Genus" %in% pathogens$Taxon[j]){
      for (i in 1:nrow(pathogen)) {
        if (pathogens$Pathogen[j] %in% pathogen$Genus[i]){
          pathogen$pathogenic[i] = 1
          if (is.na(pathogen$pathogenic_taxon[i])){
            pathogen$pathogenic_taxon[i] = pathogens$Pathogen[j]
          }}}
      # search in Species
    } else if ("Species" %in% pathogens$Taxon[j]){
      for (i in 1:nrow(pathogen)) {
        if (pathogens$Pathogen[j] %in% pathogen$Species[i]){
          pathogen$pathogenic[i] = 1
          if (is.na(pathogen$pathogenic_taxon[i])){
            pathogen$pathogenic_taxon[i] = pathogens$Pathogen[j]
          }}}
    } else {
      pathogen$pathogenic[i] = 0
      pathogen$pathogenic_taxon[i] = NA}
  }

  # Remove genus from species names
  pathogen$Species <- spec.names

  # Update the tax_table
  tax = phyloseq::tax_table(as.matrix(pathogen))
  phyloseq::tax_table(ps.species) <- tax

  # Extract pathogens only
  pathogen.only <- pathogen[pathogen$pathogenic == 1,]
  path.names <- row.names(pathogen.only)

  # Filter pathogens from ASV table
  pathogen.ps <- phyloseq::prune_taxa(path.names, ps.species)

  # Return the filtered phyloseq object
  return(pathogen.ps)
}
