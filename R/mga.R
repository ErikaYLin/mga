#' @export
# Microbial Genetic Analysis with co-occurrence network construction
mga <- function(fastq.Fs, fastq.Rs, # file paths for forward and reverse raw fastq files
                filtFs, filtRs, # file paths for filtered and trimmed sequences
                refFasta, # file path of reference FASTA for taxonomic classification
                metadata = NULL, # sample metadata in list or data frame format
                make.unique = TRUE, # defines all unclassified taxa as being unique
                group.species = TRUE, # agglomerates phyloseq taxa by species for computed metrics
                tree.args = list(k = 4,
                                 inv = 0.2,
                                 model = "GTR", # "JC", "F81", "SYM", "GTR", etc. **SEE phangorn::optim.pml FOR ALL MODELS**
                                 rearrangement = "stochastic"), # "none", "NNI", "stochastic", "ratchet"
                network = TRUE, # runs network analysis
                network.args = list(type = "taxa", # "samples"
                                    distance = "jaccard", # "jaccard", "unifrac", "wunifrac", DPCoA", "jsd"
                                    max.dist = 0.45, # default distance set by `make_network` function from `phyloseq` package is 0.4
                                    keep.isolates = TRUE),
                seed = 100,
                multithread = FALSE,
                verbose = TRUE) {


  set.seed(seed)

  # Extract sample names, assuming file paths have format: FOL_DER/SAMPLENAME_XXX.fastq.gz OR FOL_DER/SAMPLENAME_XXX.fastq
  sampleNames <- sapply(strsplit(fastq.Fs, "_"), `[`, 2)
  sampleNames <- sapply(strsplit(sampleNames, "/"), `[`, 2)

  if (verbose) {message("Dereplicating files.")}

  # Dereplicating identical sequences
  derepFs <- dada2::derepFastq(filtFs, verbose = TRUE)
  derepRs <- dada2::derepFastq(filtRs, verbose = TRUE)
  # Naming the derep-class objects by the sample names
  names(derepFs) <- sampleNames
  names(derepRs) <- sampleNames

  if (verbose) {message("Learning sequencing errors pre-DADA2.")}

  # Estimating sequencing error rates pre-DADA2
  # The application of the DADA2 method follows the workflow by Callahan et al. (2017).
  # https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
  errF <- dada2::learnErrors(filtFs, multithread = multithread)
  errR <- dada2::learnErrors(filtRs, multithread = multithread)

  if (verbose) {message("DADA2 to infer sequences.")}

  # DADA2 method to infer ASVs
  dadaFs <- dada2::dada(derepFs, err = errF, multithread = multithread)
  dadaRs <- dada2::dada(derepRs, err = errR, multithread = multithread)

  if (verbose) {message("Merging forward and reverse sequence pairs. Constructing the sequence table.")}

  # Merging sequence pairs into object "merged"
  merged <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs)
  # Constructing sequence table
  seqtabAll <- dada2::makeSequenceTable(merged)
  rownames(seqtabAll) <- sampleNames

  if (verbose) {message("Removing chimeras to build the ASV table.")}

  # Remove chimeras from the sequence table
  ASVtab <- dada2::removeBimeraDenovo(seqtabAll)

  if (verbose) {message("Taxonomic classification based on chosen reference FASTA.")}

  # Assign taxonomy to the ASVs based on chosen reference FASTA and build a taxonomy table
  taxTab2 <- dada2::assignTaxonomy(ASVtab,
                                   refFasta = refFasta,
                                   multithread = multithread)

  if (verbose) {message("Labelling unclassified taxa.")}

  # Remove the prefix from taxa names
  # Removal of prefixes follow the workflow by Hui (2021).
  # https://www.yanh.org/2021/01/01/microbiome-r/#build-phyloseq-project
  taxTab <- data.frame(row.names = row.names(taxTab2),
                       Kingdom = stringr::str_replace(taxTab2[,1], "k__",""),
                       Phylum = stringr::str_replace(taxTab2[,2], "p__",""),
                       Class = stringr::str_replace(taxTab2[,3], "c__",""),
                       Order = stringr::str_replace(taxTab2[,4], "o__",""),
                       Family = stringr::str_replace(taxTab2[,5], "f__",""),
                       Genus = stringr::str_replace(taxTab2[,6], "g__",""),
                       if (ncol(taxTab2) == 7) {
                       Species = stringr::str_replace(taxTab2[,7], "s__","")})

  # Rename "NA" elements to "Unclassified __"
  if (ncol(taxTab) == 7) {
    for (i in 1:nrow(taxTab)) {
      if (is.na(taxTab[i, 2]) || taxTab[i, 2] == "") {
        kingdom <- paste("Unclassified", taxTab[i, 1], sep = "_")
        taxTab[i, 2:7] <- kingdom
      } else if (is.na(taxTab[i, 3]) || taxTab[i, 3] == "") {
        phylum <- paste("Unclassified", taxTab[i, 2], sep = "_")
        taxTab[i, 3:7] <- phylum
      } else if (is.na(taxTab[i, 4]) || taxTab[i, 4] == "") {
        class <- paste("Unclassified", taxTab[i, 3], sep = "_")
        taxTab[i, 4:7] <- class
      } else if (is.na(taxTab[i, 5]) || taxTab[i, 5] == "") {
        order <- paste("Unclassified", taxTab[i, 4], sep = "_")
        taxTab[i, 5:7] <- order
      } else if (is.na(taxTab[i, 6]) || taxTab[i, 6] == "") {
        family <- paste("Unclassified", taxTab[i, 5], sep = "_")
        taxTab[i, 6:7] <- family
      } else if (is.na(taxTab[i, 7]) || taxTab[i, 7] == "") {
        taxTab$Species[i] <- paste("Unclassified", taxTab$Genus[i], sep = "_")
      }
    }
  } else if (ncol(taxTab) == 6) {
    for (i in 1:nrow(taxTab)) {
      if (is.na(taxTab[i, 2]) || taxTab[i, 2] == "") {
        kingdom <- paste("Unclassified", taxTab[i, 1], sep = "_")
        taxTab[i, 2:6] <- kingdom
      } else if (is.na(taxTab[i, 3]) || taxTab[i, 3] == "") {
        phylum <- paste("Unclassified", taxTab[i, 2], sep = "_")
        taxTab[i, 3:6] <- phylum
      } else if (is.na(taxTab[i, 4]) || taxTab[i, 4] == "") {
        class <- paste("Unclassified", taxTab[i, 3], sep = "_")
        taxTab[i, 4:6] <- class
      } else if (is.na(taxTab[i, 5]) || taxTab[i, 5] == "") {
        order <- paste("Unclassified", taxTab[i, 4], sep = "_")
        taxTab[i, 5:6] <- order
      } else if (is.na(taxTab[i, 6]) || taxTab[i, 6] == "") {
        family <- paste("Unclassified", taxTab[i, 5], sep = "_")
        taxTab[i, 6:6] <- family
      }
    }
  }

  # Add species column if missing
  if (ncol(taxTab) == 6) {
    taxTab$Species = paste(taxTab$Genus, "sp", sep = "_")
  }

  # Number replicates of unclassified taxa
  if (make.unique){
    for(i in 1: ncol(taxTab)){
      # finding all unique unclassified taxa per column (taxon) and assigning them to UNIQUE
      UNIQUE <- unique(gsub("Unclassified_", "",
                            taxTab[,i][grepl("Unclassified_", taxTab[,i])]))
      # nested for loop that repeats the assessment within the same column for index j
      for(j in 1:length(UNIQUE)){
        # criterion to search for the unclassified taxa and make them unique
        criterion <- paste("Unclassified_", UNIQUE[j], sep = "")
        taxTab[,i][taxTab[,i] == criterion] <- make.unique(taxTab[,i][taxTab[,i] == criterion])
      }
    }
  }

  if (verbose) {message("Multiple sequence alignment for the phylogeny.")}

  # Multiple sequence alignment pre-tree
  # Phylogenetic tree construction follows the workflow by Callahan et al. (2017).
  # https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
  seqs <- dada2::getSequences(ASVtab)
  # names(seqs) <- seqs # This propagates to the tip labels of the tree # try replacing with: paste(taxTab$Genus, taxTab$Species, sep = "_")
  names(seqs) <- paste(taxTab$Genus, taxTab$Species, sep = "_")
  alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor=NA, verbose=FALSE)

  if (verbose) {message("Fitting the phylogenetic tree.")}

  # Phylogenetic tree construction
  phangAlign <- phangorn::phyDat(methods::as(alignment, "matrix"), type = "DNA")
  dm <- phangorn::dist.ml(phangAlign)
  treeNJ <- phangorn::NJ(dm) # Note: tip order != sequence order
  fit = phangorn::pml(treeNJ, data = phangAlign)
  fitGTR <- stats::update(fit, k = tree.args$k, inv = tree.args$inv) # consider sensitivity analysis on these values
  fitGTR <- phangorn::optim.pml(fitGTR, model = tree.args$model, optInv = TRUE, optGamma = TRUE,
                                rearrangement = tree.args$rearrangement,
                                control = phangorn::pml.control(trace = 0))

  if (verbose) {message("Storing phylogenetic outputs in a list.")}

  phylo_tree <- list()
  phylo_tree$phangAlign <- phangAlign
  phylo_tree$dm <- dm
  phylo_tree$treeNJ <- treeNJ
  phylo_tree$fit <- fit
  phylo_tree$fitGTR <- fitGTR

  if (verbose) {message("Checking and reformatting the sample metadata.")}

  # Use metadata if provided, sample names only if metadata not available
  if (!is.null(metadata)){
    sample_info <- as.data.frame(metadata)
    # Use existing sample.ID if available, assign sample.ID if not
    if (sampleNames %in% metadata$sample.ID){sample_info$sample.ID = metadata$sample.ID}
    else {sample_info$sample.ID <- sampleNames}
  } else {sample_info <- data.frame(sample.ID = sampleNames)}

  # Check that all samples from the ASV table are found in the sample metadata
  all(rownames(seqtabAll) %in% sample_info$sample.ID) # TRUE
  # Set row names as sample names for the tabulated metadata file (originally in `.csv` format)
  rownames(sample_info) <- sample_info$sample.ID

  if (verbose) {message("Combining data into a phyloseq object.")}

  # Now all the sequence and taxonomic data can be combined into a single phyloseq-class object
  tax = phyloseq::tax_table(as.matrix(taxTab))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(ASVtab,
                                               taxa_are_rows = FALSE),
                           phyloseq::sample_data(sample_info),
                           tax,
                           phyloseq::phy_tree(fitGTR$tree))

  # Aggregate duplicate species
  species <- phyloseq::tax_glom(ps, taxrank = 'Species', NArm = FALSE)

  if (network) {

    if (verbose) {message("Building the network from the phyloseq object.")}

    # Construct the network using `ps`
    net <- phyloseq::make_network(if (group.species){species} else{ps},
                                  type = network.args$type,
                                  distance = network.args$distance,
                                  max.dist = network.args$max.dist,
                                  keep.isolates = network.args$keep.isolates)
    sampledata <- as.data.frame(phyloseq::sample_data(ps)) # subset sample_info as a data frame

    if (verbose) {message("Computing site-wide network diagnostics.")}

    # Count vertices and edges in site network
    v <- igraph::vcount(net)
    e <- igraph::ecount(net)

    connectivity <- e/v   # avg number associations between ASVs
    connectance <- e/v^2  # fraction of edges out of all possible edges
    # Degree of each node
    degrees <- igraph::degree(net, v = igraph::V(net), mode = "all")

  } else if (network == FALSE) {
    sampledata <- as.data.frame(phyloseq::sample_data(ps)) # subset sample_info as a data frame
  }

  if (verbose) {message("Calculating diversity measures for each sample.")}

  # Sum the presences in each sample for species richness
  rich <- nrow(phyloseq::tax_table(species)) # species count in sample
  # total individuals from all species in each sample
  n <- apply(ASVtab, 1, function(x) sum(x))
  ASVs <- ncol(ASVtab) # total number of ASVs in site

  # Alpha species diversity measures
  # Shannon index
  ps_alpha_div <- phyloseq::estimate_richness(if (group.species){species} else{ps}, split = TRUE, measure = "Shannon")
  # Simpson index
  ps_alpha_div2 <- phyloseq::estimate_richness(if (group.species){species} else{ps}, split = TRUE, measure = "Simpson")

  # Phylogenetic diversity
  # Sum all branch lengths in the phylogenetic tree
  PD <- sum(phyloseq::tree_layout(phyloseq::phy_tree(if (group.species){species} else{ps}))$edgeDT[["edge.length"]])

  if (network) {

    if (verbose) {message("Storing diversity measures and network diagnostics in data frames.")}

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

    results.samples <- cbind(results.samples, sample_info)
    # results.samples <- dplyr::relocate(sample.ID, .before = Shannon) # move Sample.ID column to leftmost

    degree.samp <- data.frame(ASVtab[1,], degrees, row.names = colnames(ASVtab))
    colnames(degree.samp)[1] <- sampledata$sample.ID

  } else if (network == FALSE) {

    if (verbose) {message("Storing diversity measures in a data frame.")}

    # Display sample measures in a data frame
    results.samples <- data.frame(ps_alpha_div,
                                  ps_alpha_div2,
                                  rich,
                                  PD,
                                  n.indivs = n,
                                  ASVs)

    results.samples <- cbind(results.samples, sample_info)
    # results.samples <- dplyr::relocate(sample.ID, .before = Shannon) # move Sample.ID column to leftmost

  }

  if (verbose) {message("Building mga list of outputs.")}

  # list of all objects to return
  ps_network <- list()
  ps_network$seqtabAll <- seqtabAll
  ps_network$ASVtab <- ASVtab
  ps_network$taxTab <- taxTab
  ps_network$phylo_tree <- phylo_tree
  ps_network$sampledata <- sampledata
  ps_network$ps <- ps

  if (group.species) {
  ps_network$ps.species <- species
  }

  ps_network$results.samples <- results.samples

  if (network) {
    ps_network$net <- net
    ps_network$degree.samp <- degree.samp
  }

  # Convert to `mga` class
  ps_network <- new.mga(ps_network)

  return(ps_network)
}




