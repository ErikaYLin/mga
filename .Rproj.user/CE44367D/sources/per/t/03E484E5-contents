# FUNCTION: Dereplicates, infers ASVs (DADA2 method), merges sequence pairs, builds sequence table, removes chimeras, assigns taxonomy, generates a phylogenetic tree, combines all data into a phyloseq-class object, constructs a network plot, extracts measures of community structure and network topology.
# Returns: List with sequence table, table of ASVs, taxonomy table, phylogenetic tree, sample metadata, phyloseq object, network plot, data frame of richness/diversity metrics and network diagnostics for each sample.
# NOTE: "fastq.fwd" and "fastq.rev" are the specific file paths for the forward and reverse reads of each sample.

ps.net <- function(fastq.Fs, fastq.Rs, # file paths for forward and reverse raw fastq files
                   filtFs, filtRs, # file paths for filtered and trimmed sequences
                   refFasta, # file path of reference FASTA for taxonomic classification
                   metadata = NULL, # sample metadata in list or data frame format
                   # filter.args = list(truncLen = c(240,160),
                   #                    maxN = 0,
                   #                    maxEE = c(2,2),
                   #                    truncQ = 2,
                   #                    rm.phix = TRUE, # defining max 2 expected errors per read
                   #                    compress = TRUE,
                   #                    multithread = FALSE), # On Windows set multithread = FALSE
                   tree.args = list(k = 4, 
                                    inv = 0.2,
                                    model = "GTR", # "JC", "F81", "SYM", "GTR", etc. **SEE phangorn::optim.pml FOR ALL MODELS**
                                    rearrangement = "stochastic"), # "none", "NNI", "stochastic", "ratchet"
                   network.args = list(type = "taxa", # "samples"
                                       distance = "jaccard", # "jaccard", "unifrac", "wunifrac", DPCoA", "jsd"
                                       max.dist = 0.45, # default distance set by `make_network` function from `phyloseq` package is 0.4
                                       keep.isolates = TRUE), 
                   seed = 100,
                   multithread = TRUE) {
  
  
  
  # ARGUMENTS FOR DEBUG
  
  # # Load sample metadata
  # samdf <- read.csv("C:/Users/elinrg26/OneDrive - UBC/Research/GitHub/Cover Crop Study/Fungal_Community_Workflow/Mapping file for ITS sequencing.csv")
  # # Define filter path for the forward and reverse reads of each sample fastq file
  # seq_path <- "./CC_Seq"  # Replace with the directory for fastq file after unzipping the folder
  # 
  # # Sort forward and reverse reads to be in the same order
  # fnFs <- sort(list.files(seq_path, pattern="_R1_001.fastq.gz"))
  # fnRs <- sort(list.files(seq_path, pattern="_R2_001.fastq.gz"))
  # # Specify the full path to the fnFs and fnRs
  # fnFs <- file.path(seq_path, fnFs)
  # fnRs <- file.path(seq_path, fnRs)
  # 
  # # Extract sample names, assuming file paths have format: FOL_DER/SAMPLENAME_XXX.fastq.gz OR FOL_DER/SAMPLENAME_XXX.fastq
  # sNames <- sapply(strsplit(fnFs, "_"), `[`, 2)
  # sNames <- sapply(strsplit(sNames, "/"), `[`, 2)
  # 
  # # Define file names for filtered fastq.gz files and assign file paths
  # filt_path <- file.path(seq_path, "filtered") # Place filtered files in new subdirectory if one does not already exist
  # if(!file_test("-d", filt_path)) dir.create(filt_path)
  # filtFs <- file.path(filt_path, paste0(sNames, "_F_filt.fastq.gz"))
  # filtRs <- file.path(filt_path, paste0(sNames, "_R_filt.fastq.gz"))
  # 
  # fastq.Fs = fnFs
  # fastq.Rs = fnRs
  # 
  # filter.args = list(truncLen = c(240,160),
  #                    maxN = 0,
  #                    maxEE = c(2,2),
  #                    truncQ = 2,
  #                    rm.phix = TRUE, # defining max 2 expected errors per read
  #                    compress = TRUE,
  #                    multithread = FALSE) # On Windows set multithread = FALSE
  # fungi <- "./sh_general_release_dynamic_29.11.2022.fasta"
  # refFasta = fungi # file path of reference FASTA for taxonomic classification
  # # metadata argument
  # for (i in 1:length(fnFs)){
  # 
  #   meta_data <- list()
  #   for (j in 1:length(fnFs)){
  #     ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
  #     ID <- sapply(strsplit(ID, "/"), `[`, 2)
  #     meta_data[[j]] <- as.list(samdf[samdf$sample.ID %in% ID,])
  #   }
  # }
  # metadata = meta_data[[j]] # samdf (sample metadata in list or data frame format)
  # tree.args = list(k = 4,
  #                  inv = 0.2,
  #                  model = "GTR", # "Jukes-Cantor", "F81", "symmetric", "GTR"
  #                  rearrangement = "stochastic") # "none", "NNI", "stochastic", "ratchet"
  # network.args = list(type = "taxa", # "samples"
  #                     distance = "jaccard", # "jaccard", "unifrac", "wunifrac", DPCoA", "jsd"
  #                     max.dist = 0.45,
  #                     keep.isolates = TRUE) # default distance set by `make_network` function from `phyloseq` package
  # seed = 100
  # multithread = TRUE

  
  
  
  
  
  
  
  set.seed(seed)
  
  # Extract sample names, assuming file paths have format: FOL_DER/SAMPLENAME_XXX.fastq.gz OR FOL_DER/SAMPLENAME_XXX.fastq
  sampleNames <- sapply(strsplit(fastq.Fs, "_"), `[`, 2)
  sampleNames <- sapply(strsplit(sampleNames, "/"), `[`, 2)
  
  # message("Filtering and trimming files.")
  
  # Filter forward and reverse reads, trim sequences to positions 240 and 160.
  # Filtering, trimming and dereplication of sequences follow the workflow by Callahan et al. (2017).
  # https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
  # out <- dada2::filterAndTrim(fastq.Fs, filtFs, fastq.Rs, filtRs,
  #                             truncLen = filter.args$truncLen,
  #                             maxN = filter.args$maxN,
  #                             maxEE = filter.args$maxEE,
  #                             truncQ = filter.args$truncQ,
  #                             rm.phix = filter.args$rm.phix,
  #                             compress = filter.args$compress,
  #                             multithread = filter.args$multithread) 
  
  message("Dereplicating files.")
  
  # Dereplicating identical sequences
  derepFs <- dada2::derepFastq(filtFs, verbose = TRUE)
  derepRs <- dada2::derepFastq(filtRs, verbose = TRUE)
  # Naming the derep-class objects by the sample names
  names(derepFs) <- sampleNames
  names(derepRs) <- sampleNames
  
  message("Learning sequencing errors pre-DADA2.")
  
  # Estimating sequencing error rates pre-DADA2
  # The application of the DADA2 method follows the workflow by Callahan et al. (2017).
  # https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
  errF <- dada2::learnErrors(filtFs, multithread = multithread)
  errR <- dada2::learnErrors(filtRs, multithread = multithread)
  
  message("DADA2 to infer sequences.")
  
  # DADA2 method to infer ASVs
  dadaFs <- dada2::dada(derepFs, err = errF, multithread = multithread)
  dadaRs <- dada2::dada(derepRs, err = errR, multithread = multithread)
  
  message("Merging forward and reverse sequence pairs. Constructing the sequence table.")
  
  # Merging sequence pairs into object "merged"
  merged <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs)
  # Constructing sequence table
  seqtabAll <- dada2::makeSequenceTable(merged)
  rownames(seqtabAll) <- sampleNames
  
  message("Removing chimeras to build the ASV table.")
  
  # Remove chimeras from the sequence table
  ASVtab <- dada2::removeBimeraDenovo(seqtabAll)
  
  message("Taxonomic classification based on chosen reference FASTA.")
  
  # Assign taxonomy to the ASVs based on chosen reference FASTA and build a taxonomy table
  taxTab <- dada2::assignTaxonomy(ASVtab,
                                  refFasta = refFasta,
                                  multithread = multithread)
  
  message("Labelling unclassified taxa.")
  
  # Remove the prefix from taxa names
  # Removal of prefixes follow the workflow by Hui (2021).
  # https://www.yanh.org/2021/01/01/microbiome-r/#build-phyloseq-project
  taxTab <- data.frame(row.names = row.names(taxTab),
                       Kingdom = stringr::str_replace(taxTab[,1], "k__",""),
                       Phylum = stringr::str_replace(taxTab[,2], "p__",""),
                       Class = stringr::str_replace(taxTab[,3], "c__",""),
                       Order = stringr::str_replace(taxTab[,4], "o__",""),
                       Family = stringr::str_replace(taxTab[,5], "f__",""),
                       Genus = stringr::str_replace(taxTab[,6], "g__",""),
                       Species = stringr::str_replace(taxTab[,7], "s__",""))
  
  # Rename "NA" elements to "Unclassified __"
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
  
  # Number replicates of unclassified taxa
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
  
  message("Multiple sequence alignment for the phylogeny.")
  
  # Multiple sequence alignment pre-tree
  # Phylogenetic tree construction follows the workflow by Callahan et al. (2017).
  # https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
  seqs <- dada2::getSequences(ASVtab)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor=NA, verbose=FALSE)
  
  message("Fitting the phylogenetic tree.")
  
  # Phylogenetic tree construction
  phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- phangorn::dist.ml(phangAlign)
  treeNJ <- phangorn::NJ(dm) # Note: tip order != sequence order
  fit = phangorn::pml(treeNJ, data = phangAlign)
  fitGTR <- stats::update(fit, k = tree.args$k, inv = tree.args$inv) # consider sensitivity analysis on these values
  fitGTR <- phangorn::optim.pml(fitGTR, model = tree.args$model, optInv = TRUE, optGamma = TRUE,
                                rearrangement = tree.args$rearrangement,
                                control = phangorn::pml.control(trace = 0))
  
  message("Storing phylogenetic outputs in a list.")
  
  phylo_tree <- list()
  phylo_tree$phangAlign <- phangAlign
  phylo_tree$dm <- dm
  phylo_tree$treeNJ <- treeNJ
  phylo_tree$fit <- fit
  phylo_tree$fitGTR <- fitGTR
  
  message("Checking and reformatting the sample metadata.")
  
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
  
  message("Combining data into a phyloseq object.")
  
  # Now all the sequence and taxonomic data can be combined into a single phyloseq-class object
  tax = phyloseq::tax_table(as.matrix(taxTab))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(ASVtab,
                                               taxa_are_rows = FALSE),
                           phyloseq::sample_data(sample_info),
                           tax,
                           phyloseq::phy_tree(fitGTR$tree))
  
  message("Building the network from the phyloseq object.")
  
  # Construct the network using `ps`
  net <- phyloseq::make_network(ps,
                                type = network.args$type,
                                distance = network.args$distance,
                                max.dist = network.args$max.dist,
                                keep.isolates = network.args$keep.isolates)
  sampledata <- as.data.frame(phyloseq::sample_data(ps)) # subset sample_info as a data frame
  
  message("Computing site-wide network diagnostics.")
  
  # Count vertices and edges in site network
  v <- igraph::vcount(net)
  e <- igraph::ecount(net)
  
  connectivity <- e/v   # avg number associations between ASVs
  connectance <- e/v^2  # fraction of edges out of all possible edges
  # Degree of each node
  degrees <- igraph::degree(net, v = igraph::V(net), mode = "all")
  
  message("Calculating diversity measures for each sample.")
  
  # Aggregate duplicate species
  species <- phyloseq::tax_glom(ps, taxrank = 'Species', NArm = FALSE)
  
  # Sum the presences in each sample for species richness
  rich <- nrow(phyloseq::tax_table(species)) # species count in sample
  # total individuals from all species in each sample
  n <- apply(ASVtab, 1, function(x) sum(x))
  ASVs <- ncol(ASVtab) # total number of ASVs in site
  
  # Alpha species diversity measures
  # Shannon index
  ps_alpha_div <- phyloseq::estimate_richness(ps, split = TRUE, measure = "Shannon")
  # Simpson index
  ps_alpha_div2 <- phyloseq::estimate_richness(ps, split = TRUE, measure = "Simpson")
  
  # Phylogenetic diversity
  # Sum all branch lengths in the phylogenetic tree
  PD <- sum(phyloseq::tree_layout(phyloseq::phy_tree(ps))$edgeDT[["edge.length"]])
  
  message("Storing diversity measures and network diagnostics in data frames.")
  
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
  
  message("Returning ps_network list of products.")
  
  # list of all objects to return
  ps_network <- list()
  ps_network$seqtabAll <- seqtabAll
  ps_network$ASVtab <- ASVtab
  ps_network$taxTab <- taxTab
  ps_network$phylo_tree <- phylo_tree
  ps_network$sampledata <- sampledata
  ps_network$ps <- ps
  ps_network$net <- net
  ps_network$results.samples <- results.samples
  ps_network$degree.samp <- degree.samp
  
  return(ps_network)
}








# function with arguments not in list form
ps.net2 <- function(fastq.Fs, fastq.Rs, # file paths for forward and reverse raw fastq files
                    filtFs, filtRs, # file paths for filtered and trimmed sequences
                    refFasta,
                    metadata = NULL, # samdf (sample metadata in list or data frame format: ex. meta <- list(sample.ID = 1, crop = "clover", site = 3))
                    truncLen = c(240,160),
                    maxN = 0, maxEE = c(2,2),
                    truncQ = 2, rm.phix = TRUE, # defining max 2 expected errors per read
                    compress = TRUE,
                    network.args = list(type = "taxa", # "samples"
                                        distance = "jaccard", # "jaccard", "unifrac", "wunifrac", DPCoA", "jsd"
                                        max.dist = 0.45,
                                        keep.isolates = TRUE),
                    seed = 100,
                    multithread = TRUE) {
  set.seed(100)
  
  # Extract sample names
  sampleNames <- sapply(strsplit(fastq.Fs, "_"), `[`, 2)
  sampleNames <- sapply(strsplit(sampleNames, "/"), `[`, 2)
  
  message("Filtering and trimming files.")
  
  # Filter forward and reverse reads, trim sequences to positions 240 and 160.
  out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                              truncLen = c(240,160),
                              maxN = 0, maxEE = c(2,2),
                              truncQ = 2, rm.phix = TRUE, # defining max 2 expected errors per read
                              compress = TRUE,
                              multithread = FALSE) # On Windows set multithread = FALSE) 
  
  message("Dereplicating files.")
  
  # Dereplicating identical sequences
  derepFs <- dada2::derepFastq(filtFs, verbose = TRUE)
  derepRs <- dada2::derepFastq(filtRs, verbose = TRUE)
  # Naming the derep-class objects by the sample names
  names(derepFs) <- sampleNames
  names(derepRs) <- sampleNames
  
  message("Learning sequencing errors pre-DADA2.")
  
  # Estimating sequencing error rates pre-DADA2
  errF <- dada2::learnErrors(filtFs, multithread=TRUE)
  errR <- dada2::learnErrors(filtRs, multithread=TRUE)
  
  message("DADA2 to infer sequences.")
  
  # DADA2 method to infer ASVs
  dadaFs <- dada2::dada(derepFs, err = errF, multithread = TRUE)
  dadaRs <- dada2::dada(derepRs, err = errR, multithread = TRUE)
  
  message("Merging forward and reverse sequence pairs. Constructing the sequence table.")
  
  # Merging sequence pairs into object "merged"
  merged <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs)
  # Constructing sequence table
  seqtabAll <- dada2::makeSequenceTable(merged)
  rownames(seqtabAll) <- sampleNames
  
  message("Removing chimeras to build the ASV table.")
  
  # Remove chimeras from the sequence table
  ASVtab <- dada2::removeBimeraDenovo(seqtabAll)
  
  message("Assigning taxonomy and building a taxonomy table.")
  
  refFASTA <- "./sh_general_release_dynamic_29.11.2022.fasta"
  # Build taxonomy table from the classified sequences
  taxTab <- dada2::assignTaxonomy(ASVtab,
                                  refFasta = refFASTA,
                                  multithread = TRUE)
  
  message("Labelling unclassified taxa.")
  
  # Remove the prefix from taxa names
  taxTab <- data.frame(row.names = row.names(taxTab),
                       Kingdom = stringr::str_replace(taxTab[,1], "k__",""),
                       Phylum = stringr::str_replace(taxTab[,2], "p__",""),
                       Class = stringr::str_replace(taxTab[,3], "c__",""),
                       Order = stringr::str_replace(taxTab[,4], "o__",""),
                       Family = stringr::str_replace(taxTab[,5], "f__",""),
                       Genus = stringr::str_replace(taxTab[,6], "g__",""),
                       Species = stringr::str_replace(taxTab[,7], "s__",""))
  
  # Rename "NA" elements to "Unclassified __"
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
  
  # Number replicates of unclassified taxa
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
  
  message("Multiple sequence alignment for the phylogeny.")
  
  # Multiple sequence alignment pre-tree
  seqs <- dada2::getSequences(ASVtab)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor=NA,verbose=FALSE)
  
  message("Fitting the phylogenetic tree.")
  
  # Phylogenetic tree construction
  phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- phangorn::dist.ml(phangAlign)
  treeNJ <- phangorn::NJ(dm) # Note: tip order != sequence order
  fit = phangorn::pml(treeNJ, data = phangAlign)
  fitGTR <- stats::update(fit, k = 4, inv = 0.2)
  fitGTR <- phangorn::optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                                rearrangement = "stochastic",
                                control = phangorn::pml.control(trace = 0))
  
  message("Storing phylogenetic outputs in a list.")
  
  phylo_tree <- list(phangAlign, dm, treeNJ, fit, fitGTR)
  
  message("Checking and reformatting the sample metadata.")
  
  metadata = as.list(read.csv(("C:/Users/elinrg26/OneDrive - UBC/Research/GitHub/MP microbiome/Network Analysis Workflow/Mapping file for ITS sequencing.csv")))
  
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
  
  message("Combining data into a phyloseq object.")
  
  # Now all the sequence and taxonomic data can be combined into a single phyloseq-class object
  tax = phyloseq::tax_table(as.matrix(taxTab))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(ASVtab,
                                               taxa_are_rows = FALSE),
                           phyloseq::sample_data(sample_info),
                           tax,
                           phyloseq::phy_tree(fitGTR$tree))
  
  message("Building the network from the phyloseq object.")
  
  # Construct the network using `ps`
  net <- phyloseq::make_network(ps, type = "taxa",
                                distance = "jaccard",
                                max.dist = 0.55,
                                keep.isolates = TRUE)
  sampledata <- as.data.frame(phyloseq::sample_data(ps)) # subset sample data as a data frame
  
  plot(net)
  
  message("Computing site-wide network diagnostics.")
  
  # Count vertices and edges in site network
  v <- igraph::vcount(net)
  e <- igraph::ecount(net)
  
  connectivity <- e/v   # avg number associations between ASVs
  connectance <- e/v^2  # fraction of edges out of all possible edges
  # Degree of each node
  degrees <- igraph::degree(net, v = igraph::V(net), mode = "all")
  
  message("Calculating diversity measures for each sample.")
  
  # Aggregate duplicate species
  species <- t(as.data.frame(ASVtab))  # new data frame for species counts
  species <- as.data.frame(species)
  species$Species <- taxTab$Species  # match species names to ASVs
  species <- stats::aggregate(x = species[, colnames(species) != "Species"],
                              by = list(species$Species), FUN = sum)
  rownames(species) <- species$Group.1
  keep.rows <- rownames(species)
  species <- as.data.frame(species[,-1], row.names = keep.rows)
  colnames(species) <- row.names(ASVtab)
  
  # Sum the presences in each sample for species richness
  rich <- nrow(species) # species count in sample
  # rich.rel <- rich/nrow(species) # species in sample/species in all samples   ### MAY NOT BE POSSIBLE (don't have all species in site so it would just be = 1)
  # total individuals from all species in each sample
  n <- apply(ASVtab, 1, function(x) sum(x))
  ASVs <- ncol(ASVtab) # total number of ASVs in site
  
  # Alpha species diversity measures
  # Shannon index
  ps_alpha_div <- phyloseq::estimate_richness(ps, split = TRUE, measure = "Shannon")
  # Simpson index
  ps_alpha_div2 <- phyloseq::estimate_richness(ps, split = TRUE, measure = "Simpson")
  
  # Phylogenetic diversity
  # Sum all branch lengths in the phylogenetic tree
  PD <- sum(phyloseq::tree_layout(phyloseq::phy_tree(ps))$edgeDT[["edge.length"]])
  
  message("Storing diversity measures and network diagnostics in data frames.")
  
  # Display sample measures in a data frame
  results.samples <- data.frame(ps_alpha_div, 
                                ps_alpha_div2, 
                                rich, #rich.rel, 
                                PD,
                                n.indivs = n, 
                                ASVs, 
                                vertices = v, 
                                edges = e, 
                                connectivity, 
                                connectance)
  if (!is.null(metadata)){
    results.samples <- merge(results.samples, sampledata)
    results.samples <- dplyr::relocate(Sample.ID, .before = Shannon) # move Sample.ID column to leftmost
  }
  
  degree.samp <- data.frame(ASVtab[1,], degrees, row.names = colnames(ASVtab))
  colnames(degree.samp)[1] <- sampledata$sample.ID
  
  message("Returning ps_network list of products.")
  
  # list of all objects to return
  ps_network <- list()
  ps_network$seqtabAll <- seqtabAll
  ps_network$ASVtab <- ASVtab
  ps_network$taxTab <- taxTab
  ps_network$phylo_tree <- phylo_tree
  ps_network$ps <- ps
  ps_network$net <- net
  ps_network$results.samples <- results.samples
  ps_network$degree.samp <- degree.samp
  ps_network$sampledata <- sampledata
  
  return(ps_network)
}




# FUNCTION: Dereplicates, infers ASVs (DADA2 method), merges sequence pairs, builds sequence table, removes chimeras, assigns taxonomy, generates a phylogenetic tree, combines all data into a phyloseq-class object, constructs a network plot, extracts measures of community structure and network topology.
# Returns: List with sequence table, table of ASVs, taxonomy table, phylogenetic tree, sample metadata, phyloseq object, network plot, data frame of richness/diversity metrics and network diagnostics for each sample.
# NOTE: "fastq.fwd" and "fastq.rev" are the specific file paths for the forward and reverse reads of each sample.

ps.net3 <- function(fastq.Fs, fastq.Rs, # file paths for forward and reverse raw fastq files
                   filtFs, filtRs, # file paths for filtered and trimmed sequences
                   refFasta, # file path of reference FASTA for taxonomic classification
                   metadata = NULL, # sample metadata in list or data frame format
                   # filter.args = list(truncLen = c(240,160),
                   #                    maxN = 0,
                   #                    maxEE = c(2,2),
                   #                    truncQ = 2,
                   #                    rm.phix = TRUE, # defining max 2 expected errors per read
                   #                    compress = TRUE,
                   #                    multithread = FALSE), # On Windows set multithread = FALSE
                   tree.args = list(k = 4, 
                                    inv = 0.2,
                                    model = "GTR", # "JC", "F81", "SYM", "GTR", etc. **SEE phangorn::optim.pml FOR ALL MODELS**
                                    rearrangement = "stochastic"), # "none", "NNI", "stochastic", "ratchet"
                   network.args = list(type = "taxa", # "samples"
                                       distance = "jaccard", # "jaccard", "unifrac", "wunifrac", DPCoA", "jsd"
                                       max.dist = 0.45, # default distance set by `make_network` function from `phyloseq` package is 0.4
                                       keep.isolates = TRUE), 
                   seed = 100,
                   multithread = TRUE) {
  
  
  
  # set.seed(seed)
  
  # Extract sample names, assuming file paths have format: FOL_DER/SAMPLENAME_XXX.fastq.gz OR FOL_DER/SAMPLENAME_XXX.fastq
  sampleNames <- sapply(strsplit(fastq.Fs, "_"), `[`, 2)
  sampleNames <- sapply(strsplit(sampleNames, "/"), `[`, 2)
  
  
  message("Dereplicating files.")
  
  # Dereplicating identical sequences
  derepFs <- dada2::derepFastq(filtFs, verbose = TRUE)
  derepRs <- dada2::derepFastq(filtRs, verbose = TRUE)
  # Naming the derep-class objects by the sample names
  names(derepFs) <- sampleNames
  names(derepRs) <- sampleNames
  
  message("Learning sequencing errors pre-DADA2.")
  
  # Estimating sequencing error rates pre-DADA2
  # The application of the DADA2 method follows the workflow by Callahan et al. (2017).
  # https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
  errF <- dada2::learnErrors(filtFs, multithread = multithread)
  errR <- dada2::learnErrors(filtRs, multithread = multithread)
  
  message("DADA2 to infer sequences.")
  
  # DADA2 method to infer ASVs
  dadaFs <- dada2::dada(derepFs, err = errF, multithread = multithread)
  dadaRs <- dada2::dada(derepRs, err = errR, multithread = multithread)
  
  message("Merging forward and reverse sequence pairs. Constructing the sequence table.")
  
  # Merging sequence pairs into object "merged"
  merged <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs)
  # Constructing sequence table
  seqtabAll <- dada2::makeSequenceTable(merged)
  rownames(seqtabAll) <- sampleNames
  
  message("Removing chimeras to build the ASV table.")
  
  # Remove chimeras from the sequence table
  ASVtab <- dada2::removeBimeraDenovo(seqtabAll)
  
  message("Taxonomic classification based on chosen reference FASTA.")
  
  # Assign taxonomy to the ASVs based on chosen reference FASTA and build a taxonomy table
  taxTab <- dada2::assignTaxonomy(ASVtab,
                                  refFasta = refFasta,
                                  multithread = multithread)
  
  message("Labelling unclassified taxa.")
  
  # Remove the prefix from taxa names
  # Removal of prefixes follow the workflow by Hui (2021).
  # https://www.yanh.org/2021/01/01/microbiome-r/#build-phyloseq-project
  taxTab <- data.frame(row.names = row.names(taxTab),
                       Kingdom = stringr::str_replace(taxTab[,1], "k__",""),
                       Phylum = stringr::str_replace(taxTab[,2], "p__",""),
                       Class = stringr::str_replace(taxTab[,3], "c__",""),
                       Order = stringr::str_replace(taxTab[,4], "o__",""),
                       Family = stringr::str_replace(taxTab[,5], "f__",""),
                       Genus = stringr::str_replace(taxTab[,6], "g__",""),
                       Species = stringr::str_replace(taxTab[,7], "s__",""))
  
  # Rename "NA" elements to "Unclassified __"
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
  
  # Number replicates of unclassified taxa
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
  
  message("Multiple sequence alignment for the phylogeny.")
  
  # Multiple sequence alignment pre-tree
  # Phylogenetic tree construction follows the workflow by Callahan et al. (2017).
  # https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
  seqs <- dada2::getSequences(ASVtab)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor=NA, verbose=FALSE)
  
  message("Fitting the phylogenetic tree.")
  
  # Phylogenetic tree construction
  phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- phangorn::dist.ml(phangAlign)
  treeNJ <- phangorn::NJ(dm) # Note: tip order != sequence order
  fit = phangorn::pml(treeNJ, data = phangAlign)
  fitGTR <- stats::update(fit, k = tree.args$k, inv = tree.args$inv) # consider sensitivity analysis on these values
  fitGTR <- phangorn::optim.pml(fitGTR, model = tree.args$model, optInv = TRUE, optGamma = TRUE,
                                rearrangement = tree.args$rearrangement,
                                control = phangorn::pml.control(trace = 0))
  
  message("Storing phylogenetic outputs in a list.")
  
  phylo_tree <- list()
  phylo_tree$phangAlign <- phangAlign
  phylo_tree$dm <- dm
  phylo_tree$treeNJ <- treeNJ
  phylo_tree$fit <- fit
  phylo_tree$fitGTR <- fitGTR
  
  message("Checking and reformatting the sample metadata.")
  
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
  
  message("Combining data into a phyloseq object.")
  
  # Now all the sequence and taxonomic data can be combined into a single phyloseq-class object
  tax = phyloseq::tax_table(as.matrix(taxTab))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(ASVtab,
                                               taxa_are_rows = FALSE),
                           phyloseq::sample_data(sample_info),
                           tax,
                           phyloseq::phy_tree(fitGTR$tree))
  
  message("Building the network from the phyloseq object.")
  
  # Construct the network using `ps`
  net <- phyloseq::make_network(ps,
                                type = network.args$type,
                                distance = network.args$distance,
                                max.dist = network.args$max.dist,
                                keep.isolates = network.args$keep.isolates)
  sampledata <- as.data.frame(phyloseq::sample_data(ps)) # subset sample_info as a data frame
  
  message("Computing site-wide network diagnostics.")
  
  # Count vertices and edges in site network
  v <- igraph::vcount(net)
  e <- igraph::ecount(net)
  
  connectivity <- e/v   # avg number associations between ASVs
  connectance <- e/v^2  # fraction of edges out of all possible edges
  # Degree of each node
  degrees <- igraph::degree(net, v = igraph::V(net), mode = "all")
  
  message("Calculating diversity measures for each sample.")
  
  # Aggregate duplicate species
  species <- phyloseq::tax_glom(ps, taxrank = 'Species', NArm = FALSE)
  
  # Sum the presences in each sample for species richness
  rich <- nrow(phyloseq::tax_table(species)) # species count in sample
  # total individuals from all species in each sample
  n <- apply(ASVtab, 1, function(x) sum(x))
  ASVs <- ncol(ASVtab) # total number of ASVs in site
  
  # Alpha species diversity measures
  # Shannon index
  ps_alpha_div <- phyloseq::estimate_richness(ps, split = TRUE, measure = "Shannon")
  # Simpson index
  ps_alpha_div2 <- phyloseq::estimate_richness(ps, split = TRUE, measure = "Simpson")
  
  # Phylogenetic diversity
  # Sum all branch lengths in the phylogenetic tree
  PD <- sum(phyloseq::tree_layout(phyloseq::phy_tree(ps))$edgeDT[["edge.length"]])
  
  message("Storing diversity measures and network diagnostics in data frames.")
  
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
  
  message("Returning ps_network list of products.")
  
  # list of all objects to return
  ps_network <- list()
  ps_network$seqtabAll <- seqtabAll
  ps_network$ASVtab <- ASVtab
  ps_network$taxTab <- taxTab
  ps_network$phylo_tree <- phylo_tree
  ps_network$sampledata <- sampledata
  ps_network$ps <- ps
  ps_network$net <- net
  ps_network$results.samples <- results.samples
  ps_network$degree.samp <- degree.samp
  
  return(ps_network)
}




ps.only <- function(fastq.Fs, fastq.Rs, # file paths for forward and reverse raw fastq files
                   filtFs, filtRs, # file paths for filtered and trimmed sequences
                   refFasta, # file path of reference FASTA for taxonomic classification
                   metadata = NULL, # sample metadata in list or data frame format
                   # filter.args = list(truncLen = c(240,160),
                   #                    maxN = 0,
                   #                    maxEE = c(2,2),
                   #                    truncQ = 2,
                   #                    rm.phix = TRUE, # defining max 2 expected errors per read
                   #                    compress = TRUE,
                   #                    multithread = FALSE), # On Windows set multithread = FALSE
                   tree.args = list(k = 4, 
                                    inv = 0.2,
                                    model = "GTR", # "JC", "F81", "SYM", "GTR", etc. **SEE phangorn::optim.pml FOR ALL MODELS**
                                    rearrangement = "stochastic"), # "none", "NNI", "stochastic", "ratchet"
                   seed = 100,
                   multithread = TRUE) {
  
  
  set.seed(seed)
  
  # Extract sample names, assuming file paths have format: FOL_DER/SAMPLENAME_XXX.fastq.gz OR FOL_DER/SAMPLENAME_XXX.fastq
  sampleNames <- sapply(strsplit(fastq.Fs, "_"), `[`, 2)
  sampleNames <- sapply(strsplit(sampleNames, "/"), `[`, 2)
  
  message("Dereplicating files.")
  
  # Dereplicating identical sequences
  derepFs <- dada2::derepFastq(filtFs, verbose = TRUE)
  derepRs <- dada2::derepFastq(filtRs, verbose = TRUE)
  # Naming the derep-class objects by the sample names
  names(derepFs) <- sampleNames
  names(derepRs) <- sampleNames
  
  message("Learning sequencing errors pre-DADA2.")
  
  # Estimating sequencing error rates pre-DADA2
  # The application of the DADA2 method follows the workflow by Callahan et al. (2017).
  # https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
  errF <- dada2::learnErrors(filtFs, multithread = multithread)
  errR <- dada2::learnErrors(filtRs, multithread = multithread)
  
  message("DADA2 to infer sequences.")
  
  # DADA2 method to infer ASVs
  dadaFs <- dada2::dada(derepFs, err = errF, multithread = multithread)
  dadaRs <- dada2::dada(derepRs, err = errR, multithread = multithread)
  
  message("Merging forward and reverse sequence pairs. Constructing the sequence table.")
  
  # Merging sequence pairs into object "merged"
  merged <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs)
  # Constructing sequence table
  seqtabAll <- dada2::makeSequenceTable(merged)
  rownames(seqtabAll) <- sampleNames
  
  message("Removing chimeras to build the ASV table.")
  
  # Remove chimeras from the sequence table
  ASVtab <- dada2::removeBimeraDenovo(seqtabAll)
  
  message("Taxonomic classification based on chosen reference FASTA.")
  
  # Assign taxonomy to the ASVs based on chosen reference FASTA and build a taxonomy table
  taxTab <- dada2::assignTaxonomy(ASVtab,
                                  refFasta = refFasta,
                                  multithread = multithread)
  
  message("Labelling unclassified taxa.")
  
  # Remove the prefix from taxa names
  # Removal of prefixes follow the workflow by Hui (2021).
  # https://www.yanh.org/2021/01/01/microbiome-r/#build-phyloseq-project
  taxTab <- data.frame(row.names = row.names(taxTab),
                       Kingdom = stringr::str_replace(taxTab[,1], "k__",""),
                       Phylum = stringr::str_replace(taxTab[,2], "p__",""),
                       Class = stringr::str_replace(taxTab[,3], "c__",""),
                       Order = stringr::str_replace(taxTab[,4], "o__",""),
                       Family = stringr::str_replace(taxTab[,5], "f__",""),
                       Genus = stringr::str_replace(taxTab[,6], "g__",""),
                       Species = stringr::str_replace(taxTab[,7], "s__",""))
  
  # Rename "NA" elements to "Unclassified __"
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
  
  # Number replicates of unclassified taxa
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
  
  message("Multiple sequence alignment for the phylogeny.")
  
  # Multiple sequence alignment pre-tree
  # Phylogenetic tree construction follows the workflow by Callahan et al. (2017).
  # https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
  seqs <- dada2::getSequences(ASVtab)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor=NA, verbose=FALSE)
  
  message("Fitting the phylogenetic tree.")
  
  # Phylogenetic tree construction
  phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- phangorn::dist.ml(phangAlign)
  treeNJ <- phangorn::NJ(dm) # Note: tip order != sequence order
  fit = phangorn::pml(treeNJ, data = phangAlign)
  fitGTR <- stats::update(fit, k = tree.args$k, inv = tree.args$inv) # consider sensitivity analysis on these values
  fitGTR <- phangorn::optim.pml(fitGTR, model = tree.args$model, optInv = TRUE, optGamma = TRUE,
                                rearrangement = tree.args$rearrangement,
                                control = phangorn::pml.control(trace = 0))
  
  message("Storing phylogenetic outputs in a list.")
  
  phylo_tree <- list()
  phylo_tree$phangAlign <- phangAlign
  phylo_tree$dm <- dm
  phylo_tree$treeNJ <- treeNJ
  phylo_tree$fit <- fit
  phylo_tree$fitGTR <- fitGTR
  
  message("Checking and reformatting the sample metadata.")
  
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
  
  message("Combining data into a phyloseq object.")
  
  # Now all the sequence and taxonomic data can be combined into a single phyloseq-class object
  tax = phyloseq::tax_table(as.matrix(taxTab))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(ASVtab,
                                               taxa_are_rows = FALSE),
                           phyloseq::sample_data(sample_info),
                           tax,
                           phyloseq::phy_tree(fitGTR$tree))
  
  message("Calculating diversity measures for each sample.")
  
  # Aggregate duplicate species
  species <- phyloseq::tax_glom(ps, taxrank = 'Species', NArm = FALSE)
  
  # Sum the presences in each sample for species richness
  rich <- nrow(phyloseq::tax_table(species)) # species count in sample
  # total individuals from all species in each sample
  n <- apply(ASVtab, 1, function(x) sum(x))
  ASVs <- ncol(ASVtab) # total number of ASVs in site
  
  # Alpha species diversity measures
  # Shannon index
  ps_alpha_div <- phyloseq::estimate_richness(ps, split = TRUE, measure = "Shannon")
  # Simpson index
  ps_alpha_div2 <- phyloseq::estimate_richness(ps, split = TRUE, measure = "Simpson")
  
  # Phylogenetic diversity
  # Sum all branch lengths in the phylogenetic tree
  PD <- sum(phyloseq::tree_layout(phyloseq::phy_tree(ps))$edgeDT[["edge.length"]])
  
  message("Storing diversity measures and network diagnostics in data frames.")
  
  # Display sample measures in a data frame
  results.samples <- data.frame(ps_alpha_div, 
                                ps_alpha_div2, 
                                rich, 
                                PD,
                                n.indivs = n, 
                                ASVs)
  
  results.samples <- cbind(results.samples, sample_info)
  
  message("Returning ps_network list of products.")
  
  # list of all objects to return
  ps_network <- list()
  ps_network$seqtabAll <- seqtabAll
  ps_network$ASVtab <- ASVtab
  ps_network$taxTab <- taxTab
  ps_network$phylo_tree <- phylo_tree
  ps_network$sampledata <- sampledata
  ps_network$ps <- ps
  ps_network$results.samples <- results.samples
  
  return(ps_network)
}
