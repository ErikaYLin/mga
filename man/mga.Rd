\name{mga}
\alias{mga}
\title{Process and analyze filtered microbial genetic sequences}
\usage{
mga(fastq.Fs, fastq.Rs, filtFs, filtRs, refFasta, metadata = NULL, make.unique = TRUE, group.taxa = "Species", tree.args = list(k = 4, inv = 0.2, model = "GTR", rearrangement = "stochastic"), network = TRUE, network.args = list(type = "taxa", distance = "jaccard", max.dist = 0.45, keep.isolates = TRUE), seed = 100, multithread = FALSE, verbose = TRUE)
}
\description{
Simple microbiome analysis by processing filtered genetic sequencing data to produce a \code{mga-class} object with tabulated OTUs/ASVs, taxonomic classification, phylogenetic trees, and sample metadata. \code{mga} also constructs co-occurrence networks, and computes metrics of community structure and network topology. \code{mga} makes use of the \code{\link[phyloseq]} and \code{\link[dada2]} workflows for sequence processing. Any primers should be trimmed from sequences prior to analysis.
}
\arguments{
  \item{fastq.Fs, fastq.Rs}{File paths for forward (Fs) and reverse (Rs) raw fastq files.}
  \item{filtFs, filtRs}{File paths for filtered and trimmed forward (Fs) and reverse (Rs) sequences.}
  \item{refFasta}{A file path of reference FASTA to be used for taxonomic classification. The reference library should be formatted for DADA2. A list of available DADA2-formatted reference FASTAs can be found at: \link{https://benjjneb.github.io/dada2/training.html}}
  \item{metadata}{Sample metadata in list or data frame format.}
  \item{make.unique}{A logical argument that determines if all unclassified taxa of the same higher-order taxon are to be considered as unique taxa. If \code{TRUE} unclassified taxa are made unique by numbering replicates.}
  \item{group.taxa}{The \code{phyloseq} object and stored ASVs are agglomerated by the chosen taxonomic level ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", or "Species") for network construction and calculation of diversity indices. If \code{NULL}, the initial phyloseq object with separate ASVs is used. See \code{\link[phyloseq]{tax_glom}} for details.}
  \item{tree.args}{List of the parameters fed into functions from the \code{phangorn} package and \code{stats::update()} for constructing phylogenetic trees for each sample. See \code{phangorn} for a list of available arguments.}
  \item{network}{A logical argument that determines if the co-occurrence network analysis is included.}
  \item{network.args}{List of the parameters fed into \code{make_network()} from the \code{phyloseq} package used in constructing the co-occurrence network for each sample. See \code{make_network()} for additional arguments available.}
  \item{seed}{Any integer value used to set the seed.}
  \item{multithread}{A logical argument that determines whether (\code{TRUE}) or not (\code{FALSE}) multithreading and parallelization is used.}
  \item{verbose}{\code{TRUE} enables \code{mga} to return messages for each processing step executed.}

}
\value{
Returns a \code{mga-class} object listing the following items:
    \item{seqtabAll}{A table containing the sequences inferred from each sample using the \code{dada2} package.}
    \item{ASVtab}{A table containing the OTU/ASV sequences in \code{seqtabAll} with chimeras removed.}
    \item{taxTab}{An object of class 'data.frame' with the full taxonomic classification assigned to each sequence.}
    \item{phylo_tree}{A list of the objects created and used in constructing the phylogeny (\code{phangAlign} of class 'phyDat', \code{dm} of class 'dist', \code{treeNJ} of class 'phylo', \code{fit} of class 'pml', \code{fitGTR} of class 'pml').}
    \item{sampledata}{A \code{phyloseq-class} list of the sample metadata provided.}
    \item{ps}{A \code{phyloseq-class} object holding the sequence table, taxonomy table, phylogenetic tree, and the sample metadata in the form of its component data classes: \code{otu_table-class}, \code{taxonomyTable-class}, \code{phylo-class}, \code{sample_data-class}. Storing these objects in a \code{phyloseq-class} object facilitates handling of the data in downstream analyses.}
    \item{ps.taxa}{The \code{phyloseq-class} object agglomerated by the selected taxonomic level. Only returned if a taxonomic level is specified for the `group.taxa` argument.}
    \item{net}{An \code{igraph} object for the co-occurrence network constructed. If `network = TRUE`. `ps.taxa` is used for network construction if `group.taxa` is specified.}
    \item{results.samples}{A 'data.frame' of sample metadata and the computed community diversity and network topology metrics for each sample: Shannon's diversity, Simpson's diversity, count/richness of taxa grouped at the level specified by `group.taxa` (rich),  phylogenetic diversity (PD), total sequence read count (reads), OTU/ASV count (ASVs), number of vertices (vertices) and edges (edges) in the network, network connectivity, and network connectance. If \code{group.taxa = NULL}, `rich` is equivalent to species count/richness. Phylogenetic diversity (PD) calculated by \code{mga} is equivalent to Faith's phylogenetic diversity, a frequently used term in similar applications, such as QIIME/QIIME2.}
    \item{degree.samp}{A 'data.frame' of the degrees of (number of edges connected to) each vertex in the network for each sample. If `network = TRUE`.}
}
\note{
To ensure that sample names are properly extracted, fastq file paths should be named following the format: FOL_DER/SAMPLENAME_XXX.fastq.gz OR FOL_DER/SAMPLENAME_XXX.fastq. Note that this may require manual curation.

Forward and reverse fastq files need to be manually filtered and trimmed before the analysis, as sequencing quality and error rates vary. See \code{\link[dada2]{filterAndTrim}} for details and examples. Additionally, sequencing companies may return sample files that still contain primers. Any primers should be removed prior to \code{mga} processing.
}

%%\references{
%% ~put references to the literature/web site here ~
%%}

\author{Erika Y. Lin
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.

%%## The function is currently defined as
%%function (x)
%%{
%%  }
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
