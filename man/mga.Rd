\name{mga}
\alias{mga}
\alias{mga.netfree}
\title{Process and analyze raw microbial genetic sequences}
\usage{
mga(fastq.Fs, fastq.Rs, filtFs, filtRs, refFasta, metadata = NULL, tree.args = list(k = 4, inv = 0.2, model = "GTR", rearrangement = "stochastic"), network.args = list(type = "taxa", distance = "jaccard", max.dist = 0.45, keep.isolates = TRUE), seed = 100, multithread = TRUE)
mga.netfree(fastq.Fs, fastq.Rs, filtFs, filtRs, refFasta, metadata = NULL, tree.args = list(k = 4, inv = 0.2, model = "GTR", rearrangement = "stochastic"), seed = 100, multithread = TRUE)
}
\description{
Simple microbiome analysis by processing raw genetic sequencing data to produce a \code{phyloseq-class} object with tabulated OTUs, taxonomic classification, phylogenetic trees, and sample metadata. \code{mga} also constructs co-occurrence networks, and computes metrics of community structure and network topology. \code{mga.netfree} excludes  co-occurrence network construction and metrics.
}
\arguments{
  \item{fastq.Fs, fastq.Rs}{File paths for forward (Fs) and reverse (Rs) raw fastq files.}
  \item{filtFs, filtRs}{File paths for filtered and trimmed forward (Fs) and reverse (Rs) sequences.}
  \item{refFasta}{A file path of reference FASTA to be used for taxonomic classification.}
  \item{metadata}{Sample metadata in list or data frame format.}
  \item{tree.args}{List of the parameters fed into functions from the \code{phangorn} package and \code{stats::update()} for constructing phylogenetic trees for each sample. See \code{phangorn} for a list of available arguments.}
  \item{network.args}{List of the parameters fed into \code{make_network()} from the \code{phyloseq} package used in constructing the co-occurrence network for each sample. See \code{make_network()} for additional arguments available.}
  \item{seed}{Any integer value used to set the seed.}
  \item{multithread}{A logical argument that determines whether (\code{TRUE}) or not (\code{FALSE}) multithreading and parallelization is used.}

}
\value{
Returns a list of list of the following objects:
    \item{seqtabAll}{A table containing the sequences inferred from each sample using the \code{dada2} package.}
    \item{ASVtab}{A table containing the OTU/ASV sequences in \code{seqtabAll} with chimeras removed.}
    \item{taxTab}{An object of class 'data.frame' with the full taxonomic classification assigned to each sequence.}
    \item{phylo_tree}{A list of the objects created and used in constructing the phylogeny (\code{phangAlign} of class 'phyDat', \code{dm} of class 'dist', \code{treeNJ} of class 'phylo', \code{fit} of class 'pml', \code{fitGTR} of class 'pml').}
    \item{sampledata}{A \code{phyloseq-class} list of the sample metadata provided.}
    \item{ps}{A \code{phyloseq-class} object holding the sequence table, taxonomy table, phylogenetic tree, and the sample metadata in the form of its component data classes: \code{otu_table-class}, \code{taxonomyTable-class}, \code{phylo-class}, \code{sample_data-class}. Storing these objects in a \code{phyloseq-class} object facilitates handling of the data in downstream analyses.}
    \item{net}{An \code{igraph} object for the co-occurrence network constructed.}
    \item{results.samples}{A 'data.frame' of sample metadata and the computed community diversity and network tolopology metrics for each sample: Shannon's diversity, Simpson's diversity, richness (rich), phylogenetic diversity (PD), total sequence count (n.indivs), OTU/ASV count (ASVs), number of vertices (vertices) and edges (edges) in the network, network connectivity, and network connectance.}
    \item{degree.samp}{A 'data.frame' of the degrees of (number of edges connected to) each vertex in the network for each sample.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{Erika Y. Lin
}
%%\note{
%%  ~~further notes~~
%%}

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
