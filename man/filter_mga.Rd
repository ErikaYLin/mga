\name{filter_mga}
\alias{filter_mga}
\title{Prune mga-class objects by key defined taxa
}
\description{
Filters \code{mga-class} objects produced by \code{\link{mga}} for select taxa that are specified in a table with columns for the taxa names and taxonomic level at which to search. Selected taxa can keep kept, dropped, or both as separate \code{mga-class} objects.
}
\usage{
filter_mga(mga, keep = TRUE, drop = TRUE, taxa = NULL, group.taxa = "Species", network = TRUE, network.args = list(type = "taxa", distance = "jaccard", max.dist = 0.45, keep.isolates = TRUE))
}
\arguments{
  \item{mga}{A \code{mga-class} object produced by \code{\link{mga}}.}
  \item{keep}{Keep only selected taxa in output.}
  \item{drop}{Drop only selected taxa in output. If both \code{keep} and \code{drop} are TRUE, both objects are listed in the returned output.}
  \item{taxa}{A vector of taxa names to filter for, with species provided as underscored binomials (e.g. \code{c("Bacteroides, Staphylococcus, Escherichia_coli")}).

  Alternatively: A table of taxa to filter for. The table must have two columns: "taxon" containing the names of taxa to keep, and "group" with the taxonomic level of each taxon listed (i.e. "Class", "Family", "Species", etc.).}
  \item{group.taxa}{The \code{phyloseq} object and stored ASVs are agglomerated by the chosen taxonomic level ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", or "Species") for network construction and calculation of diversity indices. If \code{NULL}, the initial phyloseq object with separate ASVs is used. See \code{\link[phyloseq]{tax_glom}} for details. NOTE: Taxonomic level chosen for 'group.taxa' should be equal to or broader than the finest level of taxa selected for filtering.}
  \item{network}{A logical argument that determines if the co-occurrence network analysis is included.}
  \item{network.args}{List of the parameters fed into \code{\link[phyloseq]{make_network}} from the \code{phyloseq} package used in constructing the co-occurrence network for each sample. See \code{\link[phyloseq]{make_network}} for additional arguments available.}

}
\value{
Returns the same \code{mga-class} object, filtered to keep and/or drop only the taxa specified. See \code{\link{mga}} for the full list of objects returned in a \code{mga-class} object. If both \code{keep = TRUE} and \code{drop = TRUE}, a list of both \code{mga} objects is returned as 'mga_keep' and 'mga_drop'.

Note that \code{seqtabAll}, \code{sampledata}, and \code{phylo_tree} remain unchanged from the original \code{mga} object fed into the function. However, the updated phylogenetic tree can be extracted from the pruned \code{phyloseq-class} objects ('ps' and/or 'ps.taxa') for plotting using \code{link[mga]{plot.mga}}.

}
\author{Erika Y. Lin}

\note{
If none of the taxa selected for filtering are identified, the original \code{mga} object fed into the function will be returned. If there are very few occurrences of the taxa for filtering, the phylogenetic tree may be reduced to 1 or fewer branches, hence there will be no tree and no value for PD.
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
