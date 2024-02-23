\name{drop_taxa}
\alias{drop_taxa}
\title{Prune mga-class objects by key defined taxa
}
\description{
Filters \code{mga-class} objects produced by \code{\link{mga}} for select taxa that are specified in a table with columns for the taxa names and taxonomic level at which to search.
}
\usage{
drop_taxa(mga, taxa = NULL, network = TRUE, network.args = list(type = "taxa", distance = "jaccard", max.dist = 0.45, keep.isolates = TRUE))
}
\arguments{
  \item{mga}{A \code{mga-class} object produced by \code{\link{mga}}.}
  \item{taxa}{A table of taxa to filter for. The table must have two columns: "taxon" containing the names of taxa to keep, and "group" with the taxonomic level of each taxon listed (i.e. "Class", "Family", "Species", etc.).}
  \item{network}{A logical argument that determines if the co-occurrence network analysis is included.}
  \item{network.args}{List of the parameters fed into \code{\link[phyloseq]{make_network}} from the \code{phyloseq} package used in constructing the co-occurrence network for each sample. See \code{\link[phyloseq]{make_network}} for additional arguments available.}

}
\value{
Returns the same \code{mga-class} object, filtered to keep only the taxa specified. See \code{\link{mga}} for the full list of objects returned in a \code{mga-class} object. Note that \code{seqtabAll}, \code{sampledata}, and \code{phylo_tree} remain unchanged from the original \code{mga} object fed into the function.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
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
