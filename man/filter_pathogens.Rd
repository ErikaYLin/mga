\name{filter_pathogens}
\alias{filter_pathogens}
\title{Filter phyloseq-class objects by key pathogenic taxa
}
\description{
Filters \code{phyloseq-class} objects produced by \code{mga()} and the \code{phyloseq} package for pathogenic taxa that are specified in a table with columns for the pathogenic taxa and taxonomic level at which to search.
}
\usage{
filter_pathogens(x, pathogens = NULL)
}
\arguments{
  \item{x}{A \code{phyloseq-class} object produced by \code{mga()} and the \code{phyloseq} package.}
  \item{pathogens}{A table of pathogenic taxa to filter for. The table must have two columns: "Pathogen" containing the names of pathogenic taxa, and "Taxon" with the taxonomic level of each taxon listed (i.e. "Class", "Family", "Species", etc.).}

}
\value{
Returns the same \code{phyloseq-class} object, filtered to keep only the taxa specified.
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
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
