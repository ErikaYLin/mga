\name{plot.mga}
\alias{plot.mga}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting methods for mga objects
}
\description{
Produces plots of \code{mga} objects (phylogenetic tree and co-occurrence network) using \code{\link[phyloseq]{phyloseq}} plotting methos.
}
\usage{
\S3method{plot}{mga}(mga, type = "tree", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mga}{A \code{mga} object built using \code{\link{mga}}.}
  \item{type}{The type of plot to draw. Possible types are:

  * \code{"tree"} for a phylogenetic tree pulled from \code{\link[phyloseq]{phy_tree}},
  * \code{"network"} for a co-occurrence network pulled from \code{\link[phyloseq]{make_network}} and plotted with \code{\link[phyloseq]{plot_network}}.}

  \item{...}{Additional plotting arguments and garphing parameters.}

}
}
\details{
\code{type = "network"} is only possible if \code{network = "TRUE"} when using \code{\link{mga}} to buld the plotting input.
}

%% \references{
%% ~put references to the literature/web site here ~
}
\author{Erika Y. Lin
%% }
%% \note{
%%  ~~further notes~~
%% }

\seealso{
\code{\link{mga}} for details on \code{mga-class} objects, \code{\link[phyloseq]{phy_tree}} for phylogenetic tree construction, and \code{\link[phyloseq]{make_network}} and \code{\link[phyloseq]{plot_network}} for details on making and plotting \code{igraph} objects through \code{phyloseq}.
}
%% \examples{
%% ##---- Should be DIRECTLY executable !! ----
%% ##-- ==>  Define data, use random,
%% ##--	or do  help(data=index)  for the standard data sets.

%% ## The function is currently defined as
%% function (x)
%% {
%%   }
%% }

