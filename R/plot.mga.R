#' @export
# Function for plotting mga

plot.mga <- function(mga,
                     type = "tree", # "tree", "network", and eventually "errors"
                     group.taxa = TRUE, # uses mga$ps.taxa over mga$ps for tree
                     ...) {

  # Plot phylogenetic tree
  if (type == "tree") {
    # Do not run for group.taxa if ps.taxa does not exist
    if (group.taxa & is.null(mga$ps.taxa)) {
      stop("mga object not previously grouped by taxon. Check if 'ps.taxa' exists in mga object.")
    } else {
      plot(phyloseq::phy_tree(if (group.taxa){mga$ps.taxa} else{mga$ps}), ...)
    }
  } else

  # Plot co-occurrence network
  if (type == "network")  {
    if ("net" %in% names(mga)) {
      phyloseq::plot_network(mga$net, ...)
    } else {
      warning("No network in mga-class object to plot. Try rerunning mga() with 'network = TRUE'.")
    }
}
}
