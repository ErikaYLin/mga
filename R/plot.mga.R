#' @export
# Function for plotting mga

plot.mga <- function(mga,
                     type = "tree", # "tree", "network", and eventually "errors"
                     group.species = TRUE, # uses mga$ps.species over mga$ps for tree
                     ...) {

  # Plot phylogenetic tree
  if (type == "tree") {
  plot(phyloseq::phy_tree(if (group.species){mga$ps.species} else{mga$ps}), ...)
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
