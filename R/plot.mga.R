#' @export
# Function for plotting mga

plot.mga <- function(mga,
                     type = "tree", # "tree", "network", and eventually "errors"
                     ...) {
  
  # Plot phylogenetic tree
  if (type == "tree") {
  plot(phyloseq::phy_tree(mga$ps), ...)
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
