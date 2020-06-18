####Load posterior to test convergence####
tracePlots <- function(file, burnin=0, thinning=1, plot=TRUE, display=c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")){
  require(BTRTools)
  require(coda)
  
  rjout <- loadRJ(file, burnin = burnin, thinning = thinning)
  chain_out <- type.convert(rjout$rj_output)
  rownames(chain_out) = chain_out[,"It"]
  chain_out = chain_out[,-1]
  # Retrieve numerical
  index <- sapply(chain_out,function(x) is.numeric(x))
  chain <- mcmc(chain_out[,index])
  
  # plot the trace
  if(plot){
    plot(chain[,display])
  }
  
  # Just compute some statistics (autocorrelation...)
  cat("Effective sample size:","\n")
  print(effectiveSize(chain[,display]))
  
  # return results
  invisible(chain)
  
}


####function to obtain post probabilies of rate shifts on node####
return_pprob <- function(PP, threshold=0.5, cl = "nOrgnNRate"){
  nodes <- PP$data$descNode[which((PP$data[ , cl] / PP$niter) >= threshold)]
  pprobs <- PP$data[which((PP$data[ , cl] / PP$niter) >= threshold) , cl] / PP$niter
  
  result = list(nodes=nodes, pprobs=pprobs)
  return(result)
}

#uses treeio to combine tree topology and bayestraits posterior into single S4 object for plotting.
add_rjpp_to_tree <- function(rjpp_out){
  rjpp_data <- as_tibble(rjpp_out$data)
  timetree <- rjpp_out$meantree
  timetree$edge.length <- rjpp_data$orgBL[-1]
  timetree$root.time <- max(nodeHeights(timetree))
  rjpp_data_nodes <- rjpp_data %>% rename(., node=descNode) %>% mutate(., iters = rjpp_out$niter) %>% mutate(., ppRate = round(nOrgnNRate/iters,2))
  timetree <- treeio::as.treedata(timetree)
  treedata <- treeio::full_join(timetree, y = rjpp_data_nodes, by = "node")
  return(treedata)
}
