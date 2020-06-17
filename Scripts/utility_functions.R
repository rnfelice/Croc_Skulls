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

#modified function for plotting branches colored by rate, derrived from 
#phytools function plotBranchByTrait
mytreebybranch<-function (tree, x, mode = c("edges", "tips", "nodes"), palette = "rainbow", nbin=15,
                          legend = TRUE, xlims = NULL, ...) 
{
  mode <- mode[1]
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (mode == "tips") {
    x <- c(x[tree$tip.label], fastAnc(tree, x))
    names(x)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)
    XX <- matrix(x[tree$edge], nrow(tree$edge), 2)
    x <- rowMeans(XX)
  }
  else if (mode == "nodes") {
    XX <- matrix(x[tree$edge], nrow(tree$edge), 2)
    x <- rowMeans(XX)
  }
  if (hasArg(tol)) 
    tol <- list(...)$tol
  else tol <- 1e-06
  if (hasArg(prompt)) 
    prompt <- list(...)$prompt
  else prompt <- FALSE
  if (hasArg(type)) 
    type <- list(...)$type
  else type <- "phylogram"
  if (hasArg(show.tip.label)) 
    show.tip.label <- list(...)$show.tip.label
  else show.tip.label <- TRUE
  if (hasArg(show.node.label)) 
    show.node.label <- list(...)$show.node.label
  else show.node.label <- FALSE
  if (hasArg(edge.width)) 
    edge.width <- list(...)$edge.width
  else edge.width <- 4
  if (hasArg(edge.lty)) 
    edge.lty <- list(...)$edge.lty
  else edge.lty <- 1
  if (hasArg(font)) 
    font <- list(...)$font
  else font <- 3
  if (hasArg(cex)) 
    cex <- list(...)$cex
  else cex <- par("cex")
  if (hasArg(adj)) 
    adj <- list(...)$adj
  else adj <- NULL
  if (hasArg(srt)) 
    srt <- list(...)$srt
  else srt <- 0
  if (hasArg(no.margin)) 
    no.margin <- list(...)$no.margin
  else no.margin <- TRUE
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (hasArg(label.offset)) 
    label.offset <- list(...)$label.offset
  else label.offset <- 0.01 * max(nodeHeights(tree))
  if (hasArg(underscore)) 
    underscore <- list(...)$underscore
  else underscore <- FALSE
  if (hasArg(x.lim)) 
    x.lim <- list(...)$x.lim
  else x.lim <- NULL
  if (hasArg(y.lim)) 
    y.lim <- list(...)$y.lim
  else y.lim <- if (legend && !prompt && type %in% c("phylogram", 
                                                     "cladogram")) 
    c(1 - 0.06 * length(tree$tip.label), length(tree$tip.label))
  else NULL
  if (hasArg(direction)) 
    direction <- list(...)$direction
  else direction <- "rightwards"
  if (hasArg(lab4ut)) 
    lab4ut <- list(...)$lab4ut
  else lab4ut <- NULL
  if (hasArg(tip.color)) 
    tip.color <- list(...)$tip.color
  else tip.color <- "black"
  if (hasArg(plot)) 
    plot <- list(...)$plot
  else plot <- TRUE
  if (hasArg(rotate.tree)) 
    rotate.tree <- list(...)$rotate.tree
  else rotate.tree <- 0
  if (hasArg(open.angle)) 
    open.angle <- list(...)$open.angle
  else open.angle <- 0
  if (is.function(palette)) 
    cols <- palette(n = nbin)
  else {
    if (palette == "heat.colors") 
      cols <- heat.colors(n = nbin)
    if (palette == "gray") 
      cols <- gray(nbin:1/nbin)
    if (palette == "rainbow") 
      cols <- rainbow(nbin, start = 0.7, end = 0)
  }
  if (is.null(xlims)) 
    xlims <- range(x) + c(-tol, tol)
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i + 
        1
    cols[i]
  }
  colors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  par(lend = 2)
  xx <- plot.phylo(tree, type = type, show.tip.label = show.tip.label, 
                   show.node.label = show.node.label, edge.color = colors, 
                   edge.width = edge.width, edge.lty = edge.lty, font = font, 
                   cex = cex, adj = adj, srt = srt, no.margin = no.margin, 
                   root.edge = root.edge, label.offset = label.offset, underscore = underscore, 
                   x.lim = x.lim, y.lim = y.lim, direction = direction, 
                   lab4ut = lab4ut, tip.color = tip.color, plot = plot, 
                   rotate.tree = rotate.tree, open.angle = open.angle, lend = 2, 
                   new = FALSE)
  if (legend == TRUE && is.logical(legend)) 
    legend <- round(0.3 * max(nodeHeights(tree)), 2)
  if (legend) {
    if (hasArg(title)) 
      title <- list(...)$title
    else title <- "trait value"
    if (hasArg(digits)) 
      digits <- list(...)$digits
    else digits <- 1
    if (prompt) 
      add.color.bar(legend, cols, title, xlims, digits, 
                    prompt = TRUE)
    else add.color.bar(legend, cols, title, xlims, digits, 
                       prompt = FALSE, x = par()$usr[1] + 0.05 * (par()$usr[2] - 
                                                                    par()$usr[1]), y = par()$usr[3] + 0.05 * (par()$usr[4] - 
                                                                                                                par()$usr[3]))
  }
  invisible(xx)
}