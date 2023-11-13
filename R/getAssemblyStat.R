## Helper function ##
getAssemblyStat <- function(contig.sizes, asm.index = 'assembly') {
  ## Calculate assembly statistics
  contig.sizes <- as.numeric(contig.sizes)
  contig.sizes <- sort(contig.sizes, decreasing = TRUE)
  total.contig.size <- sum(contig.sizes)
  cum.contig.size <- cumsum(contig.sizes)
  
  L50 <- sum((cum.contig.size) < (total.contig.size * 0.5)) + 1
  N50 <- contig.sizes[L50]
  
  L90 <- sum((cum.contig.size) < (total.contig.size * 0.9)) + 1
  N90 <- contig.sizes[L90]
  
  aun <- sum(contig.sizes * contig.sizes) / sum(contig.sizes)
  
  asm.metrics.df <- data.frame(ID = asm.index, N50 = N50, N90 = N90, AuN = aun, 
                               n.contigs = length(contig.sizes), 
                               min.ctg.size = min(contig.sizes),
                               max.ctg.size = max(contig.sizes))
  ## Return assembl statistics
  return( asm.metrics.df)
}