#Functions for finding network motifs
#source("./pattern_templates.R")

`%+%` <- function(x, y)  mapply(sum, x, y, MoreArgs = list(na.rm = TRUE)) #skip NAs when adding vectors

sbml2dnf <- function(SBML){
  "Converts the regulatory functions in an SBML objects attribute 'interactions' into their disjunctive normal forms. 
  Used to creative trinary adjacency matrices with trinaryadjmat(SBML). 
  For binary adjacency matrices, use only the SBML object in its unmodified form."
  saveNetwork(SBML, file = "tmpnet.txt", generateDNFs = T)
  SBML_DNF <- loadNetwork("tmpnet.txt", symbolic = T)
  file.remove("tmpnet.txt")
  return(SBML_DNF)
}

graph_to_amatrix <- function(g) {
  "Converts an igraph graph object to an adjacency matrix"
  a <- get.adjacency(g)
  mat <- matrix(as.numeric(a), ncol=ncol(a))
  for (d in 1:dim(mat)[2]){if (mat[d,d] >= 1){mat[d,d] <- 1}}
  return(mat)
}

conv2adjmat <- function(sbmlnet, inputcorrected = FALSE){
  "Converts an SBML object to adjacency matrix. 
  If inputcorrected = T, all input edges are set to zero. 
  A vertex v is said to be an input if it is only regulated by itself, 
  meaning the sum of column v in the adjacency matrix is one, with the only non-zero entry 
  being at position [v,v]."
  adjmat <- sapply(sbmlnet$interactions, function(gene) {v <- rep(0,length(sbmlnet$genes)); 
  v[gene$input] <- 1; return(v)})
  if (inputcorrected == TRUE){
    for (d in 1:dim(adjmat)[1]){
      if (adjmat[d,d] >= 1){adjmat[d,d] <- 1}
    }
    for (col in 1:dim(adjmat)[2]){
      if (sum(adjmat[,col]) == 1 & adjmat[col,col] == 1){
        adjmat[col,col] <-0
      }
    }
  }
  return(adjmat)
}

trinaryadjmatbygene <- function(SBML, genenr){
  #input from SBML
  #Max. 3 Layers in tree: Layer1 is always operator OR/AND, Layer2 can be atom or AND (mixed possible), 
  #                       Layer3 is always atom due to DNF structure
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  SBML_DNF <- sbml2dnf(SBML)
  ints <- SBML_DNF$interactions
  col <- rep(0, length(ints))
  tree <- ints[[genenr]]
  #Case I: Direct regulation
  if (tree$`type` == "atom"){
    col[tree$index] <- tree$negated + 1
  }
  #Case II: Direct AND of regulators
  else if (tree$`type` == "operator" && tree$operator == "&"){
    ops <- length(tree$operands)
    for (o in 1:ops){
      col[tree$operands[[o]]$index] <- tree$operands[[o]]$negated + 1
    }
  }
  #Case III (mixed, general case): OR of (either atoms or AND of regulators)
  else if (tree$`type` == "operator" && tree$operator == "|"){
    ops <- length(tree$operands)
    for (o in 1:ops){
      if (tree$operands[[o]]$`type` == "atom"){
        
        if ((col[tree$operands[[o]]$index] == 1 && sum(tree$operands[[o]]$negated + 1) == 2) 
            | (col[tree$operands[[o]]$index] == -1 && sum(tree$operands[[o]]$negated + 1) == 1)){
          print(paste("Warning: Position [", col[tree$operands[[o]]$index], ",", genenr, "] is regulated ambiguously.", sep=""))
        }
        
        col[tree$operands[[o]]$index] <- tree$operands[[o]]$negated + 1
      }
      if (tree$operands[[o]]$`type` == "operator"){
        nratoms <- length(tree$operands[[o]]$operands)
        for (A in 1:nratoms){
          
          if ((col[tree$operands[[o]]$operands[[A]]$index] == 1 && sum(tree$operands[[o]]$operands[[A]]$negated + 1) == 2) 
              | (col[tree$operands[[o]]$operands[[A]]$index] == -1 && tree$operands[[o]]$operands[[A]]$negated + 1)){
            print(paste("Warning: Position [", col[tree$operands[[o]]$operands[[A]]$index], ",", genenr, "] is regulated ambiguously.", sep=""))
          }
          
          col[tree$operands[[o]]$operands[[A]]$index] <- tree$operands[[o]]$operands[[A]]$negated + 1
        }
      }
    }
  }
  col[col == 2] <- -1
  #Check for discrepancies due to "redundant" regulations in SBML_DNF boolean formula
  diff <- adjmat[,genenr] - abs(col)
  if (sum(diff) == 0){return(col)}else{
    probrows <- which(diff != 0) #indices [probrows[i], genenr] must be corrected
    for (p in probrows){
      probname <- colnames(adjmat)[p]
      probnameneg <- paste("!", probname, sep="")
      if (grepl(probnameneg, SBML$interactions[[genenr]]$expression)){col[p] <- -1}
      else{col[p] <- 1}
    }
  }
  return(col)
}

trinaryadjmat <- function(SBML){
  "Takes an SBML object and returns a trinary adjacency matrix where a value of 1 at position (i,j)
  indicates that gene j is positively regulated by gene i whereas the value -1 indicates a negative regulation."
  dim <- length(SBML$genes)
  res <- matrix(NA, nrow = dim, ncol = dim)
  for (g in 1:dim){
    res[,g] <- trinaryadjmatbygene(SBML, g)
  }
  return(res)
}

generalfflfinder <- function(trinadjmat, ffltype="C1"){
  "Takes a trinary adjacency matrix and returns a vector whose length corresponds to the number of genes 
  in the network with each value indicating how many FFL motifs of the given type the gene participates in."
  possibleFFLtypes <- c("C1", "C2", "C3", "C4", "I1", "I2", "I3", "I4")
  if (ffltype %in% possibleFFLtypes == FALSE){stop("Not a valid FFL type, choose C1-C4 or I1-I4.")}
  #xyedege is positive in: 
  if (ffltype %in% c("C1", "C3", "I1", "I3")){xyedge <- +1} else {xyedge <- -1}
  #yzedege is positive in: 
  if (ffltype %in% c("C1", "C2", "I3", "I4")){yzedge <- +1} else {yzedge <- -1}
  #xzedege is positive in: 
  if (ffltype %in% c("C1", "C4", "I1", "I4")){xzedge <- +1} else {xzedge <- -1}
  
  dim <- dim(trinadjmat)[2]
  counter <- 0
  freqbygene <- rep(0, dim)
  for (g in 1:dim){
    #index g -> check outputs of g
    outputplus <- which(trinadjmat[g,] == xyedge) #FIRST EDGE X->Y
    outputplus <- outputplus[! outputplus %in% g]
    if (length(outputplus) > 0){
      for (oplus in outputplus){
        #index oplus -> check outputs of oplus
        outputminus <- which(trinadjmat[oplus,] == yzedge) #SECOND EDGE Y->Z
        outputminus <- outputminus[! outputminus %in% c(g, oplus)]
        if (length(outputminus) > 0){
          for (ominus in outputminus){
            #index ominus -> check outputs of ominus
            if (trinadjmat[g, ominus] == xzedge){ #THIRD EDGE X->Z
              counter <- counter +1
              if (length(unique(c(g, oplus, ominus))) == 3){ #3 distinct nodes in the motif
                freqbygene[g] <- freqbygene[g] + 1
                freqbygene[oplus] <- freqbygene[oplus] + 1
                freqbygene[ominus] <- freqbygene[ominus] + 1
              }
            }}}
      }
    }
  }
  return(freqbygene)
}

motif_gene_occurrence <- function(pattern_adjmat, graph_adjmat, gene){
  dim <- dim(graph_adjmat)[2]
  if (typeof(gene) %in% c("integer", "double")){
    if (gene <= dim) {gene <- colnames(graph_adjmat)[gene]} else {return("Gene number exceeds dimension of adjacency matrix.")}
  } else if (typeof(gene) == "character"){
    if (gene %in% colnames(graph_adjmat) == FALSE ){return("Gene does not exist.")}
  }
  genenr <- which(colnames(graph_adjmat) == gene)
  graph <- graph_from_adjacency_matrix(graph_adjmat)
  pattern_graph <- graph_from_adjacency_matrix(pattern_adjmat)
  counter <- 0
  patternlist <- subgraph_isomorphisms(pattern_graph, graph)
  l <- length(patternlist)
  pl <- dim(pattern_adjmat)[2]
  if (l == 0){return("No such pattern exists in the graph.")}
  else {
    motifnrvec <- c(NULL)
    for (m in 1:l){
      for (g in 1:pl){
        if (length(integer(patternlist[[m]][[g]])) == genenr){
          counter <- counter + 1
          motifnrvec <- append(motifnrvec, m)
        }
      }
    }}
  return(list(counter, motifnrvec))
  # Returns how often this motif appears overall, and the corresponding numbers as given by subraph_isomorphisms(pattern, graph)
  #Gives correct order of genes corresponding to X,Y,Z
}

motiffreqbygene <- function(pattern_adjmat, SBML){
  "Returns a vector showing how many network motifs of the given pattern any gene takes part in.
  Use generalfflfinder() function for searching feed-forward loops on a trinary matrix."
  graph_adjmat <- conv2adjmat(SBML)
  dim <- length(SBML$genes)
  freqs <- rep(NA, dim)
  if (length(subgraph_isomorphisms(graph_from_adjacency_matrix(pattern_adjmat), graph_from_adjacency_matrix(graph_adjmat))) == 0){
    {return(rep(0, dim))}}
  else {
    for (gene in 1:dim){
      freqs[gene] <- motif_gene_occurrence(pattern_adjmat, graph_adjmat, gene)[[1]]
    }}
  return(freqs)
}

totalgenes <- function(adjmat){
  #Returns total number of vertices (i.e. genes)
  return(dim(adjmat)[1])
}

generateuniquestates <- function(N=2^l, l, fixGenes, returnList = T){
  "Returns a matrix with N rows, each row being a random unique Boolean sequence of length l.
  fixGenes is an optional vector of states that should have a certain value. Needs to be of length l.
  If fixGenes[g] = 1 or 0 -> remove all states where grid[,g] is 0 or 1. If fixGenes[g] has any other value, ignore it."
  if (N > 2^l){return(paste("N needs to be smaller than 2^l=", 2^l, sep=""))}
  if (missing(fixGenes) ==F & missing(N) == T){N <- 2^(l-length(which(fixGenes %in% c(1,0) ==T)))}
  #return matrix with N random unique binary vectors of length l
  states <- rep(list(0:1), l)
  grid <- expand.grid(states)
  colnames(grid) <- seq(1,l)
  if (N==2 & l==1){
    if (returnList == F){
      return(grid)
    } else{
      return(lapply(seq_len(nrow(grid)), function(i) grid[i,]))
    }
  }
  
  if (missing(fixGenes)){
    randomstates <- sample(1:dim(grid)[1], N, replace=F)
    grid <- grid[randomstates,]
    rownames(grid) <- seq(1, dim(grid)[1])
    if (returnList == F){
      return(grid)
    } else {
      return(lapply(seq_len(nrow(grid)), function(i) grid[i,]))
    }
  }
  else {
    nrfixed <- length(which(fixGenes %in% c(1,0) ==T))
    maxpossible <- 2^(l-nrfixed)
    if (length(fixGenes) != l){return(paste("Vector fixGenes must have length l=", l, sep=""))}
    for (g in 1:l){
      if (fixGenes[g] == 1){grid <- grid[which(grid[,g]==1),]}
      else if (fixGenes[g] == 0){grid <- grid[which(grid[,g]==0),]}
    }
    if (dim(grid)[1] < N){return(paste("Error: Not enough states left. Must have N<=", maxpossible, ".", sep=""))}else{
      randomstates <- sample(1:dim(grid)[1], N, replace=F)
      grid <- grid[randomstates,]
      rownames(grid) <- seq(1, dim(grid)[1])
    }
    if (returnList == F){
      return(grid)
    } else {
      return(lapply(seq_len(nrow(grid)), function(i) grid[i,]))
    }
  }
}

returnTopRankings_readfiles <- function(measure="VBnDP", path, nets, n, P){
  #P = toppercent, adjust given network size => p
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML, inputcorrected = FALSE)
  graph <- graph_from_adjacency_matrix(adjmat)
  p <- ceiling((P/100)*dim(adjmat)[1])
  
  if (measure == "VB"){
    #VB <- VertexBetweenness(adjmat, normalized=TRUE)
    VB <- readRDS(paste0(path, "Network measures/net_",n,"_VB.RDS"))
    xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
    rankedtopnodes <- rank(-xselect_VB, ties.method="min")
  } else if (measure == "DP"){
    #DP <- DeterminativePower(SBML)
    DP <- readRDS(paste0(path, "Network measures/net_",n,"_DP.RDS"))
    xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
    rankedtopnodes <- rank(-xselect_DP, ties.method="min")
  } else if (measure == "Connectivity"){
    Connectivity <- singlemoduleZ(SBML)
    xselect_Connectivity <- sort(Connectivity, decreasing = TRUE)[1:p]
    rankedtopnodes <- rank(-xselect_Connectivity, ties.method="min")
  } else if (measure == "Hamming"){
    Hamming <- meanofminattrdistsmerged(SBML, method="trinary")
    yselect_Hamming <- sort(Hamming, decreasing = TRUE)[1:p]
    rankedtopnodes <- rank(-yselect_Hamming, ties.method="min")
  } else if (measure == "AttrLoss"){
    AttrLoss <- attrcoveragemerged(SBML, method="trinary")
    yselect_Loss <- sort(AttrLoss)[1:p]
    rankedtopnodes <- rank(yselect_Loss, ties.method="min")
  } else if (measure == "AttrGain"){
    AttrGain <- newattrgainbygenemerged(SBML, method="trinary")
    yselect_AttrGain <- sort(AttrGain, decreasing = TRUE)[1:p]
    rankedtopnodes <- rank(-yselect_AttrGain, ties.method="min")
  } else if (measure == "VBuDP"){
    #VB <- VertexBetweenness(adjmat, normalized=TRUE)
    VB <- readRDS(paste0(path, "Network measures/net_",n,"_VB.RDS"))
    xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
    rankedVB <- rank(-xselect_VB, ties.method="min")
    #DP <- DeterminativePower(SBML)
    DP <- readRDS(paste0(path, "Network measures/net_",n,"_DP.RDS"))
    xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
    rankedDP <- rank(-xselect_DP, ties.method="min")
    uniongenes <- union(names(xselect_VB), names(xselect_DP))
    rankedtopnodes <- rep(NA, length(uniongenes))
    names(rankedtopnodes) <- uniongenes
    for (g in 1:length(rankedtopnodes)){#For all nodes in selected union, average rankings if they are in both, otherwise copy single ranking
      VBindex <- which(names(rankedVB) == uniongenes[g])
      DPindex <- which(names(rankedDP) == uniongenes[g])
      print(VBindex)
      if (length(VBindex) == 0){rankedtopnodes[g] <- rankedDP[DPindex]}
      else if (length(DPindex) == 0){rankedtopnodes[g] <- rankedVB[VBindex]}
      else {rankedtopnodes[g] <- (rankedVB[VBindex] + rankedDP[DPindex])/2}
      rankedtopnodes <- sort(rankedtopnodes)
    }
  } else if (measure == "VBnDP"){
    #VB <- VertexBetweenness(adjmat, normalized=TRUE)
    VB <- readRDS(paste0(path, "Network measures/net_",n,"_VB.RDS"))
    xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
    rankedVB <- rank(-xselect_VB, ties.method="min")
    #DP <- DeterminativePower(SBML)
    DP <- readRDS(paste0(path, "Network measures/net_",n,"_DP.RDS"))
    xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
    rankedDP <- rank(-xselect_DP, ties.method="min")
    intersectgenes <- intersect(names(xselect_VB), names(xselect_DP))
    rankedtopnodes <- rep(NA, length(intersectgenes))
    names(rankedtopnodes) <- intersectgenes
    for (g in 1:length(rankedtopnodes)){#For all nodes in selected intersection, average their rankings in VB and DP
      VBindex <- which(names(rankedVB) == intersectgenes[g])
      DPindex <- which(names(rankedDP) == intersectgenes[g])
      rankedtopnodes[g] <- (rankedVB[VBindex] + rankedDP[DPindex])/2
    }
  } else if (measure == "DynAvg"){
    #Hamming <- meanofminattrdistsmerged(SBML, method="trinary")
    #AttrLoss <- attrcoveragemerged(SBML, method="trinary")
    #AttrGain <- newattrgainbygenemerged(SBML, method="trinary")
    #AvgDynRank <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
    AvgDynRank <- readRDS(paste0(path, "Network measures/net_",n,"_DynAvg.RDS"))
    rankedtopnodes <- sort(AvgDynRank)[1:p]
  }
  
  return(rankedtopnodes)
}

VertexBetweenness <- function(adjmat, normalized=T){
  " Returns the vertex betweenness for all n vertices in the network. 
  If normalized=T, each value is normalised by the factor 1/(n*n-3*n+2)."
  graph <- graph_from_adjacency_matrix(adjmat)
  #returns double with betweenness centrality for every vertex
  return(betweenness(graph, normalized = normalized))
}

flattenAttractor <- function(attractor, nr, equidistant=F){
  "Takes an attractor object from the BoolNet function getAttractors() and the index of a specific attractor 
  in this object. Returns a trinary vector whose length corresponds to the total number of genes in 
  the network where each value indicates a pattern of gene expression. 
  '1' and '0' indicate the gene being always turned on or always turned off in the attractor states respectively, 
  while '2' indicates any kind of oscillation. If \texttt{equidistant=T}, 
  the value of oscillating genes is set to '0.5' instead."
  totalattrs <- length(attractor[[2]])
  if (nr > totalattrs){return(paste("Error: There are only", totalattrs, "attractors in total."))}
  length <- length(attractor$stateInfo$genes)
  #Decimal zu Bin?r
  attr <- apply(attractor[[2]][[nr]]$involvedStates, MARGIN = 2, BoolNet:::dec2bin, len = length)
  #fixed oder oszillierend?
  fixedT <- apply(attr, MARGIN = 1, function(x) all(x == T))
  fixedF <- apply(attr, MARGIN = 1, function(x) all(x == F))
  
  attractorState <- rep(2, length)
  
  attractorState[which(fixedT)] <- 1
  attractorState[which(fixedF)] <- 0
  if (equidistant==F){return(attractorState)}
  else{
    return(replace(attractorState, attractorState==2, 0.5))
  }
  #modify to return correct result for attractor of length 1 (necessarily only 1 or 0)
}
probactivation <- function(attractor, nr, gene){
  totalattrs <- length(attractor[[2]])
  if (nr > totalattrs){return(paste("Error: There are only", totalattrs, "attractors in total."))}
  length <- length(attractor$stateInfo$genes)
  attr <- apply(attractor[[2]][[nr]]$involvedStates, MARGIN = 2, BoolNet:::dec2bin, len = length)
  genets <- attr[gene,] #time series of expression of GENE in given attractor NR
  probofexpr <- sum(genets)/length(genets)
  return(probofexpr)
}
probAttractor <- function(attractor, nr){
  "Takes an attractor object obtained using the BoolNet function getAttractors() and an attractor index. 
  Returns the probability for any given gene to be active during an attractor cycle."
  total <- length(attractor$attractors)
  dim <- length(attractor$stateInfo$genes)
  res <- rep(NA, dim)
  for (g in 1:dim){
    res[g] <- probactivation(attractor, nr, g)
  }
  return(res)
}

#Dynamic measure - Hamming distance
meanofminattrdists <- function(SBML, modification, modgene, method = "trinary"){
  if (method %in% c("trinary", "prob") == F){return("Not a valid method, choose trinary or prob.")}
  if (modification %in% c("overexpression", "knockout") == F){return("Not a valid modification, choose overexpression or knockout.")}
  attrs <- getAttractors(SBML, method = "sat.exhaustive")
  modSBML <- SBML
  if (modification == "overexpression"){modSBML$fixed[modgene] <- 1}
  else if (modification == "knockout"){modSBML$fixed[modgene] <- 0}
  modattrs <- getAttractors(modSBML, method = "sat.exhaustive")
  res <- rep(NA, length(modattrs$attractors))
  for (ma in 1:length(modattrs$attractors)){
    if (method == "trinary"){
      flattenedmod <- flattenAttractor(modattrs, ma, equidistant = T)
    } else {flattenedmod <- probAttractor(modattrs, ma)}
    flattenedmod <- flattenedmod[-modgene]
    dists <- rep(NA, length(attrs$attractors))
    for (a in 1:length(attrs$attractors)){
      if (method == "trinary"){
        flattenedorig <- flattenAttractor(attrs, a, equidistant = T)
      } else {flattenedorig <- probAttractor(attrs, a)}
      flattenedorig <- flattenedorig[-modgene]
      dist_ma_from_a <- sum(abs(flattenedorig - flattenedmod))
      dists[a] <- dist_ma_from_a
    }
    res[ma] <- min(dists)
  }
  return(mean(res))
}
meanofminattrdistsbygene <- function(SBML, modification, method = "trinary"){
  "Modifies every gene either by overexpression or knockout and returns the mean of the minimal distances 
  between a fixed modified and all unmodified attractors as measured by their overlaps, given either the 
  method 'trinary' or 'prob'."
  if (modification %in% c("overexpression", "knockout") == F){return("Not a valid modification, choose overexpression or knockout.")}
  res <- rep(NA, length(SBML$genes))
  for (g in 1:length(SBML$genes)){
    res[g] <- meanofminattrdists(SBML, modification, g, method)
  }
  return(res)
}
meanofminattrdistsmerged <- function(SBML, method = "trinary"){
  "Returns a vector choosing the worst-case output from either OE or KO for any gene, 
  i.e. the most impactful perturbation."
  res <- OE <- KO <- rep(NA, length(SBML$genes))
  for (g in 1:length(SBML$genes)){
    OE[g] <- meanofminattrdists(SBML, modification="overexpression", g, method)
    KO[g] <- meanofminattrdists(SBML, modification="knockout", g, method)
    res[g] <- max(OE[g], KO[g])
  }
  names(res) <- SBML$genes
  return(res)
}

#Dynamic measure - Attractor loss
attrcoverageaftermod <- function(SBML, modification, modgene, method = "trinary"){
  if (method %in% c("trinary", "prob") == F){return("Not a valid method, choose trinary or prob.")}
  if (modification %in% c("overexpression", "knockout") == F){return("Not a valid modification, choose overexpression or knockout.")}
  attrs <- getAttractors(SBML, method = "sat.exhaustive")
  modSBML <- SBML
  if (modification == "overexpression"){modSBML$fixed[modgene] <- 1}
  else if (modification == "knockout"){modSBML$fixed[modgene] <- 0}
  modattrs <- getAttractors(modSBML, method = "sat.exhaustive")
  #Check how many orig. attractors have a matching modattr. with distance = 0 (phenotype preserved, excluding modgene)
  losscounter <- 0 #count phenotypes that have no match in modSBML, finally divide by nr. of orig.attrs
  for (a in 1:length(attrs$attractors)){
    dists <- rep(NA, length(modattrs$attractors))
    if (method == "trinary"){
      flattenedorig <- flattenAttractor(attrs, a, equidistant = T)
    } else {flattenedorig <- probAttractor(attrs, a)}
    flattenedorig <- flattenedorig[-modgene]
    for(ma in 1:length(modattrs$attractors)){
      if (method == "trinary"){
        flattenedmod <- flattenAttractor(modattrs, ma, equidistant = T)
      } else {flattenedmod <- probAttractor(modattrs, ma)}
      flattenedmod <- flattenedmod[-modgene]
      dists[ma] <- sum(abs(flattenedorig-flattenedmod))
    }
    if (min(dists) > 0){losscounter <- losscounter + 1}
  }
  return(1 - (losscounter/length(attrs$attractors)))
}
attrcoverageaftermodbygene <- function(SBML, modification, method = "trinary"){
  "Modifies every gene either by overexpression or knockout and returns the percentage of original attractors 
  that are still existent in the modified network, that is those with a minimal distance of 0 to a modified 
  attractor, as measured by their overlaps, given either the \texttt{method} 'trinary' or 'prob'."
  if (modification %in% c("overexpression", "knockout") == F){return("Not a valid modification, choose overexpression or knockout.")}
  res <- rep(NA, length(SBML$genes))
  for (g in 1:length(SBML$genes)){
    res[g] <- attrcoverageaftermod(SBML, modification, g, method)
  }
  names(res) <- SBML$genes
  return(res)
}
attrcoveragemerged <- function(SBML, merging="min", method="trinary"){
  if (merging %in% c("min", "mean") == F){return("Error: merging must be 'min' or 'mean'.")}
  if (method %in% c("trinary", "prob") == F){return("Error: method must be 'trinary' or 'prob'.")}
  #merging=c("mean", min") - mean: average OE + KO coverage; min - worst case, pick lower coverage
  #method=c("trinary", "prob")
  res <- OE <- KO <- rep(NA, length(SBML$genes))
  for (g in 1:length(SBML$genes)){
    OE[g] <- attrcoverageaftermod(SBML, modification="overexpression", g, method=method)
    KO[g] <- attrcoverageaftermod(SBML, modification="knockout", g, method=method)
    if (merging=="mean"){res[g] <- (OE[g]+KO[g])/2}
    else if (merging=="min"){res[g] <- min(OE[g], KO[g])}else{return("Not a valid merging scheme.")}
  }
  names(res) <- SBML$genes
  return(res)
}

#Dynamic measure - Attractor gain
newattrgain <- function(SBML, modification, modgene, method="trinary"){
  if (method %in% c("trinary", "prob") == F){return("Not a valid method, choose trinary or prob.")}
  if (modification %in% c("overexpression", "knockout") == F){return("Not a valid modificiation!")}
  #mod every gene, count occurrences of NEW attractors, i.e. distance to ALL orig attrs is >0
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  size <- dim(adjmat)[1]
  attrs <- getAttractors(SBML, method = "sat.exhaustive")
  nrattrs <- length(attrs$attractors)
  modSBML <- SBML
  if (modification=="overexpression"){modSBML$fixed[modgene] <- 1} else if (modification=="knockout"){modSBML$fixed[modgene] <- 0}
  modattrs <- getAttractors(modSBML, method = "sat.exhaustive")
  nrmodattrs <- length(modattrs$attractors)
  res <- nrmodattrs #start by assuming all modattrs are new, -1 if a modattr has an orig. equivalent with distance =0
  #Check how many orig. attractors have a matching modattr. with distance = 0 (phenotype preserved, excluding modgene)
  for(ma in 1:nrmodattrs){
    dists <- rep(NA, length(attrs$attractors)) #distances of given ma to all available a
    if (method == "trinary"){
      flattenedmod <- flattenAttractor(modattrs, ma, equidistant = T)
    } else {flattenedmod <- probAttractor(modattrs, ma)}
    flattenedmod <- flattenedmod[-modgene]
    for (a in 1:nrattrs){
      if (method == "trinary"){
        flattenedorig <- flattenAttractor(attrs, a, equidistant = T)
      } else {flattenedorig <- probAttractor(attrs, a)}
      flattenedorig <- flattenedorig[-modgene]
      dists[a] <- sum(abs(flattenedorig-flattenedmod))
    }
    if (min(dists) == 0){res <- res - 1} #there exists an orig attr with dist=0 to ma -> ma is not a newly gained attr
  }
  return(res)
}
newattrgainbygene <- function(SBML, modification, method="trinary"){
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  size <- dim(adjmat)[1]
  res <- rep(NA, size)
  for (g in 1:size){
    res[g] <- newattrgain(SBML=SBML, modification=modification, modgene = g, method=method)
  }
  names(res) <- SBML$genes
  return(res)
}
newattrgainbygenemerged <- function(SBML, method="trinary"){
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  size <- dim(adjmat)[1]
  res <- rep(NA, size)
  OEimpact <- newattrgainbygene(SBML, modification="overexpression", method=method)
  KOimpact <- newattrgainbygene(SBML, modification="knockout", method=method)
  for (g in 1:size){
    res[g] <- max(OEimpact[g], KOimpact[g])
  }
  names(res) <- colnames(adjmat)
  return(res)
}

#Static measure - Determinative Power
binaryShannonEntropy <- function(p){
  "Get binary Shannon entropy h(p), probability p must be [0,1]."
  if (p < 0 | p > 1){return("Error: Probability p must be in range [0,1].")}
  if (p == 0 | p == 1){return(0)}
  return(-p*log2(p)-(1-p)*log2(1-p))
}
suppBoolFunc <- function(SBML, bfunc){
  size <- length(SBML$genes)
  if (bfunc > size | bfunc < 1){return(paste("Error: bfunc must be in range [1,", size, "].", sep=""))}
  #generate all relevant 2^input combinations 
  inputs <- SBML$interactions[[bfunc]]$input
  #print(paste("Function", bfunc, "has", length(inputs), "inputs."))
  allinputcombis <- generateuniquestates(N=2^length(inputs), l=length(inputs), returnList = F)
  colnames(allinputcombis) <- inputs
  #feed into bfunc, get mapping to 1 or 0
  #stringFunc <- SBML$interactions[[bfunc]]$expression
  outputs <- rep(NA, dim(allinputcombis)[1])
  #take allinputcombis[1,] -> map to vector of zeros except with indices "inputs" replaced with these values?
  for (i in 1:dim(allinputcombis)[1]){
    initialvec <- rep(0, size)
    initialvec[c(inputs)] <- allinputcombis[i,]
    initialvec <- unlist(initialvec)
    newstate <- unlist(stateTransition(SBML, state = initialvec, type="synchronous"))
    #Does bfunc(inputs==newstate) = 1?
    outputs[i] <- newstate[bfunc]
  }
  #delete all input combination from matrix which mapped to 0
  #remainder is supp f_bfunc (input1, input2,...)
  support <- allinputcombis[which(outputs == 1),]
  return(support)
}
sumsupp <- function(support, j, xj){
  if (xj %in% c(1,0) == F){return("Error: xj must be binary.")}
  #Get sum_x € supp f_i P(X=x | Xj=xj)
  #P(X=x | Xj=xj) = P(X=x, xj=xj)/P(xj=xj)
  Pxj <- 0.5 #by assumption, denominator of Bayesian expression always 1/2 (no bias)
  inputs <- dim(support)[2]
  #For how many B^n states in the support of bfunc is the value at position xj=xj?
  xjmatched <- length(which(support[, which(colnames(support) == j)] == xj))/(2^inputs)
  return(xjmatched/Pxj)
}
MutInfo <- function(SBML, fi, j){
  #MI(fi(X),Xj) = h(sumsuppfi p_x) - 
  # 0.5*h(sumsuppfi P(X=x|Xj=1)) - 
  # 0.5*h(sumsuppfi P(X=x|Xj=0))
  support <- suppBoolFunc(SBML, fi)
  if (typeof(support) == "integer"){ #what happens if BoolFunc fi has only one input, either pos. or neg.?
    term1 <- binaryShannonEntropy(1/2)
    #sum over suppfi px is 1/2 as there is one state in the support (either (1) or (0)) out of 2 possible states -> h(1/2) = 1
    term2 <- term3 <- binaryShannonEntropy(1)
    #either term2 (for neg. reg.) or term (for pos.reg.) don't exist anyway, only one element in support, only one input j
    #this input is either 1 or 0 once, the other case Xj=0 or Xj=1 never occurs
    #for the remaining term, e.g. term2: P(X=(1)|Xj=1) = P(1)/P(Xj=1) = 0.5/0.5 = 1 -> h(1)=0
    #--> single input BoolFuncs fi with only one regulator j always yield MutInfo=1-0.5*0-0.5*0=1
  } else {
    term1 <- binaryShannonEntropy(dim(support)[1]/2^(dim(support)[2]))
    #term1 = sum x € suppfi p_x = total nr of elements in suppfi / 2^inputs of fi
    term2 <- binaryShannonEntropy(sumsupp(support, j=j, xj=1))
    term3 <- binaryShannonEntropy(sumsupp(support, j=j, xj=0))
  }
  return(term1-0.5*term2-0.5*term3)
}
DP_j <- function(SBML, j){
  #get sum of MutInfo(fi(X), Xj) for all i that are regulated by j, i.e. all nodes having j as an input
  adjmat <- conv2adjmat(SBML)
  inputs <- which(adjmat[j,] == 1) #which nodes have j as a common input? need fi from these, given Xj=1 or Xj=0...
  #print(inputs)
  #print(paste("Function j=", j, " has inputs: ", toString(inputs), ".", sep=""))
  res <- 0
  for (i in inputs){
    MI_fi_Xj <- MutInfo(SBML, fi = i, j = j)
    res <- res +  MI_fi_Xj
  }
  return(res)
}
DeterminativePower <- function(SBML){
  size <- length(SBML$genes)
  DPvec <- rep(NA, size)
  for (j in 1:size){
    DPvec[j] <- DP_j(SBML, j)
  }
  names(DPvec) <- SBML$genes
  return(DPvec)
}

#Node-wise static measures
ShimbelIndex <- function(SBML){
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML)
  graph <- graph_from_adjacency_matrix(adjmat)
  size <- dim(adjmat)[1]
  si <- rep(NA, size)
  for (i in 1:size){
    Shimbel <- igraph::shortest_paths(graph, from=i, to=V(graph))
    si[i] <- sum(unlist(lapply(Shimbel$vpath, FUN = length))) #sum of lengths of shortest paths from node i to all nodes
  }
  names(si) <- colnames(adjmat)
  return(si)
}
kCorebyvertex <- function(adjmat, mode="all"){
  "Returns a named vector, showing the highest value of k for which each gene belongs to a k-core. 
  Depending on whether mode is set to 'in', 'out', or 'all', only incoming, outgoing or all edges 
  are considered for the vertices' ranks. Additionally, if an integer value is given for plotk the corresponding 
  k-core will be highlighted in a plot."
  graph <- graph_from_adjacency_matrix(adjmat)
  if (mode %in%  c("all", "out", "in")){
    return(coreness(graph, mode))}
  else {return("Not a valid mode. Please choose among 'all', 'out' or 'in'.")}
}
vertexdegrees <- function(adjmat, mode = "all"){
  "Returns a vector containing the ranks of all vertices. 
  Depending on whether 'mode' is chosen as 'in' or 'out', only incoming or outcoming edges are counted. 
  If 'mode' is chosen as 'all', both will be counted."
  dim <- dim(adjmat)[2]
  if (is.null(dim)){return(0)}
  degrees <- c(NA, dim)
  if (mode == "in"){
    for (j in 1:dim){degrees[j] <- sum(adjmat[,j])}
    return(degrees)
  } else if (mode == "out"){
    for (j in 1:dim){degrees[j] <- sum(adjmat[j,])}
    return(degrees)
  } else if (mode == "all"){
    for (j in 1:dim){degrees[j] <- sum(adjmat[,j])}
    for (j in 1:dim){degrees[j] <- degrees[j] + sum(adjmat[j,])}
    return(degrees)
  } else {return("Not a valid mode. Choose in, out or all")}
}
laplacian <- function(adjmat, rank = "all"){
  "Returns the Laplacian matrix of the graph, using the total vertex ranks. 
  If rank is chosen as 'in' or 'out', only the in- or out-degrees are chosen as the ranks."
  vertdegs <- vertexdegrees(adjmat, mode = rank)
  #uses in + out degree as total vertex degree for directed graph
  #
  degreematrix <- matrix(0, dim(adjmat)[1], dim(adjmat)[2])
  for (i in 1:dim(adjmat)[1]){
    degreematrix[i,i] <- vertdegs[i]
  }
  L <- degreematrix - adjmat
  return(L)
}
resistancedistance <- function(adjmat, v1, v2){
  "Returns the resistance between vertices v1 and v2 if the graph is treated as
  an electrical circuit with each edge being a 1 Ohm resistor."
  L <- laplacian(adjmat)
  i <- v1
  j <- v2
  Linv <- solve(L)
  r_ij <- Linv[i,i] + Linv[j,j] - Linv[i,j] - Linv[j,i]
  names(r_ij) <- NULL
  return(r_ij)
}
resistancedistancematrix <- function(adjmat){
  "Returns the resistance distance matrix for all genes."
  L <- laplacian(adjmat)
  d <- dim(adjmat)[2]
  Linv <- solve(L)
  R <- matrix(0, d, d)
  for (i in 1:d){
    for (j in 1:d){
      R[i,j] <- Linv[i,i] + Linv[j,j] - Linv[i,j] - Linv[j,i]
    }
  }
  return(R)
}
meanresistancesbygene <- function(adjmat){
  "Returns a vector containing the mean resistances of every gene to all other genes."
  rdmat <- resistancedistancematrix(adjmat)
  means <- apply(rdmat, 1, mean)
  names(means) <- colnames(adjmat)
  return(means)
}

nroutlinks <- function(SBML){
  size <- length(SBML$genes)
  outlinks <- rep(NA, size)
  adjmat <- conv2adjmat(SBML)
  for (g in 1:size){
    outlinks[g] <- sum(adjmat[g,])
  }
  names(outlinks) <- SBML$genes
  return(outlinks)
}

#Generate z-scores (connectivity) normed by the network itself (not over all networks)
totaldegree <- function(SBML){
  adjmat <- conv2adjmat(SBML)
  size <- length(SBML$genes)
  k <- rep(NA, size)
  z <- rep(NA, size)
  for (i in 1:size){
    k[i] <- sum(adjmat[i,]) + sum(adjmat[,i])
    if (adjmat[i,i] == 1){k[i] <- k[i] - 1} 
  }
  names(k) <- colnames(adjmat)
  return(k)
}
singlemoduleZ <- function(SBML){
  adjmat <- conv2adjmat(SBML)
  size <- length(SBML$genes)
  k <- rep(NA, size)
  z <- rep(NA, size)
  for (i in 1:size){
    k[i] <- sum(adjmat[i,]) + sum(adjmat[,i])
    if (adjmat[i,i] == 1){k[i] <- k[i] - 1} 
    #Use -1 to only count self-loop once, not as both an outgoing and incoming edge
  }
  Kavg <- mean(k)
  Ksigma <- sd(k)
  
  for (i in 1:size){z[i] <- (k[i]-Kavg)/Ksigma}
  names(z) <- SBML$genes
  return(z)
}

#Generate z-scores (connectivity) normed over all networks
globalNormZ <- function(pathtonets, n){
  nets <- mixedsort(dir(path, pattern = ".txt"))
  alltotaldegrees <- c()
  for (N in 1:length(nets)){
    netfile <- paste(pathtonets, nets[N], sep = "")
    SBML <- loadNetwork(netfile)
    adjmat <- conv2adjmat(SBML)
    size <- length(SBML$genes)
    k <- rep(NA, size)
    z <- rep(NA, size)
    for (i in 1:size){
      k[i] <- sum(adjmat[i,]) + sum(adjmat[,i])
      if (adjmat[i,i] == 1){k[i] <- k[i] - 1} 
    }
    alltotaldegrees <- append(alltotaldegrees, k)
    #collect node degrees of all nodes across networks, get mean and sd of this for normalisation
  }
  
  Kavg <- mean(alltotaldegrees)
  Ksigma <- sd(alltotaldegrees)
  
  print(c(Kavg, Ksigma))
  
  netfile <- paste(pathtonets, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML)
  size <- length(SBML$genes)
  k <- rep(NA, size)
  z <- rep(NA, size)
  for (i in 1:size){
    k[i] <- sum(adjmat[i,]) + sum(adjmat[,i])
  }
  
  for (i in 1:size){z[i] <- (k[i]-Kavg)/Ksigma}
  names(z) <- SBML$genes
  return(z)
}

# Comparison of scores in VBnDP and connectivity
percentile <- function(x, perc, fullrange=TRUE){
  #fullrange=F => works as before but with sign changed! --> ecdf(-x)(-perc)
  #fullrange=T => highest value in x or above always yields 1, lowest value in x or below always yields 0 ==> NEW SCALING!
  if (fullrange==F){
    return(ecdf(-x)(-perc))
  }
  else {
    #get number of different values in vector x, find position # pos where perc falls in here -> 1 - (#pos - 1)/(nrunique-1)
    nrunique <- length(unique(x))
    if (nrunique == 1){
      print("Warning: No variation is given vector, only 1 unique score.")
      return(1)
    }
    sortedx <- sort(x)
    if (perc %in% sortedx){
      nrpos <- which(unique(sortedx) == perc)
    }
    else {
      perc <- max(x[x <= perc]) #whatever actual entry in the vector is closest to the given value of perc (round down)
      nrpos <- which(unique(sortedx) == perc)
    }
    percscore <- 1 - (nrpos - 1)/(nrunique - 1)
    #print(percscore)
    return(percscore)
  }
}

#Suggest targets, return PM sorted by mismatch
suggestTargets <- function(n, threshold, loaddata=TRUE){
  #get percentiles, rankings
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML, inputcorrected = FALSE)
  p <- ceiling((threshold/100)*dim(adjmat)[1])
  
  if (loaddata==TRUE){
    VBresults <- DPresults <- Zresults <- rep(list(NA), length(nets))
    VB <- readRDS(paste0(path, "/Network measures/", "net_", n, "_VB.RDS"))
    DP <- readRDS(paste0(path, "/Network measures/", "net_", n, "_DP.RDS"))
    Z <- readRDS(paste0(path, "/Network measures/", "net_", n, "_Z.RDS"))
    # Hamming <- readRDS(paste0(path, "/Network measures/", "net_", n, "_Hamming.RDS")) #Hamming distance
    # AttrLoss <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrLoss.RDS")) #Attractor loss
    # AttrGain <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrGain.RDS")) #Attractor gain
    # Hamming <- Hammingresults[[n]]
    # AttrLoss <- AttrLossresults[[n]]
    # AttrGain <- AttrGainresults[[n]]
    # AvgDynRank <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
  } else {
    VB <- VBresults[[n]]
    DP <- DPresults[[n]]
    Z <- Zresults[[n]]
  }
  xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
  xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
  VBRank <- rank(-VB, ties.method="min")
  DPRank <- rank(-DP, ties.method="min")
  name <- intersect(names(xselect_VB), names(xselect_DP))
  VBnDP <- rep(NA, length(name))
  names(VBnDP) <- name
  Z_inVBnDP <- Z[which(names(Z) %in% names(VBnDP))]
  for (g in 1:length(VBnDP)){
    RankA <- VBRank[which(names(VBRank) == names(VBnDP)[g])]
    RankB <- DPRank[which(names(DPRank) == names(VBnDP)[g])]
    VBnDP[g] <- mean(c(RankA, RankB))
  }
  ZmatchVBnDP <- Z_inVBnDP[match(names(VBnDP), names(Z_inVBnDP))]
  StatZQuantRatio <- rep(NA, length(VBnDP))
  names(StatZQuantRatio) <- names(VBnDP)
  for (g in 1:length(VBnDP)){
    StatZQuantRatio[g] <- percentile(VBnDP, VBnDP[g]) - percentile(-ZmatchVBnDP, -ZmatchVBnDP[g]) #low number=high rank
  }
  gatekeepers <- which(StatZQuantRatio > 0)
  
  #define sets
  NSindices <- which(colnames(adjmat) %in% setdiff(names(VB), names(VBnDP)))
  hubindices <- setdiff(which(Z >= 2.5), NSindices) #Guimera et al.
  nonhubindices <- setdiff(seq(1:length(SBML$genes)), hubindices)
  PMindices <- which(colnames(adjmat) %in% setdiff(names(gatekeepers), colnames(adjmat)[c(NSindices)]))
  NMindices <- setdiff(seq(1:length(SBML$genes)), c(NSindices, PMindices))
  
  targets <- colnames(adjmat)
  names(targets)[hubindices] <- "Hub"
  names(targets)[nonhubindices] <- "Non-Hub"
  targets <- targets[PMindices]
  decreasingMismatchOrder <- names(sort(StatZQuantRatio[gatekeepers], decreasing = TRUE))
  
  
  #Result: Recommended intervention targets, highest mismatch between VBnDP and connectivity rankings first
  targets <- targets[match(decreasingMismatchOrder, targets)]
  if (length(targets)>0){
    print(paste0("The recommended intervention targets (ordered by highest mismatch first) are:"))
    return(targets)
  } else {
    return("There are no gatekeeper nodes in this network using the provided selection threshold.")
  }
}

#Calculate all static and dynamic measures for all networks in the "Networks" folder, save the results in "Networks/Network measures/"
calculateAllMeasures <- function(nets, path, measure2calc="VB", saveResults=TRUE, savepath=paste0(path, "/Network measures/SensitivitySpecificity/")){
  results <- rep(list(NA),length(nets))
  #list with entry for every net, save results of measures here, return a list of lists of these values
  for (n in 1:length(nets)){
    #Check if subdirectory Network measures already exists, if not create it
    ifelse(!dir.exists(file.path(savepath)), dir.create(file.path(savepath)), FALSE)
    print(paste0("Net nr. ", n))
    netfile <- paste(path, nets[n], sep = "")
    SBML <- loadNetwork(netfile)
    adjmat <- conv2adjmat(SBML, inputcorrected = F)
    graph <- graph_from_adjacency_matrix(adjmat)
    
    if (measure2calc=="VB"){
      VB <- VertexBetweenness(adjmat, normalized=T)
      results[[n]] <- VB
      if (saveResults==T){saveRDS(VB, file=paste0(savepath, "net_", n, "_VB.RDS"))}
    } else if (measure2calc=="DP"){
      DP <- DeterminativePower(SBML)
      results[[n]] <- DP
      if (saveResults==T){saveRDS(DP, file=paste0(savepath, "net_", n, "_DP.RDS"))}
    } else if (measure2calc=="Z"){
      Z <- singlemoduleZ(SBML)
      results[[n]] <- Z
      if (saveResults==T){saveRDS(Z, file=paste0(savepath, "net_", n, "_Z.RDS"))}
    } else if (measure2calc=="Hamming"){
      Hamming <- meanofminattrdistsmerged(SBML, method="trinary")
      results[[n]] <- Hamming
      if (saveResults==T){saveRDS(Hamming, file=paste0(savepath, "net_", n, "_Hamming.RDS"))}
    } else if (measure2calc=="AttrLoss"){
      AttrLoss <- attrcoveragemerged(SBML, merging="min", method="trinary")
      results[[n]] <- AttrLoss
      if (saveResults==T){saveRDS(AttrLoss, file=paste0(savepath, "net_", n, "_AttrLoss.RDS"))}
    } else if (measure2calc=="AttrGain"){
      AttrGain <- newattrgainbygenemerged(SBML, method="trinary")
      results[[n]] <- AttrGain
      if (saveResults==T){saveRDS(AttrGain, file=paste0(savepath, "net_", n, "_AttrGain.RDS"))}
    } else if (measure2calc=="DynAvg"){
      Hamming <- meanofminattrdistsmerged(SBML, method="trinary")
      AttrLoss <- attrcoveragemerged(SBML, merging="min", method="trinary")
      AttrGain <- newattrgainbygenemerged(SBML, method="trinary")
      DynAvg <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), 
                               rank(-AttrGain, ties.method="min")), na.rm=TRUE)
      results[[n]] <- DynAvg
      if (saveResults==T){saveRDS(DynAvg, file=paste0(savepath, "net_", n, "_DynAvg.RDS"))}
    }
    
  }#end loop over nets
  return(results)
}

#Return nodes scoring in the top P percent of rankings in a given network according to some measure
returnTopRankings <- function(measure="VBnDP", path, nets, n, P, loaddata=TRUE){
  #P = toppercent, adjust given network size => p
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML, inputcorrected = FALSE)
  graph <- graph_from_adjacency_matrix(adjmat)
  p <- ceiling((P/100)*dim(adjmat)[1])
  
  if (measure == "VB"){
    if (loaddata==TRUE){VB <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_VB.RDS"))}
    else {VB <- VertexBetweenness(adjmat, normalized=TRUE)}
    xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
    rankedtopnodes <- rank(-xselect_VB, ties.method="min")
  } else if (measure == "DP"){
    if (loaddata==TRUE){DP <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_DP.RDS"))}
    else {DP <- DeterminativePower(SBML)}
    xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
    rankedtopnodes <- rank(-xselect_DP, ties.method="min")
  } else if (measure == "Connectivity"){
    if (loaddata==TRUE){Connectivity <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_Z.RDS"))}
    else {Connectivity <- singlemoduleZ(SBML)}
    xselect_Connectivity <- sort(Connectivity, decreasing = TRUE)[1:p]
    rankedtopnodes <- rank(-xselect_Connectivity, ties.method="min")
  } else if (measure == "Hamming"){
    if (loaddata==TRUE){Hamming <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_Hamming.RDS"))}
    else{Hamming <- meanofminattrdistsmerged(SBML, method="trinary")}
    yselect_Hamming <- sort(Hamming, decreasing = TRUE)[1:p]
    rankedtopnodes <- rank(-yselect_Hamming, ties.method="min")
  } else if (measure == "AttrLoss"){
    if (loaddata==TRUE){AttrLoss <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_AttrLoss.RDS"))}
    else {AttrLoss <- attrcoveragemerged(SBML, method="trinary")}
    yselect_Loss <- sort(AttrLoss)[1:p]
    rankedtopnodes <- rank(yselect_Loss, ties.method="min")
  } else if (measure == "AttrGain"){
    if (loaddata==TRUE){AttrGain <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_AttrGain.RDS"))}
    else {AttrGain <- newattrgainbygenemerged(SBML, method="trinary")}
    yselect_AttrGain <- sort(AttrGain, decreasing = TRUE)[1:p]
    rankedtopnodes <- rank(-yselect_AttrGain, ties.method="min")
  } else if (measure == "DynAvg"){
    if (loaddata==TRUE){
      AttrGain <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_AttrGain.RDS"))
      AttrLoss <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_AttrLoss.RDS"))
      Hamming <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_Hamming.RDS"))
    } else {
      AttrGain <- newattrgainbygenemerged(SBML, method="trinary")
      AttrLoss <- attrcoveragemerged(SBML, method="trinary")
      Hamming <- meanofminattrdistsmerged(SBML, method="trinary")
    }
    RankedGain <- rank(-AttrGain, ties.method="min")
    RankedLoss <- rank(AttrLoss, ties.method="min")
    RankedHamming <- rank(-Hamming, ties.method="min")
    DynAvg <- colMeans(rbind(RankedGain, RankedLoss, RankedHamming))
    yselect_DynAvg <- sort(DynAvg)[1:p]
    rankedtopnodes <- yselect_DynAvg
  }else if (measure == "VBuDP"){
    if (loaddata==TRUE){
      VB <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_VB.RDS"))
      DP <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_DP.RDS"))
    } else {
      VB <- VertexBetweenness(adjmat, normalized=TRUE)
      DP <- DeterminativePower(SBML)
    }
    xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
    rankedVB <- rank(-xselect_VB, ties.method="min")
    xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
    rankedDP <- rank(-xselect_DP, ties.method="min")
    uniongenes <- union(names(xselect_VB), names(xselect_DP))
    rankedtopnodes <- rep(NA, length(uniongenes))
    names(rankedtopnodes) <- uniongenes
    for (g in 1:length(rankedtopnodes)){#For all nodes in selected union, average rankings if they are in both, otherwise copy single ranking
      VBindex <- which(names(rankedVB) == uniongenes[g])
      DPindex <- which(names(rankedDP) == uniongenes[g])
      print(VBindex)
      if (length(VBindex) == 0){rankedtopnodes[g] <- rankedDP[DPindex]}
      else if (length(DPindex) == 0){rankedtopnodes[g] <- rankedVB[VBindex]}
      else {rankedtopnodes[g] <- (rankedVB[VBindex] + rankedDP[DPindex])/2}
      rankedtopnodes <- sort(rankedtopnodes)
    }
  } else if (measure == "VBnDP"){
    if (loaddata==TRUE){
      VB <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_VB.RDS"))
      DP <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_DP.RDS"))
    } else {
      VB <- VertexBetweenness(adjmat, normalized=TRUE)
      DP <- DeterminativePower(SBML)
    }
    xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
    rankedVB <- rank(-xselect_VB, ties.method="min")
    xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
    rankedDP <- rank(-xselect_DP, ties.method="min")
    intersectgenes <- intersect(names(xselect_VB), names(xselect_DP))
    rankedtopnodes <- rep(NA, length(intersectgenes))
    names(rankedtopnodes) <- intersectgenes
    for (g in 1:length(rankedtopnodes)){#For all nodes in selected intersection, average their rankings in VB and DP
      VBindex <- which(names(rankedVB) == intersectgenes[g])
      DPindex <- which(names(rankedDP) == intersectgenes[g])
      rankedtopnodes[g] <- (rankedVB[VBindex] + rankedDP[DPindex])/2
    }
  } else if (measure == "DynAvg"){
    if (loaddata==TRUE){AvgDynRank <- readRDS(paste0("/home/weidner/Desktop/BN Statics/plotscripts_data/Networks/Network measures/net_",n,"_DynAvg.RDS"))}
    else {
      Hamming <- meanofminattrdistsmerged(SBML, method="trinary")
      AttrLoss <- attrcoveragemerged(SBML, method="trinary")
      AttrGain <- newattrgainbygenemerged(SBML, method="trinary")
      AvgDynRank <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
    }
    rankedtopnodes <- sort(AvgDynRank)[1:p]
  }
  
  return(rankedtopnodes)
}

#Calculate sensitivity and specificity for a given static measure xs, saves the results in "Networks/Network measures/SensitivitySpecificity/"
#Returns list of vectors of length 2*N for N networks, first N entries  contains sensitivity vectors, second N contain specificity vectors
calculateSensSpec <- function(statmeasure="VBnDP", dynfunc="Hamming", nets, path, loaddata=FALSE, saveResults=TRUE, savepath=paste0(path, "/Network measures/")){
  Sensresults <- Specresults <- rep(list(NA), length(nets))
  if (statmeasure == "VB"){
    xs <- 1
  } else if (statmeasure == "DP"){
    xs <- 2
  } else if (statmeasure == "VBuDP"){
    xs <- 3
  } else if (statmeasure == "VBnDP"){
    xs <- 4
  } else {print("Not a valid measure.")}
  
  if (dynfunc == "Hamming"){
    ys <- 1
  } else if (dynfunc == "AttrLoss"){
    ys <- 2
  } else if (dynfunc == "AttrGain"){
    ys <- 3
  } else if (dynfunc == "DynAvg"){
    ys <- 4
  } else {print("Not a valid measure.")}
  
  print("Check availability of directories")
  ifelse(!dir.exists(file.path(savepath)), dir.create(file.path(savepath)), FALSE)
  VBresults <- DPresults <- Zresults <- Hammingresults <- AttrLossresults <- AttrGainresults <- rep(list(NA), 35)
  
  if (loaddata==FALSE){
    print("Calculating VB")
    VBresults <- calculateAllMeasures(nets, path, measure2calc="VB", saveResults=FALSE)
    print("Calculating DP")
    DPresults <- calculateAllMeasures(nets, path, measure2calc="DP", saveResults=FALSE)
    print("Calculating Connectivity")
    Zresults <- calculateAllMeasures(nets, path, measure2calc="Z", saveResults=FALSE)
    
    if (ys==1){
      print("Calculating Hamming distance")
      Hammingresults <- calculateAllMeasures(nets, path, measure2calc="Hamming", saveResults=FALSE) 
    } else if (ys==2){
      print("Calculating Attractor loss")
      AttrLossresults <- calculateAllMeasures(nets, path, measure2calc="AttrLoss", saveResults=FALSE) 
    } else if (ys==3){
      print("Calculating Attractor gain")
      AttrGainresults <- calculateAllMeasures(nets, path, measure2calc="AttrGain", saveResults=FALSE)
    } else if (ys==4){
      #DynAvg - need all individual dynamic measures
      print("Calculating Hamming distance")
      Hammingresults <- calculateAllMeasures(nets, path, measure2calc="Hamming", saveResults=FALSE)  
      print("Calculating Attractor loss")
      AttrLossresults <- calculateAllMeasures(nets, path, measure2calc="AttrLoss", saveResults=FALSE) 
      print("Calculating Attractor gain")
      AttrGainresults <- calculateAllMeasures(nets, path, measure2calc="AttrGain", saveResults=FALSE)
    }
  } else {
    for (n in 1:length(nets)){
      VBresults[[n]] <- readRDS(paste0(path, "Network measures/net_", n, "_VB.RDS"))
      DPresults[[n]] <- readRDS(paste0(path, "Network measures/net_", n, "_DP.RDS"))
      Zresults[[n]] <- readRDS(paste0(path, "Network measures/net_", n, "_Z.RDS"))
      Hammingresults[[n]] <- readRDS(paste0(path, "Network measures/net_", n, "_Hamming.RDS")) #Hamming distance
      AttrLossresults[[n]] <- readRDS(paste0(path, "Network measures/net_", n, "_AttrLoss.RDS")) #Attractor loss
      AttrGainresults[[n]] <- readRDS(paste0(path, "Network measures/net_", n, "_AttrGain.RDS")) #Attractor gain
    }
  }#end loaddata
  
  for (n in 1:length(nets)){
    print(paste0("Net ", n))
    print(paste0("ys=",ys))
    sensit <- rep(NA, 100)
    specif <- rep(NA, 100)
    
    VB <- VBresults[[n]]
    DP <- DPresults[[n]]
    Z <- Zresults[[n]]
    if (ys==1){
      Hamming <- Hammingresults[[n]]
    } else if (ys==2){
      AttrLoss <- AttrLossresults[[n]]
    } else if (ys==3){
      AttrGain <- AttrGainresults[[n]] 
    } else if (ys==4){
      Hamming <- Hammingresults[[n]]
      AttrLoss <- AttrLossresults[[n]]
      AttrGain <- AttrGainresults[[n]] 
      DynAvg <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
    }
    
    for (toppercent in seq(1:100)){
      #print(paste("p=",toppercent,sep=""))
      netfile <- paste(path, nets[n], sep = "")
      SBML <- loadNetwork(netfile)
      adjmat <- conv2adjmat(SBML, inputcorrected = FALSE)
      genenames <- colnames(adjmat)
      size <- dim(adjmat)[1]
      colnames(adjmat)
      
      #if picking given percentage of high scoring genes instead of mean or mean+sd as threshold
      p <- ceiling((toppercent/100)*dim(adjmat)[1])
      
      #VB
      xselect_VB <- names(sort(VB, decreasing = TRUE)[1:p])
      #DP
      xselect_DP <- names(sort(DP, decreasing = TRUE)[1:p])
      #Z
      xselect_Z <- names(sort(Z, decreasing = TRUE)[1:p])
      if (ys==1){
        #Hamming
        yselect_Hamming <- names(sort(Hamming, decreasing = TRUE)[1:p])
      }
      if (ys==2){
        #AttrLoss
        yselect_AttrLoss <- names(sort(AttrLoss)[1:p])
      }
      if (ys==3){
        #AttrGain
        yselect_AttrGain <- names(sort(AttrGain, decreasing = TRUE)[1:p])
      }
      if (ys==4){
        #DynAvg
        yselect_DynAvg <- names(sort(DynAvg)[1:p])
      }
      
      ##### CHOOSE RELEVANT DYNAMIC MEASURE #####
      if (ys == 1){yselect <- yselect_Hamming}
      if (ys == 2){yselect <- yselect_AttrLoss}
      if (ys == 3){yselect <- yselect_AttrGain}
      if (ys == 4){yselect <- yselect_DynAvg}
      ##### CHOOSE RELEVANT STATIC MESAURE #####
      if (xs == 1){xselect <- xselect_VB}  #VB
      if (xs == 2){xselect <- xselect_DP}  #DP
      if (xs == 3){xselect <- union(xselect_VB, xselect_DP)}  #VB v DP
      if (xs == 4){xselect <- intersect(xselect_VB, xselect_DP)}  #VB n DP
      ##### Calculate sensitivity and specificity #####
      s_high <- xselect
      s_low <- setdiff(genenames, xselect)
      d_high <- yselect
      d_low <- setdiff(genenames, yselect)
      
      sensit[toppercent] <- length(intersect(s_high, d_high))/length(d_high)
      specif[toppercent] <- length(intersect(s_low, d_low))/length(d_low)
      sensit[which(is.nan(sensit))] <- NA
      specif[which(is.nan(specif))] <- NA
      
    }#loop over threshold T
    #Option: Save generated vector sensit, specif
    if (saveResults==TRUE){
      saveRDS(sensit, file=paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", ys ,".RDS"))
      saveRDS(specif, file=paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", ys ,".RDS"))
    }
    #Save 100er-vecs for every network and dynamic function ys
    #list of length n for following plot functions to use
    Sensresults[[n]] <- sensit
    Specresults[[n]] <- specif
    
  }#loop over nets
  return(append(Sensresults, Specresults))
}

#Calculate maximal Mutual Information along all simple paths to hubs, saves the results in "Networks/Network measures/MaxMI/
#saves after individual network if saveaftereverynet=T, otherwise only saves once all networks have been calculated
calcMaxMItoHubs <- function(nets, path, threshold=73, loaddata=TRUE, saveaftereverynet=FALSE, saveResults=TRUE){
  result <- rep(list(NA), 3)
  HubDefinition <- "Guimera"
  PMvalues <- NMvalues <- NSvalues <- c()
  maxMIbynet_PM <- maxMIbynet_NM <- maxMIbynet_NS <- rep(NA, length(nets))
  ifelse(!dir.exists(file.path(path, "/Network measures/MaxMI/")), 
         dir.create(file.path(path, "/Network measures/MaxMI/")), FALSE)
  
  if (loaddata==FALSE){
    VBresults <- calculateAllMeasures(nets, path, measure2calc="VB", saveResults=FALSE)
    DPresults <- calculateAllMeasures(nets, path, measure2calc="DP", saveResults=FALSE)
    Zresults <- calculateAllMeasures(nets, path, measure2calc="Z", saveResults=FALSE)
    Hammingresults <- calculateAllMeasures(nets, path, measure2calc="Hamming", saveResults=FALSE)
    AttrLossresults <- calculateAllMeasures(nets, path, measure2calc="AttrLoss", saveResults=FALSE)
    AttrGainresults <- calculateAllMeasures(nets, path, measure2calc="AttrGain", saveResults=FALSE)
  }
  
  for (n in 1:length(nets)){
    print(n)
    ##### get percentiles, VBnDP rankings etc #####
    netfile <- paste(path, nets[n], sep = "")
    SBML <- loadNetwork(netfile)
    adjmat <- conv2adjmat(SBML, inputcorrected = FALSE)
    graph <- graph_from_adjacency_matrix(adjmat)
    p <- ceiling((threshold/100)*dim(adjmat)[1])
    if (loaddata==TRUE){
      VB <- readRDS(paste0(path, "/Network measures/", "net_", n, "_VB.RDS"))
      DP <- readRDS(paste0(path, "/Network measures/", "net_", n, "_DP.RDS"))
      Z <- readRDS(paste0(path, "/Network measures/", "net_", n, "_Z.RDS"))
      Hamming <- readRDS(paste0(path, "/Network measures/", "net_", n, "_Hamming.RDS")) #Hamming distance
      AttrLoss <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrLoss.RDS")) #Attractor loss
      AttrGain <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrGain.RDS")) #Attractor gain
    } else {
      VB <- VBresults[[n]]
      DP <- DPresults[[n]]
      Z <- Zresults[[n]]
      Hamming <- Hammingresults[[n]]
      AttrLoss <- AttrLoss[[n]]
      AttrGain <- AttrGain[[n]]
    }
    
    AvgDynRank <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
    xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
    xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
    VBRank <- rank(-VB, ties.method="min")
    DPRank <- rank(-DP, ties.method="min")
    name <- intersect(names(xselect_VB), names(xselect_DP))
    VBnDP <- rep(NA, length(name))
    names(VBnDP) <- name
    Z_inVBnDP <- Z[which(names(Z) %in% names(VBnDP))]
    for (g in 1:length(VBnDP)){
      RankA <- VBRank[which(names(VBRank) == names(VBnDP)[g])]
      RankB <- DPRank[which(names(DPRank) == names(VBnDP)[g])]
      VBnDP[g] <- mean(c(RankA, RankB))
    }
    ZmatchVBnDP <- Z_inVBnDP[match(names(VBnDP), names(Z_inVBnDP))]
    StatZQuantRatio <- rep(NA, length(VBnDP))
    names(StatZQuantRatio) <- names(VBnDP)
    for (g in 1:length(VBnDP)){
      StatZQuantRatio[g] <- percentile(VBnDP, VBnDP[g]) - percentile(-ZmatchVBnDP, -ZmatchVBnDP[g]) #low number=high rank
    }
    gatekeepers <- which(StatZQuantRatio > 0)
    
    ##### define sets to compare for given property #####
    NSindices <- which(colnames(adjmat) %in% setdiff(names(AvgDynRank), names(VBnDP)))
    hubindices <- which(Z >= 2.5) #Z>=2.5 hub definition by Guimera et al.
    nonhubindices <- setdiff(seq(1:length(SBML$genes)), hubindices)
    PMindices <- which(colnames(adjmat) %in% setdiff(names(gatekeepers), colnames(adjmat)[c(NSindices)]))
    NMindices <- setdiff(seq(1:length(SBML$genes)), c(NSindices, PMindices))
    
    ##### Max MI along simple paths normed by path length from Node->Hub #####
    print(paste(length(hubindices), "hubs in NET", n))
    if (length(hubindices)>0){
      ##### PM #####
      if (length(PMindices)>0){
        #sumMI <- 0
        for (g in PMindices){
          #print(paste0("Gene index ", g))
          paths <- all_simple_paths(graph, from=g, to=hubindices, mode="out")
          #get last element for every [[p]] paths entry, only retain those where this is index=h
          
          for (h in hubindices){
            if (g == h){next}
            lastelems <- which(sapply(paths, tail, 1) == h)
            pathstoh <- rep(list(NA), length(lastelems))
            #print(paste(length(lastelems), "paths to hub", h, "out of", length(paths), "paths to all hubs"))
            print("PM")
            if (length(pathstoh)> 0){
              
              for (el in 1:length(lastelems)){
                pathstoh[[el]] <- paths[[lastelems[el]]]
              }
              
              pathMIs <- rep(NA, length(pathstoh))
              print(paste(length(pathstoh), "paths from index", g, "to", h))
              #print(pathMIs)
              for (pa in 1:length(pathstoh)){ #p=list of gene names in path
                p <- pathstoh[[pa]]
                MIpath <- 0
                nodenames <- names(p)
                nrnodes <- length(nodenames)
                #which(colnames(adjmat) %in% nodenames)
                nodenames
                nodeindices <- rep(NA, length(nodenames))
                for (v in 1:length(nodeindices)){
                  nodeindices[v] <- which(colnames(adjmat) == nodenames[v])
                }# edge indices going from g->h
                for (e in 1:(length(nodeindices)-1)){
                  MI <- MutInfo(SBML, fi=nodeindices[e+1], j=nodeindices[e])
                  #sum up along path, avg by pathlength, then add to sumMI
                  MIpath <- MIpath + MI
                }
                MIpath <- MIpath/(length(nodeindices)-1)
                pathMIs[pa] <- MIpath
                #print(MIpath)
              }#paths loop
              
              #print(pathMIs)
              
              #sumMI <- sumMI + max(pathMIs)
              PMvalues <- append(PMvalues, max(pathMIs))
            }#if length paths>0
          }#end hub loop
        }#end PM loop
        #sumMI <- sumMI/(length(PMindices)*length(hubindices))
        #maxMIbynet_PM[n] <- sumMI
      }#if PMs exist
      
      ##### NM #####
      if (length(NMindices)>0){
        #sumMI <- 0
        for (g in NMindices){
          print(paste0("Gene index ", g))
          paths <- all_simple_paths(graph, from=g, to=hubindices, mode="out")
          #get last element for every [[p]] paths entry, only retain those where this is index=h
          
          for (h in hubindices){
            if (g == h){next}
            lastelems <- which(sapply(paths, tail, 1) == h)
            pathstoh <- rep(list(NA), length(lastelems))
            #print(paste(length(lastelems), "paths to hub", h, "out of", length(paths), "paths to all hubs"))
            print("NM")
            if (length(pathstoh)> 0){
              
              for (el in 1:length(lastelems)){
                pathstoh[[el]] <- paths[[lastelems[el]]]
              }
              
              pathMIs <- rep(NA, length(pathstoh))
              print(paste(length(pathstoh), "paths from index", g, "to", h))
              #print(pathMIs)
              for (pa in 1:length(pathstoh)){ #p=list of gene names in path
                p <- pathstoh[[pa]]
                MIpath <- 0
                nodenames <- names(p)
                nrnodes <- length(nodenames)
                #which(colnames(adjmat) %in% nodenames)
                nodenames
                nodeindices <- rep(NA, length(nodenames))
                for (v in 1:length(nodeindices)){
                  nodeindices[v] <- which(colnames(adjmat) == nodenames[v])
                }# edge indices flowing from g->h
                for (e in 1:(length(nodeindices)-1)){
                  MI <- MutInfo(SBML, fi=nodeindices[e+1], j=nodeindices[e])
                  #sum these up along the path, avg by pathlength, then add to sumMI
                  MIpath <- MIpath + MI
                }
                MIpath <- MIpath/(length(nodeindices)-1)
                pathMIs[pa] <- MIpath
                #print(MIpath)
              }#paths loop
              #print(pathMIs)
              #sumMI <- sumMI + max(pathMIs)
              NMvalues <- append(NMvalues, max(pathMIs))
            }#if length paths>0
          }#end hub loop
        }#end NM loop
        #sumMI <- sumMI/(length(NMindices)*length(hubindices))
        #maxMIbynet_NM[n] <- sumMI
      }#if NMs exist
      
      ##### NS #####
      if (length(NSindices)>0){
        #sumMI <- 0
        for (g in NSindices){
          print(paste0("Gene index ", g))
          paths <- all_simple_paths(graph, from=g, to=hubindices, mode="out")
          #get last element for every [[p]] paths entry, only retain those where this is index=h
          
          for (h in hubindices){
            if (g == h){next}
            lastelems <- which(sapply(paths, tail, 1) == h)
            pathstoh <- rep(list(NA), length(lastelems))
            #print(paste(length(lastelems), "paths to hub", h, "out of", length(paths), "paths to all hubs"))
            print("NS")
            if (length(pathstoh)> 0){
              
              for (el in 1:length(lastelems)){
                pathstoh[[el]] <- paths[[lastelems[el]]]
              }
              
              pathMIs <- rep(NA, length(pathstoh))
              print(paste(length(pathstoh), "paths from index", g, "to", h))
              #print(pathMIs)
              for (pa in 1:length(pathstoh)){ #p=list of gene names in path
                p <- pathstoh[[pa]]
                MIpath <- 0
                nodenames <- names(p)
                nrnodes <- length(nodenames)
                #which(colnames(adjmat) %in% nodenames)
                nodenames
                nodeindices <- rep(NA, length(nodenames))
                for (v in 1:length(nodeindices)){
                  nodeindices[v] <- which(colnames(adjmat) == nodenames[v])
                }# edge indices flowing from g->h
                for (e in 1:(length(nodeindices)-1)){
                  MI <- MutInfo(SBML, fi=nodeindices[e+1], j=nodeindices[e])
                  #sum these up along the path, avg by pathlength, then add to sumMI
                  MIpath <- MIpath + MI
                }
                MIpath <- MIpath/(length(nodeindices)-1)
                pathMIs[pa] <- MIpath
                #print(MIpath)
              }#paths loop
              #print(pathMIs)
              #sumMI <- sumMI + max(pathMIs)
              NSvalues <- append(NSvalues, max(pathMIs))
            }#if length paths>0
          }#end hub loop
        }#end NS loop
        #sumMI <- sumMI/(length(NSindices)*length(hubindices))
        #maxMIbynet_NS[n] <- sumMI
      }#if NS exist
      
    }#end loop if hubs exist in network
    #
    ##### save #####
    if (saveaftereverynet==TRUE){
      # saveRDS(maxMIbynet_PM, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_PM.RDS"))
      # saveRDS(maxMIbynet_NM, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_NM.RDS"))
      # saveRDS(maxMIbynet_NS, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_NS.RDS"))
      saveRDS(PMvalues, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_PM_bynode.RDS"))
      saveRDS(NMvalues, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_NM_bynode.RDS"))
      saveRDS(NSvalues, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_NS_bynode.RDS"))
    }
    
  }#end loop over nets
  if (saveResults==TRUE){
    #Save final results
    # saveRDS(maxMIbynet_PM, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_PM.RDS"))
    # saveRDS(maxMIbynet_NM, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_NM.RDS"))
    # saveRDS(maxMIbynet_NS, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_NS.RDS"))
    saveRDS(PMvalues, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_PM_bynode.RDS"))
    saveRDS(NMvalues, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_NM_bynode.RDS"))
    saveRDS(NSvalues, file=paste0(path, "Network measures/MaxMI/", HubDefinition, "_maxMIbynet_NS_bynode.RDS"))
  }
  # result[[1]] <- maxMIbynet_PM
  # result[[2]] <- maxMIbynet_NM
  # result[[3]] <- maxMIbynet_NS
  result[[1]] <- PMvalues
  result[[2]] <- NMvalues
  result[[3]] <- NSvalues
  return(result)
}

#Get number of nodes selected across all networks
selectionsize <- function(nets, path, measure="VBnDP", P=73){
  total <- 0
  for (n in 1:length(nets)){
    VB <- readRDS(paste0(path, "Network measures/net_",n,"_VB.RDS"))
    DP <- readRDS(paste0(path, "Network measures/net_",n,"_DP.RDS"))
    l <- length(VB)
    p <- ceiling(l*(P/100))
    sortedVB <- sort(VB, decreasing = T)
    sortedDP <- sort(DP, decreasing = T)
    selectedVB <- sortedVB[1:p]
    selectedDP <- sortedDP[1:p]
    VBuDP <- union(names(selectedVB), names(selectedDP))
    VBnDP <- intersect(names(selectedVB), names(selectedDP))
    
    if (measure=="VB"){total <- total + length(selectedVB)}
    else if (measure=="DP"){total <- total + length(selectedDP)}
    else if (measure=="VBuDP"){total <- total + length(VBuDP)}
    else if (measure=="VBnDP"){total <- total + length(VBnDP)}
  }
  return(total)
}

#Generate plot for Figure 1A: Sensitivity-specificity curves for a given static measure, averaged over dynamic measures
plotAvgDynSensSpec <- function(nets, path, measure="VBnDP", loaddata=TRUE){
  if (measure == "VB"){
    titlestr <- "(A) Average of dynamic measures - VB"
    xs <- 1
  } else if (measure == "DP"){
    titlestr <- "(B) Average of dynamic measures - DP"
    xs <- 2
  } else if (measure == "VBuDP"){
    titlestr <- paste0("(C) Average of dynamic measures - VB", expression("\u222A"), "DP")
    xs <- 3
  } else if (measure == "VBnDP"){
    titlestr <- paste0("(D) Average of dynamic measures - VB",expression("\u2229"), "DP")
    xs <- 4
  } else {print("Not a valid measure.")}
  if (loaddata==FALSE){
    print("Calculating sensitivity and specificity for all 3 dynamic measures across networks")
    DynAvgresults <- calculateSensSpec(statmeasure=measure, dynfunc="DynAvg", nets, path, saveResults=FALSE)
  }
  ysens_DynAvg <- yspec_DynAvg <- ysens_sd_DynAvg <- yspec_sd_DynAvg <- rep(0,100)
  for (n in 1:length(nets)){
    #add up all gain*(1/length_nets)
    if (loaddata==TRUE){
      DynAvg_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 4 ,".RDS"))
      DynAvg_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 4 ,".RDS"))
    } else {
      DynAvg_sens <- DynAvgresults[[n]]
      DynAvg_spec <- DynAvgresults[[length(nets)+n]]
    }
    
    ysens_DynAvg <- ysens_DynAvg %+% (DynAvg_sens*(1/length(nets)))
    yspec_DynAvg <- yspec_DynAvg %+% (DynAvg_spec*(1/length(nets)))
  }
  
  ysens_sd_DynAvg <- yspec_sd_DynAvg <- rep(0,100)
  ysens_DynAvg_bynet <- yspec_DynAvg_bynet <- rep(0, length(nets))
  
  for (T in 1:100){
    print(T)
    for (n in 1:length(nets)){
      if (loaddata==TRUE){
        DynAvg_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 4 ,".RDS"))
        DynAvg_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 4 ,".RDS"))
      } else {
        DynAvg_sens <- DynAvgresults[[n]]
        DynAvg_spec <- DynAvgresults[[length(nets)+n]]
      }
      
      #loop over T: collect sensgains at every T, calc sd over all length_nets
      ysens_DynAvg_bynet[n] <- DynAvg_sens[T]
      yspec_DynAvg_bynet[n] <- DynAvg_spec[T]
      
      ysens_sd_DynAvg[T] <- sd(ysens_DynAvg_bynet, na.rm=T)
      yspec_sd_DynAvg[T] <- sd(yspec_DynAvg_bynet, na.rm=T)
    }
  }
  
  ysens <- ysens_DynAvg
  yspec <- yspec_DynAvg
  ysens_sd <- ysens_sd_DynAvg
  yspec_sd <- yspec_sd_DynAvg
  
  x <- seq(1:100)
  df <- as.data.frame(cbind(x,ysens, yspec, ysens_sd, yspec_sd))
  plot(
    ggplot(df, aes(x, y = ysens, color = "")) + 
      geom_point(aes(y = ysens, col = "Sensitivity")) + labs(shape = "shape legend title", colour = "colour legend title") +
      geom_point(aes(y = yspec, col = "Specificity")) + 
      geom_ribbon(data=df,aes(ymin=ysens-ysens_sd,ymax=ysens+ysens_sd, fill='SD (Sensitivity)'),alpha=0.2, color=NA) + 
      geom_ribbon(data=df,aes(ymin=yspec-yspec_sd,ymax=yspec+yspec_sd, fill='SD (Specificity)'),alpha=0.2, color=NA) + 
      scale_x_continuous(labels=c(0,"",20,"",40,"",60,"",80,"",100), breaks=c(0,10,20,30,40,50,60,70,80,90,100)) + 
      scale_y_continuous(labels=c(0,"",0.2,"",0.4,"",0.6,"",0.8,"",1), breaks=c(0,0.1,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1)) + 
      coord_cartesian(ylim = c(0, 1), xlim =c(0,100)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme_grey(base_size = 25) + theme(legend.position = "none", legend.title=element_blank(),
                                         axis.text=element_text(size=17.5), axis.title=element_text(size=16), 
                                         plot.title = element_text(size=20)) +
      theme(panel.grid.major = element_line(colour="grey", size = (0.5)),
            panel.grid.minor = element_line(size = (0.0), colour="grey"), panel.background = element_blank()) + 
      theme(axis.line = element_line(color="black", size=1)) +
      labs(title=titlestr,
           x = expression(paste("Threshold percentage ", italic(T)," of genes labelled as 'high impact'")), 
           y = "Sensitivity & Specificity") +
      scale_colour_discrete("") + theme(plot.title = element_text(hjust = 0.5))
  )
  print(paste0("Selection threshold: T=", which(ysens>yspec)[1]))
  print(paste0("Value of sensitivity/specificity at selection threshold: ", round(ysens[which(ysens>yspec)[1]],3)))
}

#Generate plot for Figure 1B: Sensitivity-specificity curves for a given static measure, showing all three individual dynamic measures
plotIndivDynSensSpec <- function(nets, path, measure="VBnDP", loaddata=TRUE){
  # Load data: sensitivity and specificity for dynamic measures at every threshold
  ysens_gain <- yspec_gain <- ysens_sd_gain <- yspec_sd_gain <- rep(0,100)
  ysens_loss <- yspec_loss <- ysens_sd_loss <- yspec_sd_loss <- rep(0,100)
  ysens_hd <- yspec_hd <- ysens_sd_hd <- yspec_sd_hd <- rep(0,100)
  if (measure == "VB"){
    titlestr <- "(A) Individual dynamic measures - VB"
    xs <- 1
  } else if (measure == "DP"){
    titlestr <- "(B) Individual dynamic measures - DP"
    xs <- 2
  } else if (measure == "VBuDP"){
    titlestr <- paste0("(C) Individual dynamic measures - VB", expression("\u222A"), "DP")
    xs <- 3
  } else if (measure == "VBnDP"){
    titlestr <- paste0("(D) Individual dynamic measures - VB", expression("\u2229"),"DP")
    xs <- 4
  }
  
  if (loaddata==FALSE){
    print("Calculating sensitivity and specificity for all 3 dynamic measures across networks")
    Hammingresults <- calculateSensSpec(statmeasure=measure, dynfunc="Hamming", nets, path, saveResults=FALSE)
    AttrLossresults <- calculateSensSpec(statmeasure=measure, dynfunc="AttrLoss", nets, path, saveResults=FALSE)
    AttrGainresults <- calculateSensSpec(statmeasure=measure, dynfunc="AttrGain", nets, path, saveResults=FALSE)
  }
  
  for (n in 1:length(nets)){
    if (loaddata==TRUE){
      Hamming_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
      Hamming_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
      Loss_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
      Loss_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
      Gain_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
      Gain_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
    } else {
      Hamming_sens <- Hammingresults[[n]]
      Hamming_spec <- Hammingresults[[length(nets)+n]]
      Loss_sens <- AttrLossresults[[n]]
      Loss_spec <- AttrLossresults[[length(nets)+n]]
      Gain_sens <- AttrGainresults[[n]]
      Gain_spec <- AttrGainresults[[length(nets)+n]]
    }
    ysens_hd <- ysens_hd %+% (Hamming_sens*(1/length(nets)))
    yspec_hd <- yspec_hd %+% (Hamming_spec*(1/length(nets)))
    
    ysens_loss <- ysens_loss %+% (Loss_sens*(1/length(nets)))
    yspec_loss <- yspec_loss %+% (Loss_spec*(1/length(nets)))
    
    ysens_gain <- ysens_gain %+% (Gain_sens*(1/length(nets)))
    yspec_gain <- yspec_gain %+% (Gain_spec*(1/length(nets)))
  }
  ysens_gain
  yspec_gain
  
  ysens_sd_gain <- yspec_sd_gain <- ysens_sd_loss <- yspec_sd_loss <- ysens_sd_hd <- yspec_sd_hd <- rep(0,100)
  ysens_gain_bynet <- yspec_gain_bynet <- ysens_loss_bynet <- yspec_loss_bynet <- ysens_hd_bynet <- yspec_hd_bynet <- rep(0, length(nets))
  for (T in 1:100){
    print(T)
    for (n in 1:length(nets)){
      if (loaddata==TRUE){
        Hamming_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
        Hamming_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
        Loss_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
        Loss_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
        Gain_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
        Gain_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
      } else {
        Hamming_sens <- Hammingresults[[n]]
        Hamming_spec <- Hammingresults[[length(nets)+n]]
        Loss_sens <- AttrLossresults[[n]]
        Loss_spec <- AttrLossresults[[length(nets)+n]]
        Gain_sens <- AttrGainresults[[n]]
        Gain_spec <- AttrGainresults[[length(nets)+n]]
      }
      #loop over T: collect sensgains at every T, calc sd over all length_nets
      ysens_hd_bynet[n] <- Hamming_sens[T]
      yspec_hd_bynet[n] <- Hamming_spec[T]
      ysens_loss_bynet[n] <- Loss_sens[T]
      yspec_loss_bynet[n] <- Loss_spec[T]
      ysens_gain_bynet[n] <- Gain_sens[T]
      yspec_gain_bynet[n] <- Gain_spec[T]
      
      ysens_sd_hd[T] <- sd(ysens_hd_bynet, na.rm=T)
      yspec_sd_hd[T] <- sd(yspec_hd_bynet, na.rm=T)
      ysens_sd_loss[T] <- sd(ysens_loss_bynet, na.rm=T)
      yspec_sd_loss[T] <- sd(yspec_loss_bynet, na.rm=T)
      ysens_sd_gain[T] <- sd(ysens_gain_bynet, na.rm=T)
      yspec_sd_gain[T] <- sd(yspec_gain_bynet, na.rm=T)
    }
  }
  
  x <- seq(1:100)
  df <- as.data.frame(cbind(x, ysens_hd, yspec_hd, ysens_loss, yspec_loss, ysens_gain, yspec_gain))
  
  plot(
    ggplot(df, aes(x, y = ysens, color = "")) + 
      #gain
      geom_point(aes(y = ysens_gain, col = "Sensitivity"), col="#237d02", alpha=1, shape=19) +
      geom_point(aes(y = yspec_gain, col = "Specificity"), col="#237d02", alpha=1, shape=17) +
      #loss
      geom_point(aes(y = ysens_loss, col = "Sensitivity"), col="#ab7105", alpha=1, shape=19) +
      geom_point(aes(y = yspec_loss, col = "Specificity"), col="#ab7105", alpha=1, shape=17) +
      #Hamming
      geom_point(aes(y = ysens_hd, col = "Sensitivity"), col="purple", alpha=1, shape=19) +
      geom_point(aes(y = yspec_hd, col = "Specificity"), col="purple", alpha=1, shape=17) +
      
      scale_x_continuous(labels=c(0,"",20,"",40,"",60,"",80,"",100), breaks=c(0,10,20,30,40,50,60,70,80,90,100)) + 
      scale_y_continuous(labels=c(0,"",0.2,"",0.4,"",0.6,"",0.8,"",1), breaks=c(0,0.1,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1)) + 
      coord_cartesian(ylim = c(0, 1), xlim =c(0,100)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme_grey(base_size = 25) + theme(legend.position = "bottom", 
                                         axis.text=element_text(size=17.5), axis.title=element_text(size=16), 
                                         plot.title = element_text(size=20)) + 
      theme(panel.grid.major = element_line(colour="grey", size = (0.5)),
            panel.grid.minor = element_line(size = (0.0), colour="grey"), panel.background = element_blank()) + 
      theme(axis.line = element_line(color="black", size=1)) +
      labs(title=titlestr,
           x = "Threshold percentage T of genes labelled as 'high impact'", y = "Sensitivity & Specificity") + 
      scale_colour_discrete("") + theme(plot.title = element_text(hjust = 0.5))
  )
  
  Hammingthreshold <- which(ysens_hd>yspec_hd)[1]
  Lossthreshold <- which(ysens_loss>yspec_loss)[1]
  Gainthreshold <- which(ysens_gain>yspec_gain)[1]
  print("Thresholds for Hamming distance, Attractor loss and Attractor gain:")
  print(c(Hammingthreshold, Lossthreshold, Gainthreshold))
}

#Plot figures comparing network properties (Fig3)
plotFig3ABC <- function(nets, path, Figure, threshold=73, loaddata=TRUE){
  if (Figure %in% c("A","B","C") == FALSE){return("Not a valid figure, choose A-C.")}
  l <- length(nets)
  toppercent <- threshold
  Havgperc <- nonHavgperc <- PMavgperc <- NMavgperc <- NSavgperc <- rep(NA, l)
  Hcount <- nonHcount <- PMcount <- NMcount <- NScount <- 0
  Hvalues <- PMvalues <- NMvalues <- NSvalues <- c()
  
  if (loaddata==FALSE){
    VBresults <- calculateAllMeasures(nets, path, measure2calc="VB", saveResults=FALSE)
    DPresults <- calculateAllMeasures(nets, path, measure2calc="DP", saveResults=FALSE)
    Zresults <- calculateAllMeasures(nets, path, measure2calc="Z", saveResults=FALSE)
    Hammingresults <- calculateAllMeasures(nets, path, measure2calc="Hamming", saveResults=FALSE)
    AttrLossresults <- calculateAllMeasures(nets, path, measure2calc="AttrLoss", saveResults=FALSE)
    AttrGainresults <- calculateAllMeasures(nets, path, measure2calc="AttrGain", saveResults=FALSE)
  }
  
  for (n in 1:length(nets)){
    print(n)
    ##### get measures #####
    # get percentiles, VBnDP rankings etc 
    netfile <- paste(path, nets[n], sep = "")
    SBML <- loadNetwork(netfile)
    adjmat <- conv2adjmat(SBML, inputcorrected = F)
    graph <- graph_from_adjacency_matrix(adjmat)
    size <- dim(adjmat)[1]
    p <- ceiling((toppercent/100)*size)
    if (loaddata==TRUE){
      VB <- readRDS(paste0(path, "Network measures/net_", n, "_VB.RDS"))
      DP <- readRDS(paste0(path, "Network measures/net_", n, "_DP.RDS"))
      Z <- readRDS(paste0(path, "Network measures/net_", n, "_Z.RDS"))
      Hamming <- readRDS(paste0(path, "Network measures/net_", n, "_Hamming.RDS")) #Hamming distance
      AttrLoss <- readRDS(paste0(path, "Network measures/net_", n, "_AttrLoss.RDS")) #Attractor loss
      AttrGain <- readRDS(paste0(path, "Network measures/net_", n, "_AttrGain.RDS")) #Attractor gain
    } else {
      VB <- VBresults[[n]]
      DP <- DPresults[[n]]
      Z <- Zresults[[n]]
      Hamming <- Hammingresults[[n]]
      AttrLoss <- AttrLossresults[[n]]
      AttrGain <- AttrGainresults[[n]]
    }
    AvgDynRank <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
    xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
    xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
    
    VBRank <- rank(-VB, ties.method="min")
    DPRank <- rank(-DP, ties.method="min")
    name <- intersect(names(xselect_VB), names(xselect_DP))
    all_VBnDP <- colMeans(rbind(VBRank, DPRank))
    
    VBnDP <- rep(NA, length(name))
    names(VBnDP) <- name
    Z_inVBnDP <- Z[which(names(Z) %in% names(VBnDP))]
    for (g in 1:length(VBnDP)){
      RankA <- VBRank[which(names(VBRank) == names(VBnDP)[g])]
      RankB <- DPRank[which(names(DPRank) == names(VBnDP)[g])]
      VBnDP[g] <- mean(c(RankA, RankB))
    }
    ZmatchVBnDP <- Z_inVBnDP[match(names(VBnDP), names(Z_inVBnDP))]
    StatZQuantRatio <- rep(NA, length(VBnDP))
    names(StatZQuantRatio) <- names(VBnDP)
    for (g in 1:length(VBnDP)){
      StatZQuantRatio[g] <- percentile(VBnDP, VBnDP[g], fullrange=T) - percentile(-ZmatchVBnDP, -ZmatchVBnDP[g], fullrange=T) 
      #low number=high rank
    }
    gatekeepers <- which(StatZQuantRatio > 0)
    
    ##### define sets to compare for given property #####
    NSindices <- which(colnames(adjmat) %in% setdiff(names(AvgDynRank), names(VBnDP)))
    hubindices <- which(Z >= 2.5) #Z>=2.5 hub definition by Guimera et al.
    nonhubindices <- setdiff(seq(1:length(SBML$genes)), hubindices)
    PMindices <- which(colnames(adjmat) %in% setdiff(names(gatekeepers), colnames(adjmat)[c(NSindices)]))
    NMindices <- setdiff(seq(1:length(SBML$genes)), c(NSindices, PMindices))
    
    Hcount <- Hcount + length(hubindices)
    nonHcount <- nonHcount + length(nonhubindices)
    PMcount <- PMcount + length(PMindices)
    NMcount <- NMcount + length(NMindices)
    NScount <- NScount + length(NSindices)
    
    ##### write to file #####
    if (Figure == "A"){
      #resvec <- sort(VBnDP)
      resvec <- sort(all_VBnDP)
    } else if (Figure %in% c("B","C")){
      resvec <- sort(AvgDynRank)
    }
    labelvec <- rep(NA, length(resvec))
    #get average percentile position of Hs, PMs, NMs in resvec -> print to three avgvecs
    Havgs <- nonHavgs <- PMavgs <- NMavgs <- NSavgs <- c()
    for (g in 1:length(resvec)){
      index <- which(colnames(adjmat) == names(resvec)[g])
      if (index %in% hubindices){
        labelvec[g] <- "Hub"
        if (Figure %in% c("A", "B")){
          Havgs <- append(Havgs, percentile(resvec, resvec[g]))
        } else if (Figure == "C"){
          Havgs <- append(Havgs, Z[index])
        }
      }#if in hubindices
      if (index %in% nonhubindices){
        labelvec[g] <- "Non-Hub"
        if (Figure %in% c("A", "B")){
          nonHavgs <- append(nonHavgs, percentile(resvec, resvec[g]))
        } else if (Figure == "C"){
          nonHavgs <- append(nonHavgs, Z[index])
        }
      }#if in nonhubindices
      if (index %in% PMindices){
        labelvec[g] <- "PM"
        if (Figure %in% c("A", "B")){
          PMavgs <- append(PMavgs, percentile(resvec, resvec[g]))
        } else if (Figure == "C"){
          PMavgs <- append(PMavgs, Z[index])
        }
      }
      if (index %in% NMindices){
        labelvec[g] <- "NM"
        if (Figure %in% c("A", "B")){
          NMavgs <- append(NMavgs, percentile(resvec, resvec[g]))
        } else if (Figure == "C"){
          NMavgs <- append(NMavgs, Z[index])
        }
      }
      if (index %in% NSindices){
        labelvec[g] <- "NS"
        if (Figure %in% c("A", "B")){
          NSavgs <- append(NSavgs, percentile(resvec, resvec[g]))
        } else if (Figure == "C"){
          NSavgs <- append(NSavgs, Z[index])
        }
      }
    }
    # Havgperc[n] <- mean(Havgs)
    # nonHavgperc[n] <- mean(nonHavgs)
    # PMavgperc[n] <- mean(PMavgs)
    # NMavgperc[n] <- mean(NMavgs)
    # NSavgperc[n] <- mean(NSavgs)
    
    Hvalues <- append(Hvalues, Havgs)
    PMvalues <- append(PMvalues, PMavgs)
    NMvalues <- append(NMvalues, NMavgs)
    NSvalues <- append(NSvalues, NSavgs)
  }
  
  if (Figure == "A"){
    l <- (Hcount+PMcount+NMcount+NScount)
    df <- matrix(0, nrow = l, ncol=2)
    df[1:Hcount,1] <- Hvalues
    df[(Hcount+1):(Hcount+PMcount),1] <- PMvalues
    df[(Hcount+PMcount+1):(Hcount+PMcount+NMcount),1] <- NMvalues
    df[(Hcount+PMcount+NMcount+1):l,1] <- NSvalues
    df[1:Hcount,2] <- "Hub"
    df[(Hcount+1):(Hcount+PMcount),2] <- "PM \n(Gatekeepers)"
    df[(Hcount+PMcount+1):(Hcount+PMcount+NMcount),2] <- "NM"
    df[(Hcount+PMcount+NMcount+1):l,2] <- "NS"
    df <- as.data.frame(df)
    df
    
    names(df) <- c("percentiles", "group")
    df$percentiles <- as.numeric(as.character(df$percentiles))
    df$group <- factor(df$group, levels = c("Hub", "PM \n(Gatekeepers)", "NM", "NS"))
    
    padjvals <- compare_means(percentiles ~ group, df, method="wilcox.test", paired=F, p.adjust.method = "bonferroni") #get p.adj
    #padjvals <- padjvals[-c(2),] #remove some comparisons for clean plot
    maxlimitindex <- which(padjvals$p.adj == 1)
    minlimitindex <- which(padjvals$p.adj < 0.00001)
    for (p in maxlimitindex){
      padjvals$p.adj[p] <- ">0.99"
    }
    for (p in minlimitindex){
      padjvals$p.adj[p] <- "<1e-05"
    }
    padjvals <- padjvals[-c(2,3),] #remove some comparisons for clean plot
    plabel <- c(1.15,1.3,1.5,1.7)
    plot(
      ggboxplot(df, x = "group" , y = "percentiles", xlab="Classification", ylab="Percentile score in static impact ranking",
                font.x=20, font.y=20, font.subtitle=20, font.caption=20, font.tickslab=20
      ) +
        stat_pvalue_manual(padjvals, label = "p.adj", y.position = plabel, size=7) + ggtitle("A") +
        theme(plot.title = element_text(size = 25, face = "bold")) +
        annotate("text", x=c(1,2,3,4), y=-0.1, size=7,
                 label=c(paste0("n=", Hcount),paste0("n=", PMcount),paste0("n=", NMcount),paste0("n=", NScount)))
    )
    
    # boxplot(Hvalues, PMvalues, NMvalues, xlab="Classification", ylab="Percentile score in static impact ranking",
    #         names = c("Hub", "PM/Gatekeepers", "NM"))
  }#Fig A
  
  if (Figure == "B"){
    l <- (Hcount+PMcount+NMcount+NScount)
    df <- matrix(0, nrow = l, ncol=2)
    df[1:Hcount,1] <- Hvalues
    df[(Hcount+1):(Hcount+PMcount),1] <- PMvalues
    df[(Hcount+PMcount+1):(Hcount+PMcount+NMcount),1] <- NMvalues
    df[(Hcount+PMcount+NMcount+1):l,1] <- NSvalues
    df[1:Hcount,2] <- "Hub"
    df[(Hcount+1):(Hcount+PMcount),2] <- "PM \n(Gatekeepers)"
    df[(Hcount+PMcount+1):(Hcount+PMcount+NMcount),2] <- "NM"
    df[(Hcount+PMcount+NMcount+1):l,2] <- "NS"
    df <- as.data.frame(df)
    df
    
    # #plabel <- c(1.15,1.3,1.5,1.7,1.9,2.1)
    plabel <- c(1.15,1.3,1.5,1.7)
    
    names(df) <- c("percentiles", "group")
    df$percentiles <- as.numeric(as.character(df$percentiles))
    df$group <- factor(df$group, levels = c("Hub", "PM \n(Gatekeepers)", "NM", "NS"))
    padjvals <- compare_means(percentiles ~ group, df, method="wilcox.test", paired=F, p.adjust.method = "bonferroni") #get p.adj
    maxlimitindex <- which(padjvals$p.adj == 1)
    minlimitindex <- which(padjvals$p.adj < 0.00001)
    for (p in maxlimitindex){
      padjvals$p.adj[p] <- ">0.99"
    }
    for (p in minlimitindex){
      padjvals$p.adj[p] <- "<1e-05"
    }
    padjvals <- padjvals[-c(2,3),] #remove some comparisons for clean plot
    plot(
      ggboxplot(df, x = "group" , y = "percentiles", xlab="Classification", ylab="Percentile score in dynamic impact ranking",
                font.x=20, font.y=20, font.subtitle=20, font.caption=20, font.tickslab=20
      ) +
        stat_pvalue_manual(padjvals, label = "p.adj", y.position = plabel, size=7) + ggtitle("B") +
        theme(plot.title = element_text(size = 25, face = "bold")) +
        annotate("text", x=c(1,2,3,4), y=-0.25, size=7,
                 label=c(paste0("n=", Hcount),paste0("n=", PMcount),paste0("n=", NMcount),paste0("n=", NScount)))
    )
  }#Fig B
  
  
  if (Figure == "C"){
    l <- (Hcount+PMcount+NMcount+NScount)
    df <- matrix(0, nrow = l, ncol=2)
    df[1:Hcount,1] <- Hvalues
    df[(Hcount+1):(Hcount+PMcount),1] <- PMvalues
    df[(Hcount+PMcount+1):(Hcount+PMcount+NMcount),1] <- NMvalues
    df[(Hcount+PMcount+NMcount+1):l,1] <- NSvalues
    df[1:Hcount,2] <- "Hub"
    df[(Hcount+1):(Hcount+PMcount),2] <- "PM \n(Gatekeepers)"
    df[(Hcount+PMcount+1):(Hcount+PMcount+NMcount),2] <- "NM"
    df[(Hcount+PMcount+NMcount+1):l,2] <- "NS"
    df <- as.data.frame(df)
    df
    
    #plabel <- c(4.2,4.6,5.0,1.8,2.1,2.6)
    plabel <- c(5.5,4.5,5.0,5.8)
    
    names(df) <- c("percentiles", "group")
    df$percentiles <- as.numeric(as.character(df$percentiles))
    df$group <- factor(df$group, levels = c("Hub", "PM \n(Gatekeepers)", "NM", "NS"))
    padjvals <- compare_means(percentiles ~ group, df, method="wilcox.test", paired=F, p.adjust.method = "bonferroni") #get p.adj
    maxlimitindex <- which(padjvals$p.adj == 1)
    minlimitindex <- which(padjvals$p.adj < 0.00001)
    for (p in maxlimitindex){
      padjvals$p.adj[p] <- ">0.99"
    }
    for (p in minlimitindex){
      padjvals$p.adj[p] <- "<1e-05"
    }
    padjvals <- padjvals[-c(2,3),] #remove some comparisons for clean plot
    plot(
      ggboxplot(df, x = "group" , y = "percentiles", xlab="Classification", ylab="Connectivity (z-score)",
                font.x=20, font.y=20, font.subtitle=20, font.caption=20, font.tickslab=20
      ) +
        stat_pvalue_manual(padjvals, label = "p.adj", y.position = plabel, size=7) + ggtitle("C") +
        theme(plot.title = element_text(size = 25, face = "bold")) +
        annotate("text", x=c(1,2,3,4), y=-2, size=7,
                 label=c(paste0("n=", Hcount),paste0("n=", PMcount),paste0("n=", NMcount),paste0("n=", NScount)))
    )
    # boxplot(Hvalues, PMvalues, NMvalues, NSvalues, xlab="Classification", ylab="Connectivity (z-score)",
    #         names = c("Hub", "PM/Gatekeepers", "NM", "NS"))
  }#Fig C
  
}

#Plot figure 3D, maximal MI along paths to hubs
plotFig3D <- function(nets, path, threshold=73, loaddata=TRUE, calculatepaths=FALSE, saveaftereverynet=FALSE, saveResults=TRUE){
  l <- length(nets)
  toppercent <- threshold
  Havgperc <- PMavgperc <- NMavgperc <- NSavgperc <- rep(NA, length(nets))
  PMcount <- NMcount <- NScount <- 0
  
  if (loaddata==FALSE){
    VBresults <- calculateAllMeasures(nets, path, measure2calc="VB", saveResults=FALSE)
    DPresults <- calculateAllMeasures(nets, path, measure2calc="DP", saveResults=FALSE)
    Zresults <- calculateAllMeasures(nets, path, measure2calc="Z", saveResults=FALSE)
    Hammingresults <- calculateAllMeasures(nets, path, measure2calc="Hamming", saveResults=FALSE)
    AttrLossresults <- calculateAllMeasures(nets, path, measure2calc="AttrLoss", saveResults=FALSE)
    AttrGainresults <- calculateAllMeasures(nets, path, measure2calc="AttrGain", saveResults=FALSE)
  }
  
  #Count nodes across nets in which there are hubs
  for (n in 1:length(nets)){
    print(n)
    ##### get measures #####
    # get percentiles, VBnDP rankings etc 
    netfile <- paste(path, nets[n], sep = "")
    SBML <- loadNetwork(netfile)
    adjmat <- conv2adjmat(SBML, inputcorrected = F)
    graph <- graph_from_adjacency_matrix(adjmat)
    size <- dim(adjmat)[1]
    p <- ceiling((toppercent/100)*size)
    if (loaddata==TRUE){
      VB <- readRDS(paste0(path, "Network measures/net_", n, "_VB.RDS"))
      DP <- readRDS(paste0(path, "Network measures/net_", n, "_DP.RDS"))
      Z <- readRDS(paste0(path, "Network measures/net_", n, "_Z.RDS"))
      Hamming <- readRDS(paste0(path, "Network measures/net_", n, "_Hamming.RDS")) #Hamming distance
      AttrLoss <- readRDS(paste0(path, "Network measures/net_", n, "_AttrLoss.RDS")) #Attractor loss
      AttrGain <- readRDS(paste0(path, "Network measures/net_", n, "_AttrGain.RDS")) #Attractor gain
    } else {
      VB <- VBresults[[n]]
      DP <- DPresults[[n]]
      Z <- Zresults[[n]]
      Hamming <- Hammingresults[[n]]
      AttrLoss <- AttrLoss[[n]]
      AttrGain <- AttrGain[[n]]
    }
    AvgDynRank <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
    xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
    xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
    
    VBRank <- rank(-VB, ties.method="min")
    DPRank <- rank(-DP, ties.method="min")
    name <- intersect(names(xselect_VB), names(xselect_DP))
    
    VBnDP <- rep(NA, length(name))
    names(VBnDP) <- name
    Z_inVBnDP <- Z[which(names(Z) %in% names(VBnDP))]
    for (g in 1:length(VBnDP)){
      RankA <- VBRank[which(names(VBRank) == names(VBnDP)[g])]
      RankB <- DPRank[which(names(DPRank) == names(VBnDP)[g])]
      VBnDP[g] <- mean(c(RankA, RankB))
    }
    ZmatchVBnDP <- Z_inVBnDP[match(names(VBnDP), names(Z_inVBnDP))]
    StatZQuantRatio <- rep(NA, length(VBnDP))
    names(StatZQuantRatio) <- names(VBnDP)
    for (g in 1:length(VBnDP)){
      StatZQuantRatio[g] <- percentile(VBnDP, VBnDP[g], fullrange=T) - percentile(-ZmatchVBnDP, -ZmatchVBnDP[g], fullrange=T) 
      #low number=high rank
    }
    gatekeepers <- which(StatZQuantRatio > 0)
    
    ##### define sets to compare for given property #####
    NSindices <- which(colnames(adjmat) %in% setdiff(names(AvgDynRank), names(VBnDP)))
    hubindices <- which(Z >= 2.5) #Z>=2.5 hub definition by Guimera et al.
    nonhubindices <- setdiff(seq(1:length(SBML$genes)), hubindices)
    PMindices <- which(colnames(adjmat) %in% setdiff(names(gatekeepers), colnames(adjmat)[c(NSindices)]))
    NMindices <- setdiff(seq(1:length(SBML$genes)), c(NSindices, PMindices))
    if (length(hubindices) > 0){
      PMcount <- PMcount + length(PMindices)
      NMcount <- NMcount + length(NMindices)
      NScount <- NScount + length(NSindices)
    }
  }
  
  if (calculatepaths==FALSE){
    #Load results if they have already been calculated previously
    # PMavgperc <- readRDS(paste0(path, "Network measures/MaxMI/Guimera_maxMIbynet_PM.RDS"))
    # NMavgperc <- readRDS(paste0(path, "Network measures/MaxMI/Guimera_maxMIbynet_NM.RDS"))
    # NSavgperc <- readRDS(paste0(path, "Network measures/MaxMI/Guimera_maxMIbynet_NS.RDS"))
    PMavgperc <- readRDS(paste0(path, "Network measures/MaxMI/Guimera_maxMIbynet_PM_bynode.RDS"))
    NMavgperc <- readRDS(paste0(path, "Network measures/MaxMI/Guimera_maxMIbynet_NM_bynode.RDS"))
    NSavgperc <- readRDS(paste0(path, "Network measures/MaxMI/Guimera_maxMIbynet_NS_bynode.RDS"))
  } else {
    #Otherwise, calculate them:
    maxMIs <- calcMaxMItoHubs(nets, path, saveaftereverynet = saveaftereverynet, saveResults = saveResults)
    PMavgperc <- maxMIs[[1]]
    NMavgperc <- maxMIs[[2]]
    NSavgperc <- maxMIs[[3]]
  }
  PMcount <- length(PMavgperc)
  NMcount <- length(NMavgperc)
  NScount <- length(NSavgperc)
  l <- PMcount + NMcount + NScount
  
  df <- matrix(0, nrow = l, ncol=2)
  df[1:PMcount,2] <- "PM \n(Gatekeepers)"
  df[(PMcount+1):(PMcount+NMcount),2] <- "NM"
  df[(PMcount+NMcount+1):(l),2] <- "NS"
  df[1:PMcount,1] <- PMavgperc
  df[(PMcount+1):(PMcount+NMcount),1] <- NMavgperc
  df[(PMcount+NMcount+1):(l),1] <- NSavgperc
  df <- as.data.frame(df)
  df
  names(df) <- c("percentiles", "group")
  df$percentiles <- as.numeric(as.character(df$percentiles))
  df$group <- factor(df$group, levels = c("PM \n(Gatekeepers)", "NM", "NS")) 
  
  padjvals <- compare_means(percentiles ~ group, df, method="wilcox.test", paired=F, p.adjust.method = "bonferroni") #get p.adj
  
  maxlimitindex <- which(padjvals$p.adj == 1)
  minlimitindex <- which(padjvals$p.adj < 0.00001)
  for (p in maxlimitindex){
    padjvals$p.adj[p] <- ">0.99"
  }
  for (p in minlimitindex){
    padjvals$p.adj[p] <- "<1e-05"
  }
  
  plot(
    ggboxplot(df, x = "group" , y = "percentiles", xlab="Classification", ylab="", 
              font.x=20, font.y=20, font.subtitle=20, font.caption=20, font.tickslab=20) + 
      stat_pvalue_manual(padjvals, label = "p.adj", y.position = c(0.85,0.95,1.05), size=7) +
      ylab("Averaged maximal MI along paths to hubs") + 
      theme(plot.title = element_text(size = 25, face = "bold")) + ggtitle("   D") + 
      annotate("text", x=c(1,2,3), y=-0.1, size=7,
               label=c(paste0("P=", PMcount),paste0("P=", NMcount),paste0("P=", NScount)))
  )
  
}

#Generate bootstrap data from seed, plot threshold robustness
plotThresholdBootstrap <- function(nets, path, reps=10000, generateNewSeed=FALSE, loaddata=TRUE, 
                                   saveResults=FALSE){
  netseq <- seq(1:length(nets))
  if (generateNewSeed==T){
    #If seed does not yet exist, generate & save it
    seedmatrix <- matrix(NA, nrow = reps, ncol = length(nets))
    for (r in 1:reps){
      randomselection <- sample(netseq, size=length(nets), replace = T)
      seedmatrix[r,] <- randomselection
    }
    seed <- seedmatrix
    saveRDS(seedmatrix, file=paste0(pathtoscripts, "bootstrap/seed.RDS"))
  } else {
    #If seed is already available for chosen networks, load here:
    seed <- readRDS(paste0(pathtoscripts, "bootstrap/seed.RDS"))
  }
  
  nprime <- length(nets)
  replacing = T #networks can be selected multiple times, nprime=n=34, replacing=T, select=F for bootstrap
  select = F #== F if selecting, ==T if skipping the randomly chosen nets
  gainT <- lossT <- hammingT <- balancedT <- rep(NA,reps)
  saving <- F #set TRUE if intermediate results should be saved
  for (rand in 1:reps){
    print(paste0(rand, "/", reps))
    randvec <- seed[rand,]
    gain_avg_sens <- gain_avg_spec <- loss_avg_sens <- loss_avg_spec <- hamming_avg_sens <- hamming_avg_spec <- rep(0, 100)
    DynAvg_avg_sens <- DynAvg_avg_spec <- rep(0,100)
    if (loaddata == FALSE){#need to generate data for Sensitivity and specificity first
      HammingSensSpec <- calculateSensSpec(statmeasure = "VBnDP", dynfunc="Hamming", nets, path)
      AttrLossSensSpec <- calculateSensSpec(statmeasure = "VBnDP", dynfunc="AttrLoss", nets, path)
      AttrGainSensSpec <- calculateSensSpec(statmeasure = "VBnDP", dynfunc="AttrGain", nets, path)
      DynAvgSensSpec <- calculateSensSpec(statmeasure = "VBnDP", dynfunc="DynAvg", nets, path)
    }
    for (n in randvec){
      #ys: 1=Hamming, 2=AttrLoss, 3=AttrGain, 4=DynAvg
      if (loaddata == TRUE){
        dynavg_netn_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net",n, "_xs4_ys4.RDS"))
        dynavg_netn_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net",n, "_xs4_ys4.RDS"))
        gain_netn_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net",n, "_xs4_ys3.RDS"))
        gain_netn_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net",n, "_xs4_ys3.RDS"))
        loss_netn_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net",n, "_xs4_ys2.RDS"))
        loss_netn_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net",n, "_xs4_ys2.RDS"))
        hamming_netn_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net",n, "_xs4_ys1.RDS"))
        hamming_netn_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net",n, "_xs4_ys1.RDS"))
      } else {
        dynavg_netn_sens <- DynAvgSensSpec[[n]]
        dynavg_netn_spec <- DynAvgSensSpec[[length(nets)+n]]
        gain_netn_sens <- AttrGainSensSpec[[n]]
        gain_netn_spec <- AttrGainSensSpec[[length(nets)+n]]
        loss_netn_sens <- AttrLossSensSpec[[n]]
        loss_netn_spec <- AttrLossSensSpec[[length(nets)+n]]
        hamming_netn_sens <- HammingSensSpec[[n]]
        hamming_netn_spec <- HammingSensSpec[[length(nets)+n]]
      }
      DynAvg_avg_sens <- DynAvg_avg_sens %+% dynavg_netn_sens
      DynAvg_avg_spec <- DynAvg_avg_spec %+% dynavg_netn_spec
      
      gain_avg_sens <- gain_avg_sens %+% gain_netn_sens
      gain_avg_spec <- gain_avg_spec %+% gain_netn_spec
      
      loss_avg_sens <- loss_avg_sens %+% loss_netn_sens
      loss_avg_spec <- loss_avg_spec %+% loss_netn_spec
      
      hamming_avg_sens <- hamming_avg_sens %+% hamming_netn_sens
      hamming_avg_spec <- hamming_avg_spec %+% hamming_netn_spec
      
    }
    gainT[rand] <- which(gain_avg_sens > gain_avg_spec)[1]
    lossT[rand] <- which(loss_avg_sens > loss_avg_spec)[1]
    hammingT[rand] <- which(hamming_avg_sens > hamming_avg_spec)[1]
    balancedT[rand] <- which(DynAvg_avg_sens > DynAvg_avg_spec)[1]
    #balancedT[rand] <- mean(c(gainT[rand], lossT[rand], hammingT[rand]))
    if (saveResults == T | (rand==reps)){
      saveRDS(gainT, file=paste0(pathtoscripts, "bootstrap/gainT.RDS"))
      saveRDS(lossT, file=paste0(pathtoscripts, "bootstrap/lossT.RDS"))
      saveRDS(hammingT, file=paste0(pathtoscripts, "bootstrap/hammingT.RDS"))
      saveRDS(balancedT, file=paste0(pathtoscripts, "bootstrap/balancedT.RDS"))
    }
  }
  
  df <- as.matrix(c(gainT,lossT,hammingT,balancedT))
  groups <- c(rep("gainT", reps), rep("lossT", reps), rep("hammingT", reps), rep("balancedT", reps))
  df <- as.data.frame(cbind(df, groups))
  names(df) <- c("Thresholds", "groups")
  
  # gainT <- readRDS(paste0(pathtoscripts, "bootstrap/gainT.RDS"))
  # lossT <- readRDS(paste0(pathtoscripts, "bootstrap/lossT.RDS"))
  # hammingT <- readRDS(paste0(pathtoscripts, "bootstrap/hammingT.RDS"))
  # balancedT <- readRDS(paste0(pathtoscripts, "bootstrap/balancedT.RDS"))
  boxplot(gainT, lossT, hammingT, balancedT, names=c("Attractor \ngain", "Attractor \nloss", "Hamming \ndistance", "Dynamic \naverage"),
          main="Bootstrap analysis of selection thresholds", cex.axis=0.8)
  print(median(balancedT))
  median(IQR(balancedT))
}

#Generate z-scores from a given dataset of static and dynamic measures, providing indices of the test set of nodes that is to be excluded
normaliseDataset <- function(Dataset, test.ids){
  #Throw out test data
  traindataset <- Dataset$data[,-test.ids]
  #Turn all 8 static measures in $data matrix into z-scores, normed over training set
  for (r in 1:dim(Dataset$data)[1]){
    Dataset$data[r,] <- (Dataset$data[r,]-mean(traindataset[r,], na.rm=T))/sd(traindataset[r,], na.rm=T)
    #Transforms all rows into z-scores (including the test set), based on normalisation factors from the training set
  }
  #Turn the three individual dynamic measures into z-scores as well, normalised using trainingset
  HammingMean <- mean(Dataset$hamming[-test.ids],na.rm=T)
  HammingSD <- sd(Dataset$hamming[-test.ids],na.rm=T)
  LossMean <- mean(Dataset$loss[-test.ids],na.rm=T)
  LossSD <- sd(Dataset$loss[-test.ids],na.rm=T)
  GainMean <- mean(Dataset$gain[-test.ids],na.rm=T)
  GainSD <- sd(Dataset$gain[-test.ids],na.rm=T)
  Dataset$hamming <- (Dataset$hamming - HammingMean)/HammingSD
  Dataset$loss <- (Dataset$loss - LossMean)/LossSD
  Dataset$gain <- (Dataset$gain - GainMean)/GainSD
  
  #Overwrite $dynavg as being the mean of three z-scores
  Dataset$dynavg <- colMeans(rbind(Dataset$hamming, Dataset$loss, Dataset$gain))
  
  #Criterion for labelling as high or low impact, save this for all folds and runs! -> threshold boxplot
  #Overwrite old threshold
  Dataset$threshold <- median(Dataset$dynavg, na.rm=T)
  labs <- rep(0, dim(Dataset$data)[2])
  oldLabs <- rep("low_impact", dim(Dataset$data)[2])
  
  for (g in 1:dim(Dataset$data)[2]){
    if (Dataset$dynavg[g] > Dataset$threshold){
      labs[g] <- 1
      oldLabs[g] <- "high_impact"
    }
  }
  Dataset$labs <- labs
  Dataset$oldLabs <- oldLabs
  
  return(Dataset)
}


