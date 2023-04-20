#SIM
create_SIM_adjmat <- function(n){
  #creates SIM adjmat with n outputs
  adjmat <- matrix(0, n+1,n+1)
  adjmat[1,] <- 1
  return(adjmat)
} 
SIM_adjmat <- create_SIM_adjmat(3)
SIM_graph <- graph_from_adjacency_matrix(SIM_adjmat)

#Multi-output FFL
create_MOFFL_adjmat <- function(n){
  #creates MOFFL adjmat with n outputs
  adjmat <- matrix(0, n+2,n+2)
  adjmat[1,2] <- 1
  adjmat[1:2,3:dim(adjmat)[2]] <- 1
  return(adjmat)
}
MOFFL_adjmat <- create_MOFFL_adjmat(5)
MOFFL_graph <- graph_from_adjacency_matrix(MOFFL_adjmat)

#Bifan
bifan_adjmat <- matrix(c(0,0,1,1,
                         0,0,1,1,
                         0,0,0,0,
                         0,0,0,0), nrow = 4, byrow = T)
bifan_graph <- graph_from_adjacency_matrix(bifan_adjmat)

#DOR (Extended bifan)
create_DOR_adjmat <- function(n,m){
  #creates DOR adjmat with n inputs and m outputs
  adjmat <- matrix(0, n+m,n+m)
  adjmat[1:n,(n+1):(n+m)] <- 1
  return(adjmat)
}
DOR_adjmat <- create_DOR_adjmat(3, 3)
DOR_graph <- graph_from_adjacency_matrix(DOR_adjmat)


### MOTIFS INCLUDING INHIBITION (CONVERT TO TRINARY ADJMAT) ###
C1FFL_adjmat <- matrix(c(0,1,1,
                         0,0,1,
                         0,0,0), nrow = 3, byrow = T)
C1FFL_graph <- graph_from_adjacency_matrix(C1FFL_adjmat)

C2FFL_adjmat <- matrix(c(0,-1,-1,
                         0,0,1,
                         0,0,0), nrow = 3, byrow = T)
C2FFL_graph <- graph_from_adjacency_matrix(C2FFL_adjmat)

C3FFL_adjmat <- matrix(c(0,1,-1,
                         0,0,-1,
                         0,0,0), nrow = 3, byrow = T)
C3FFL_graph <- graph_from_adjacency_matrix(C3FFL_adjmat)

C4FFL_adjmat <- matrix(c(0,-1,1,
                         0,0,-1,
                         0,0,0), nrow = 3, byrow = T)
C4FFL_graph <- graph_from_adjacency_matrix(C4FFL_adjmat)

I1FFL_adjmat <- matrix(c(0,1,1,
                         0,0,-1,
                         0,0,0), nrow = 3, byrow = T)
I1FFL_graph <- graph_from_adjacency_matrix(I1FFL_adjmat)

I2FFL_adjmat <- matrix(c(0,-1,-1,
                         0,0,-1,
                         0,0,0), nrow = 3, byrow = T)
I2FFL_graph <- graph_from_adjacency_matrix(I2FFL_adjmat)
I3FFL_adjmat <- matrix(c(0,1,-1,
                         0,0,1,
                         0,0,0), nrow = 3, byrow = T)
I3FFL_graph <- graph_from_adjacency_matrix(I3FFL_adjmat)
I4FFL_adjmat <- matrix(c(0,-1,1,
                         0,0,1,
                         0,0,0), nrow = 3, byrow = T)
I4FFL_graph <- graph_from_adjacency_matrix(I4FFL_adjmat)