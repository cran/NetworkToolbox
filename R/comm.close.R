#' Community Closeness Centrality
#' @description Computes the community closeness centrality measure of each
#' community in a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param comm A vector or matrix corresponding to the
#' community each node belongs to
#' 
#' @param weighted Is the network weighted?
#' Defaults to TRUE.
#' Set to FALSE for weighted measures
#' 
#' @return A vector of community closeness centrality values for each specified
#' community in the network
#' (larger values suggest more central positioning)
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' comm <- igraph::walktrap.community(convert2igraph(abs(A)))$membership
#' 
#' #Weighted
#' result <- comm.close(A, comm)
#' 
#' #Unweighted
#' result <- comm.close(A, comm, weighted = FALSE)
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Community Closeness Centrality----
comm.close <- function (A, comm, weighted = TRUE)
{
    if(is.null(comm))
    {stop("comm must be input")}
    
    comm <- as.vector(comm)
    
    if(ncol(A)!=length(comm))
    {stop("length of comm does not match nodes in matrix")}
    
    uniq <- unique(comm)
    len <- length(uniq)
    
    allP <- pathlengths(A, weighted = weighted)$ASPLi
    mean.allP <- mean(allP)
    remove <- matrix(0,nrow=len,ncol=1)
    
    for(j in 1:len)
    {
        rem <- which(comm==uniq[j])
        
        remove[j,] <- (mean(allP[-rem]))-mean.allP
    }
    
    norm <- remove
    
    for(k in 1:len)
    {norm[k,] <- (remove[k,] - min(remove))/(max(remove)-min(remove))}
    
    norm <- as.vector(round(norm,3))
    
    names(norm) <- uniq
    
    return(norm)
}
#----