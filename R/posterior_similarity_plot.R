#' Posterior Similarity plot
#'

plotpsm <- function(psm.ind,
                    psm.tot,
                    method="complete",
                    ...,
                    plot.type = c('ind','tot','both')){

  if(any(psm.tot !=t(psm.tot)) | any(psm.tot >1) | any(psm.tot < 0) | sum(diag(psm.tot)) != nrow(psm.tot) ){
    stop("psm.tot must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}

  

  ##------------------------------ without distinguishing between datasets ------------------------------
  hc=hclust(as.dist(1-psm.tot), method = method, members = NULL)
  psm_hc=psm.tot
  n=nrow(psm.tot)
  psm_hc[1:n,]=psm_hc[rev(hc$order),]
  psm_hc[,1:n]=psm_hc[,hc$order]

  if(plot.type %in% c('tot', 'both')){
    
    
    par(mfrow=c(1,1))
    
    image.plot(1:n,
               1:n,
               psm_hc,
               col=rev(heat.colors(100)),
               axis.args=list(yaxt='n', xaxt='n'),
               ...)
    
    
    
  }


  # Dimension of datasets
  M <- length(psm.ind)
  n <- rep(0,M)
  for(m in 1:M){
    n[m] <- nrow(psm.ind[[m]][[m]])
  }
  n <- c(0,cumsum(n))

  ##----------------------------------- Distinguishing between datasets ----------------------------------
  psm_matrix_output <- matrix(0, nrow=max(n), ncol=max(n))
  hc_diag <- NULL

  ##-- For diagonal
  for(m in 1:M){
    hc_diag[[m]] <- hclust(as.dist(1-psm.ind[[m]][[m]]), method=method, members=NULL)
    psm_hc_diag <- psm.ind[[m]][[m]]
    psm_hc_diag[1:nrow(psm_hc_diag),] <- psm_hc_diag[rev(hc_diag[[m]]$order),]
    psm_hc_diag[,1:nrow(psm_hc_diag)] <- psm_hc_diag[,hc_diag[[m]]$order]
    psm_matrix_output[(n[m]+1):n[m+1], (max(n)-n[m+1]+1):(max(n)-n[m])] <- psm_hc_diag
  }

  ##-- For off-diagonal
  for(m1 in 2:M){
    for(m2 in 1:(m1-1)){
      psm_m1_m2 <- t(psm.ind[[m1]][[m2]])
      psm_m1_m2[1:nrow(psm_m1_m2),] <- psm_m1_m2[rev(hc_diag[[m2]]$order),]
      psm_m1_m2[,1:ncol(psm_m1_m2)] <- psm_m1_m2[,hc_diag[[m1]]$order]
      psm_matrix_output[(n[m2]+1):n[m2+1], (max(n)-n[m1+1]+1):(max(n)-n[m1])] <- psm_m1_m2
    }
  }

  for(m1 in 1:(M-1)){
    for(m2 in (m1+1):M){
      psm_m1_m2 <- psm.ind[[m2]][[m1]]
      psm_m1_m2[1:nrow(psm_m1_m2),] <- psm_m1_m2[rev(hc_diag[[m2]]$order),]
      psm_m1_m2[,1:ncol(psm_m1_m2)] <- psm_m1_m2[,hc_diag[[m1]]$order]
      psm_matrix_output[(n[m2]+1):n[m2+1], (max(n)-n[m1+1]+1):(max(n)-n[m1])] <- psm_m1_m2
    }
  }

  n.max <- max(n)

  if(plot.type %in% c('ind', 'both')){
    
    image.plot(1:n.max, 
               1:n.max,
               psm_matrix_output, 
               col=rev(heat.colors(100)),
               axis.args=list(yaxt='n', xaxt='n'),
               ...)
    
    
    
    abline(v=n[-c(1,M+1)], lwd=3)
    abline(h=n.max - n[-c(1,M+1)], lwd=3)
    
    
  }

}

