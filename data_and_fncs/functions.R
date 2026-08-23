chebc <- function (aproxspace, stock, w){
  deg <- aproxspace[["degree"]]
  lb <- aproxspace[["lowerB"]]
  ub <- aproxspace[["upperB"]]
  delta <- aproxspace[["delta"]]
  dd <- length(deg)
  
  sdata <- as.matrix(cbind(stock,w),ncol=(dd+1))
  if (dd > 1) {
    ordername <- paste0("sdata[,", dd, "]")
    for (di in 2:dd) {
      odtemp <- paste0("sdata[,", dd - di + 1, "]")
      ordername <- paste(ordername, odtemp, sep = ",")
    }
    ordername <- paste("sdata[order(", ordername, "),]", 
                       sep = "")
    sdata <- eval(parse(text = ordername))
  }
  else sdata <- sdata[order(sdata[, 1]), ]
  
  st <- lapply(1:dd, function(k) unique(sdata[, k]))
  w <- sdata[, (dd + 1)]
  
  d0nodes <- lapply(1:dd, function(k) chebbasisgen(st[[k]], deg[k], lb[k], ub[k]))
  fphi <- matrix(1)
  
  for (di in dd:1) {
    ftemp <- d0nodes[[di]]
    fphi <- kronecker(fphi, ftemp)
  }
  nsqr <- fphi
  
  
  if (dim(fphi)[1] == dim(fphi)[2]) {
    coeff <- solve(nsqr, w)
    res <- list(degree = deg, lowerB = lb, upperB = ub, delta = delta, 
                coefficient = coeff)
  }
  else if (dim(fphi)[1] != dim(fphi)[2]) {
    coeff <- solve(t(nsqr) %*% nsqr, t(nsqr) %*% w)
    res <- list(degree = deg, lowerB = lb, upperB = ub, delta = delta, 
                coefficient = coeff)
  }
  return(res)
}