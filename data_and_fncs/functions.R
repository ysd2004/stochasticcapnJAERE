aproxdef <- function(deg,lb,ub,delta){
  if (delta > 1 | delta < 0)
    stop("delta should be in [0,1]!")
  
  dn <- length(deg)
  dl <- length(lb)
  du <- length(ub)
  
  if (dn != dl)
    stop("Dimension Mismatch: Stock dimension != lower bounds dimension")
  else if (dn != du)
    stop("Dimension Mismatch: Stock dimension != upper bounds dimension")
  else
    param <- list(degree = deg, lowerB = lb, upperB = ub, delta = delta)
  return(param)
}


chebbasisgen <- function(stock,npol,a,b,dorder=NULL){
  nknots <- length(stock)
  z <- (2*stock-b-a)/(b-a)
  if (npol < 4)
    stop("Degree of Chebyshev polynomial must be greater than 3!")
  
  bvec <- cbind(matrix(rep(1,nknots),ncol=1), matrix(z, ncol=1))
  colnames(bvec) <- c("j=1","j=2")
  
  for (j in 2:(npol-1)){
    Tj <- matrix(2*z*bvec[,j] - bvec[,j-1],ncol=1)
    colnames(Tj) <- paste0("j=",j+1)
    bvec <- cbind(bvec,Tj)
  }
  
  if (is.null(dorder)){
    res <- bvec
  }
  else if (dorder==1){
    bvecp <- cbind(matrix(rep(0,nknots),ncol=1), matrix(rep(2/(b-a),nknots), ncol=1))
    colnames(bvecp) <- c("j=1","j=2")
    for (j in 2:(npol-1)){
      Tjp <- matrix((4/(b-a))*bvec[,j]+2*z*bvecp[,j] - bvecp[,j-1],ncol=1)
      colnames(Tjp) <- paste0("j=",j+1)
      bvecp <- cbind(bvecp,Tjp)
    }
    res <- bvecp
  }
  else if (dorder>=2){
    bvecp <- cbind(matrix(rep(0,nknots),ncol=1), matrix(rep(2/(b-a),nknots), ncol=1))
    colnames(bvecp) <- c("j=1","j=2")
    for (j in 2:(npol-1)){
      Tjp <- matrix((4/(b-a))*bvec[,j]+2*z*bvecp[,j] - bvecp[,j-1],ncol=1)
      colnames(Tjp) <- paste0("j=",j+1)
      bvecp <- cbind(bvecp,Tjp)
    }
    bvec <- bvecp
    for (d in 2:dorder){
      bvecp <- cbind(matrix(rep(0,2*nknots),ncol=2))
      colnames(bvecp) <- c("j=1","j=2")
      for (j in 2:(npol-1)){
        Tjp <- matrix(((4*d)/(b-a))*bvec[,j]+2*z*bvecp[,j] - bvecp[,j-1],ncol=1)
        colnames(Tjp) <- paste0("j=",j+1)
        bvecp <- cbind(bvecp,Tjp)
      }
      bvec <- bvecp
    }
    res <- bvecp
  }
  else
    stop("dorder must be a positive integer or NULL!")
  
  return(res)
}

chebnodegen <- function(n,a,b){
  d1 <- length(a)
  d2 <- length(b)
  
  si <- (a+b)*0.5 + (b-a)*0.5*cos(pi*((seq((n-1),0,by=-1)+0.5)/n))
  
  if (d1 != d2){
    print("Dimension Mismatch: dim(upper) != dim(lower)")
  }
  return(si)
}

vaprox <- function (aproxspace, stock, sdot, w, covmat=NULL){
  deg <- aproxspace[["degree"]]
  lb <- aproxspace[["lowerB"]]
  ub <- aproxspace[["upperB"]]
  delta <- aproxspace[["delta"]]
  dd <- length(deg)
  m.type <- "deterministic"
  if (!is.null(covmat)){
    m.type <- "stochastic"
  }
  
  sdata <- as.matrix(cbind(stock,sdot,w),ncol=(dd*2)+1)
  if (dim(sdata)[2] != (2 * dd + 1)) 
    stop("The number of stocks are not the same with the number of sdots!")
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
  sdot <- lapply((dd + 1):(2 * dd), function(k) sdata[, k])
  w <- sdata[, (2 * dd + 1)]
  
  d0nodes <- lapply(1:dd, function(k) chebbasisgen(st[[k]], deg[k], lb[k], ub[k]))
  d1nodes <- lapply(1:dd, function(k) chebbasisgen(st[[k]], deg[k], lb[k], ub[k],dorder = 1))
  
  fphi <- matrix(1)
  sphi <- matrix(rep(0, prod(sapply(st, length)) * prod(deg)), 
                 ncol = prod(deg))
  
  for (di in dd:1) {
    ftemp <- d0nodes[[di]]
    fphi <- kronecker(fphi, ftemp)
    stempi <- d1nodes[[di]]
    sphitemp <- matrix(1)
    for (dj in dd:1) {
      if (dj != di) {
        stemp <- d0nodes[[dj]]
      }
      else stemp <- stempi
      sphitemp <- kronecker(sphitemp, stemp)
    }
    sphi <- sphi + sphitemp * sdot[[di]]
  }
  nsqr <- delta * fphi - sphi
  
  if (m.type == "stochastic"){
    
    d2nodes <- lapply(1:dd, function(k) chebbasisgen(st[[k]], deg[k], lb[k], ub[k],dorder = 2))
    
    tphi <- matrix(rep(0, prod(sapply(st, length)) * prod(deg)), 
                   ncol = prod(deg))
    cn <- (dd^2)+1
    
    if ((cn-1) != dim(covmat)[2]) 
      stop("The number of columns in covmat is not correct!")
    
    for (vi in dd:1){
      for (vj in dd:1){
        ttemp <- matrix(1)
        for (vk in dd:1){
          ttempi <- d0nodes[[vk]]
          if (vk %in% c(vi,vj)){
            ttempi <- d1nodes[[vk]]
          }
          if ((vk == vi) & (vk == vj)){
            ttempi <- d2nodes[[vk]]
          }
          ttemp <- kronecker(ttemp,ttempi)
        }
        cn <- cn-1
        tphi <- tphi + diag(covmat[,cn])%*%ttemp
      }
    }
    nsqr <- delta * fphi - sphi -(1/2)*tphi
  }
  
  if (dim(fphi)[1] == dim(fphi)[2]) {
    coeff <- solve(nsqr, w)
    res <- list(degree = deg, lowerB = lb, upperB = ub, delta = delta, 
                coefficient = coeff, model.type = m.type)
  }
  else if (dim(fphi)[1] != dim(fphi)[2]) {
    coeff <- solve(t(nsqr) %*% nsqr, t(nsqr) %*% w)
    res <- list(degree = deg, lowerB = lb, upperB = ub, delta = delta, 
                coefficient = coeff, model.type = m.type)
  }
  return(res)
}

vsim <- function (vcoeff, stock, wval = NULL) 
{
  deg <- vcoeff[["degree"]]
  lb <- vcoeff[["lowerB"]]
  ub <- vcoeff[["upperB"]]
  delta <- vcoeff[["delta"]]
  coeff <- vcoeff[["coefficient"]]
  m.type <- vcoeff[["model.type"]]
  dd <- length(deg)
  
  if (is.data.frame(stock)) {
    st <- as.matrix(stock)
  }
  else if (is.matrix(stock)) {
    st <- stock
  }
  else if (is.numeric(stock)) {
    st <- as.matrix(stock,ncol=dd)
  }
  else stop("st is not a matrix, numeric or data.frame!")
  
  nnodes <- dim(st)[1]
  
  accp <- matrix(rep(0, nnodes * dd), ncol = dd)
  Bmat <- matrix(rep(0, nnodes * prod(deg)), ncol = prod(deg))
  for (di in 1:dd) {
    Bprime <- matrix(rep(0, nnodes * prod(deg)), ncol = prod(deg))
    for (ni in 1:nnodes) {
      sti <- st[ni, ]
      dk <- dd - di + 1
      fphi <- matrix(1)
      ftemp <- chebbasisgen(sti[di], deg[di], lb[di], ub[di])
      sphi <- matrix(1)
      stempd <- chebbasisgen(sti[di], deg[di], lb[di], 
                             ub[di], dorder = 1)
      for (dj in 1:dd) {
        dk2 <- dd - dj + 1
        if (dk2 != di) {
          ftemp <- chebbasisgen(sti[dk2], deg[dk2], lb[dk2], 
                                ub[dk2])
          stemp <- ftemp
        }
        else stemp <- stempd
        fphi <- kronecker(fphi, ftemp)
        sphi <- kronecker(sphi, stemp)
      }
      Bmat[ni, ] <- fphi
      Bprime[ni, ] <- sphi
    }
    accp[, di] <- Bprime %*% coeff
  }
  iwhat <- accp * st
  iw <- matrix(rowSums(iwhat), ncol = 1)
  vhat <- Bmat %*% coeff
  colnames(accp) <- paste("acc.price", 1:dd, sep = "")
  colnames(iwhat) <- paste("iw", 1:dd, sep = "")
  colnames(iw) <- c("iw")
  if (is.null(wval) == 1) {
    wval <- "wval is not provided"
  }
  res <- list(shadowp = accp, iweach = iwhat, iw = iw, vfun = vhat, 
              stock = st, wval = wval, model.type = m.type)
  return(res)
}

unigrids <- function (nnodes, lb, ub, rtype = NULL) 
{
  dn <- length(nnodes)
  dl <- length(lb)
  du <- length(ub)
  if (dn != dl) 
    stop("Dimension Mismatch: Stock dimension != lower bounds dimension")
  else if (dn != du) 
    stop("Dimension Mismatch: Stock dimension != upper bounds dimension")
  else uniknots <- mapply(seq, from = lb, to = ub, len = nnodes, 
                          SIMPLIFY = FALSE)
  if (is.null(rtype)) {
    return(uniknots)
  }
  else if (rtype == "list") {
    return(uniknots)
  }
  else if (rtype == "grid") {
    gridcomb <- as.matrix(do.call(expand.grid, uniknots), 
                          ncol = dn)
    return(gridcomb)
  }
  else stop("type should be: NULL, list, or grid!")
}


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


vaproxsc <- function(aproxspace, stock, growthfun, covmat=NULL, itermax=300, tol=10^(-8)){
  ## profit function
  wfun <- function(s,q){
    if (eta == 1){
      wout <- b*log(q)-(c/(s^gamma))*q
    } else {
      wout <- (b^(1/eta))/(1-(1/eta))*(q^(1-(1/eta)))-(c/(s^gamma))*q
    }
    return(wout)
  }
  ## drift term
  mufun <- function(s,q){
    muout <- growthfun(s) - q
    return(muout)
  }
  ## q-optimum
  qfun <- function(s,Vs){
    qout <- b*((Vs+c/(s^gamma))^(-eta))
    return(qout)
  }
  ## initialization
  iter <- 0
  error <- 1
  q <- qfun(stock,rep(0,length(stock)))
  w <- wfun(stock,q)
  mus <- mufun(stock,q)
  cv <- vaprox(Aspace,stock,mus,w,covmat)
  cv$coefficient <- rep(0,length(cv$coefficient))
  vout <- vsim(cv,stock)
  
  while(error > tol && iter < itermax){
    vin <- vout
    Vs <- vin$shadowp
    v0 <- vin$vfun
    q <- qfun(stock,Vs)
    w <- wfun(stock,q)
    mus <- r*stock*(1-(stock/K)) - q
    cv <- vaprox(Aspace,stock,mus,w,covmat)  ## do not put sigs value here!!
    vout <- vsim(cv,stock)
    v <- vout$vfun 
    error <- max(abs(v-v0))
    iter <- iter+1
    message('iteration:',iter,'   error:',error)
  }
  return(cv)
}

vaproxpss <- function (aproxspace, stock, sdot, w, hs){
  deg <- aproxspace[["degree"]]
  lb <- aproxspace[["lowerB"]]
  ub <- aproxspace[["upperB"]]
  delta <- aproxspace[["delta"]]
  dd <- length(deg)
  m.type <- "Poisson Jump"
  
  sdata <- as.matrix(cbind(stock,sdot,w),ncol=(dd*2)+1)
  if (dim(sdata)[2] != (2 * dd + 1)) 
    stop("The number of stocks are not the same with the number of sdots!")
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
  sdot <- lapply((dd + 1):(2 * dd), function(k) sdata[, k])
  w <- sdata[, (2 * dd + 1)]
  
  d0nodes <- lapply(1:dd, function(k) chebbasisgen(st[[k]], deg[k], lb[k], ub[k]))
  d1nodes <- lapply(1:dd, function(k) chebbasisgen(st[[k]], deg[k], lb[k], ub[k],dorder = 1))
  
  fphi <- matrix(1)
  sphi <- matrix(rep(0, prod(sapply(st, length)) * prod(deg)), 
                 ncol = prod(deg))
  
  for (di in dd:1) {
    ftemp <- d0nodes[[di]]
    fphi <- kronecker(fphi, ftemp)
    stempi <- d1nodes[[di]]
    sphitemp <- matrix(1)
    for (dj in dd:1) {
      if (dj != di) {
        stemp <- d0nodes[[dj]]
      }
      else stemp <- stempi
      sphitemp <- kronecker(sphitemp, stemp)
    }
    sphi <- sphi + sphitemp * sdot[[di]]
  }
  nsqr <- matrix(diag(as.vector(delta+hs)),ncol=deg)%*%fphi - sphi
  
  if (dim(fphi)[1] == dim(fphi)[2]) {
    coeff <- solve(nsqr, w)
    res <- list(degree = deg, lowerB = lb, upperB = ub, delta = delta, 
                coefficient = coeff, model.type = m.type)
  }
  else if (dim(fphi)[1] != dim(fphi)[2]) {
    coeff <- solve(t(nsqr) %*% nsqr, t(nsqr) %*% w)
    res <- list(degree = deg, lowerB = lb, upperB = ub, delta = delta, 
                coefficient = coeff, model.type = m.type)
  }
  return(res)
}

