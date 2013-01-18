`recalculateLogLik` <- 
   function(model, fixef = model@fixef, vcor = VarCorr(model)) {
        if (class(model) != 'mer') stop("model should be an object of the class mer")
        sig2 <- attr(vcor, "sc")
        D = matrix(0, max(model@Gp), max(model@Gp))
        npos <- diff(model@Gp)
        for (i in seq_along(vcor)) {
          Dpart <- kronecker(vcor[[i]], diag(npos[i]/ncol(vcor[[i]])))
          inds <- (model@Gp[i]+1):model@Gp[i+1]
          D[inds, inds] <- Dpart
        }
        V1 <- Matrix:::t(model@Zt) %*% Matrix(D) %*% model@Zt
        V2 <- sig2 * Diagonal(ncol(model@Zt))
        V <- V1 + V2
        res <- model@y - model@X %*% model@fixef
        z <- t(res) %*% solve(V,res)
        as.numeric(- Matrix:::determinant(V,TRUE)$modulus/2 - ncol(model@Zt)*log(2*pi)/2 - z/2)
}

`groupDisp` <-
  function(formula, data, var) {
    model0 <- lmer(formula, data, REML=FALSE)
    lev <- levels(data[,var])
    logs1 <- sapply(lev, function(l) {
      model1 <- lmer(formula, data[data[,var] != l, ], REML=FALSE)
      recalculateLogLik(model0, model1@fixef, VarCorr(model1))
    })
    logLik(model0) - logs1
}
  
`obsDisp` <-
  function(formula, data, inds=1:nrow(data)) {
    model0 <- lmer(formula, data, REML=FALSE)
    logs1 <- sapply(inds, function(l) {
      model1 <- lmer(formula, data[-l, ], REML=FALSE)
      recalculateLogLik(model0, model1@fixef, VarCorr(model1))
    })
    names(logs1) = inds
    logLik(model0) - logs1
  }

