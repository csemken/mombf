######################################################################################
## DIRICHLET
######################################################################################


#Density of a Dirichlet(q) ditribution evaluated at x
ddir= function (x, q, logscale=TRUE) {
    if (is.vector(x)) x= matrix(x,nrow=1)
    if (length(q)==1) q= rep(q,ncol(x))
    sel= rowSums((x<0) | (x>1))>0
    ans= double(nrow(x))
    ans[sel]= -Inf
    xnorm= t(x[!sel,,drop=FALSE])/rowSums(x[!sel,,drop=FALSE])
    ans[!sel]= colSums((q - 1) * log(xnorm)) + lgamma(sum(q)) - sum(lgamma(q))
    if (!logscale) ans= exp(ans)
    return(ans)
}


######################################################################################
## INVERSE GAMMA
######################################################################################

#Adapted from LaplacesDemon
dinvgamma= function (x, shape = 1, scale = 1, log = FALSE) {
    x <- as.vector(x)
    shape <- as.vector(shape)
    scale <- as.vector(scale)
    if (any(shape <= 0) | any(scale <= 0)) stop("The shape and scale parameters must be positive.")
    NN <- max(length(x), length(shape), length(scale))
    x <- rep(x, len = NN)
    shape <- rep(shape, len = NN)
    scale <- rep(scale, len = NN)
    alpha <- shape
    beta <- scale
    ans <- alpha * log(beta) - lgamma(alpha) - {alpha + 1} * log(x) - {beta/x}
    if (log == FALSE) ans <- exp(ans)
    return(ans)
}


#Inverse gamma quantiles (not implemented currently, returning NAs)
qinvgamma= function(p, shape = 1, scale = 1, log = FALSE) {
  return(rep(NA,length(p)))
}



#Inverse gamma CDF
pinvgamma= function(x, shape = 1, scale = 1, log = FALSE) {
    shape <- as.vector(shape)
    scale <- as.vector(scale)
    if (any(shape <= 0) | any(scale <= 0)) stop("The shape and scale parameters must be positive.")
    NN <- max(length(x), length(shape), length(scale))
    x <- rep(x, len = NN)
    shape <- rep(shape, len = NN)
    scale <- rep(scale, len = NN)
    alpha <- shape
    beta <- scale
    ans= gammainc(alpha, beta/x)[3]
    if (log) ans= log(ans)
    return(ans)
}


#Incomplete gamma function and regularized incomplete gamma function, adapted from R package pracma
gammainc <-  function (x, a) {
    if (!is.numeric(a) || !is.numeric(x)) stop("All arguments must be real numbers.")
    if (length(a) > 1 || length(x) > 1)  stop("Arguments must be of length 1; function is not vectorized.")
    if (x == 0 && a == 0) return(1)
    if (x == 0) return(gamma(a))
    if (a < 0) stop("Argument 'a' must be real and nonnegative.")
    if (x > 0) { xam <- -x + a * log(x) } else { xam <- -x + a * log(x + (0+0i)) }
    if (abs(xam) > 700 || abs(a) > 170) {
        warning("Arguments 'x' and/or 'a' are too large.")
        return(NA)
    }
    gin <- gim <- gip <- 0
    if (x == 0) {
        ga <- gamma(a)
        gim <- ga
        gip <- 0
    } else if (x <= 1 + a) {
        s <- 1/a
        r <- s
        for (k in 1:60) {
            r <- r * x/(a + k)
            s <- s + r
            if (abs(r/s) < 1e-15) break
        }
        gin <- exp(xam) * s
        ga <- gamma(a)
        gip <- gin/ga
        gim <- ga - gin
    } else if (x > 1 + a) {
        t0 <- 0
        for (k in 60:1) { t0 <- (k - a)/(1 + k/(x + t0)) }
        gim <- exp(xam)/(x + t0)
        ga <- gamma(a)
        gin <- ga - gim
        gip <- 1 - gim/ga
    }
    return(c(lowinc = Re(gin), uppinc = Re(gim), reginc = Re(gip)))
}



######################################################################################
## INVERSE WISHART
######################################################################################

#log of the Multivariate gamma function
lmgamma= function(p,a) { 0.25*p*(p-1)*log(pi) + sum(lgamma(a + 0.5*(1-(1:p)))) }

diwish <- function(Sigma, nu, S, logscale=FALSE) {
    #Inverse Wishart density, adapted from LaplacesDemon package
    if (!is.matrix(Sigma)) Sigma <- matrix(Sigma)
    if (!is.matrix(S)) S <- matrix(S)
    if (!identical(dim(S), dim(Sigma))) stop("The dimensions of Sigma and S differ.")
    if (nu < nrow(S)) stop("The nu parameter is less than the dimension of S.")
    p <- nrow(Sigma)
    detSigma= as.numeric(determinant(Sigma,logarithm=TRUE)$modulus)
    detS= as.numeric(determinant(S,logarithm=TRUE)$modulus)
    ans <- -(nu * p/2) * log(2) - lmgamma(p,.5*nu) + (nu/2) * detS - ((nu + p + 1)/2) * detSigma - 0.5 * sum(diag((S %*% solve(Sigma))))
    if (!logscale) ans <- exp(ans)
    return(ans)
}


riwish <- function(nu, S, Sinv) {
    if (missing(Sinv)) Sinv <- solve(S)
    solve(rWishart(1, df=nu, Sigma=Sinv)[,,1])
}
