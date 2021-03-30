
## Create the scaled precision of the structured component of the BYM2 model
#
# adj.mat: matrix defining neighbors
#          NOTE: that the diag(adj.mat) is all 0s
# sparse.mat: LOGICAL. if T, return sparse matrix object
#
# returns matrix of the scaled precision of the structured portion of the BYM2 model
make_BYM2_struct_scaled_prec <- function(adj.mat,
                                         sparse.mat = TRUE){

  if(sparse.mat){
    # first, make the unscaled ICAR precision
    # NOTE: this is singular!
    Q.u <- Diagonal(n = nrow(adj.mat), x = rowSums(adj.mat)) - adj.mat

    # add a small jitter to the diagonal for numerical stability (and invertibility)
    # not necessary, but recommended
    Q.u <- Q.u + Diagonal(nrow(adj.mat)) * max(diag(Q.u)) * sqrt(.Machine$double.eps)
  }else{
    # first, make the unscaled ICAR precision
    # NOTE: this is singular!
    Q.u <- diag(nrow = nrow(adj.mat),
                ncol = nrow(adj.mat), x = rowSums(adj.mat)) - adj.mat

    # add a small jitter to the diagonal for numerical stability (and invertibility)
    # not necessary, but recommended
    Q.u <- Q.u + diag(nrow = nrow(adj.mat),
                      ncol = nrow(adj.mat),
                      x = 1) * max(diag(Q.u)) * sqrt(.Machine$double.eps)
  }

  ## # Compute the geometric mean of the variances, which are on the
  ## # diagonal of Q.inv
  ## scaling_factor <- exp(mean(log(diag(solve(Q.u)))))

  ## # and scale
  ## Q.s2 <- Q.u / scaling_factor

  Q.s <- inla.scale.model(Q.u,
                          constr = list(A = matrix(1, nrow = 1, ncol = ncol(Q.u)), e = 0))

  # for some reason, the Matrix pkg doesn't always think this is
  # symmetric (even though all ij entries equal ji entries as far as I
  # can tell), so we force symmetry

  if(!isSymmetric(Q.s)){
    Q.s <- forceSymmetric(Q.s)
  }

  return(Q.s)
}

## Calculate and return the precision matrix for a (scaled) BYM2 model
#
# adj.mat: adjacency matrix
# phi: proportion of total var that is structured
# tau: total precision of struct + unstruct
# sparse.mat: LOGICAL. if T, return sparse matrix object
#
# returns the (sparse) precision matrix for w = (w1, w2)
# w1 = total effect
# w2 = scaled structured effect
make_BYM2_joint_unstruct_scaled_struc_prec_mat <- function(adj.mat,
                                                           phi,
                                                           tau,
                                                           sparse.mat = TRUE
                                                           ){

  # get the number of regions
  n <- nrow(adj.mat)

  # first, make the scaled precision matrix of the structured component
  Q.s <- make_BYM2_struct_scaled_prec(adj.mat = adj.mat, sparse.mat = sparse.mat)

  ## we make the joint prec (unstruct + struct) and then populate the
  ## 4 (2x2) blocks of the matrix

  if (sparse.mat) {
    # joint mat
    prec.mat <- Matrix(0, nrow = 2 * n, ncol = 2 * n)

    # blocks
    block11 <- diag(ncol = n, nrow = n, x = tau / (1 - phi))
    block12 <- diag(ncol = n, nrow = n, x = -sqrt(phi * tau) / (1 - phi))
    # block21 == block12
    block22 <- Q.s + diag(ncol = n, nrow = n, x = phi / (1 - phi))
  }else{

    # joint mat
    prec.mat <- matrix(0, nrow = 2 * n, ncol = 2 * n)

    # blocks
    block11 <- Diagonal(n = n, x = tau / (1 - phi))
    block12 <- Diagonal(n = n, x = -sqrt(phi * tau) / (1 - phi))
    # block21 == block12
    block22 <- Q.s + Diagonal(n = n, x = phi / (1 - phi))
  }

  # populate
  prec.mat[1:n, 1:n] <- block11
  prec.mat[1:n, (n + 1):(2 * n)] <- block12
  prec.mat[(n + 1):(2 * n), 1:n] <- block12
  prec.mat[(n + 1):(2 * n), (n + 1):(2 * n)] <- block22

  return(prec.mat)
}

## Simulate a multivar norm given a precision mat
#
# mu: mean vec
# prec: precision matrix
# n.sims: numer of draws
# sumtozero: logical, if T, modify draws to sum-to-zero

# returns a matrix of draws (columns) across the dimension of the RV (rows)
rmvnorm_prec <- function(mu, prec, n.sims, sumtozero = FALSE) {

  if(length(mu) == 1){
    mu <- rep(mu, nrow(prec))
  }

  z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L <- Matrix::Cholesky(prec, super=TRUE)
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  x <- mu + z

  if(!sumtozero){

    return(x)

  }else{

    # we need to adjust each draw to constrain it
    # (efficiently) make some relevant objects
    A <- rep(1, nrow(prec))
    Qinv.A <- solve(prec, A)

    # this is how we do it for 1 draw
    app.constr <- function(x){
      x - Qinv.A %*% ((A %*%Qinv.A) ^ (-1)) %*% (A %*% x - 0)
    }

    if(n.sims == 1){
      x.c <- app.constr(x)
    }else{
      x.c <- do.call("cbind", apply(x, 2, app.constr))
    }

    return(x.c)

  } # else(!sumtozero)

}

## Simulate a BYM2 RV, sampling simultaneously the total (unstruct +
## struct) and the structured components
#
# mu: mean vec
# prec: joint precision matrix of bym2 (first block is unstructured,
#    second block is structured)
# n.sims: numer of draws
# sumtozero: logical, if T, modify draws to sum-to-zero ONLY the
#    structured component
# constrain.idx vec. if sumtozero, which indices should be forced to sum-to-zero?

# returns a matrix of draws (columns) across the dimension of the RV (rows)
rbym2_simul_1constr_prec <- function(mu, prec, n.sims,
                                     sumtozero = FALSE,
                                     constrain.idx = NULL) {

  if(length(mu) == 1){
    mu <- rep(mu, nrow(prec))
  }

  x <- rmvnorm_prec(mu = mu, prec = prec, n.sims = n.sims, sumtozero = sumtozero)

  if(!sumtozero){

    return(x)

  }else{

    # we need to adjust each draw to constrain it s.t.
    # Ax.c = e
    # (efficiently) make some relevant objects

    # by default, constrain the second half of the vector - which I
    # use to be the structured components of the bym2
    if(is.null(constrain.idx)) constrain.idx <- tail(1:length(mu), length(mu) / 2)

    # only applying the sum - to - zero constraint on the second half
    # of the effects (the structured portion)
    A <- rep(0, length(mu))
    A[constrain.idx] <- 1
    A <- matrix(A,
                nrow = 1, byrow = T)
    e <- matrix(rep(0, 1), ncol = 1)

    Qinv.A <- Matrix::solve(prec, t(A))

    # this is how we do it for 1 draw
    app.constr <- function(x){
      x - Qinv.A %*% ((A %*% Qinv.A) ^ (-1)) %*% (A %*% x - e)
    }

    if(n.sims == 1){
      x.c <- app.constr(x)
    }else{
      x.c <- do.call("cbind", apply(x, 2, app.constr))
    }

    return(list(x = x, x.c = x.c))

  } # else(!sumtozero)

}

## function to simulate data given a population and latent field
## values
#
# pop: population to observe
# int: intercept (added to field before sampling)
# field: GMRfield values (latent)
# lik: likelihhod for the data. one of:
#  "poisson"
#  "normal" # TODO
#  "binom" # TODO
#
# returns a data.table with two columns: obs (data) and risk (field)

sim_field_data <- function(pop, int, field, lik = "poisson"){

  if(lik != "poisson"){
    stop("only poisson has been implemented")
  }

  # convert the latent field using the canonical link for the selected
  # likelihood

  lik.link.dict <- list("poisson" = "exp",
                        "binom"   = "plogis",
                        "normal"  = "identity")

  invlink.func <- lik.link.dict[[lik]]

  risk.field <- do.call(invlink.func, list(int + field))

  # simulate data

  if(lik == "poisson"){
    obs <- rpois(n = length(risk.field), lambda = risk.field * pop)
  }

  if (lik == "binom") {
    TODO
  }

  if (lik == "normal") {
    TODO
  }

  return(data.table('obs'=obs, 'risk'=risk.field))

}












#######################################
## APPENDIX (TEST/UNUSED) CODE BELOW ##
#######################################

## ## Simulate a BYM2 RV, sampling the unstructured and structured
## ## components separately - and constraining them separately to sum to
## ## zero
## #
## # tau: total precision
## # phi: proportion of total that is structured
## # *.mu: mean vec for *
## # struct.prec: scaled precision matrix of bym2
## # n.sims: numer of draws
## # sumtozero: logical, if T, modify draws to sum-to-zero both the total
## # and the structured components
## #
## # returns a matrix of draws (columns) across the dimension of the RV (rows)
## rbym2_prec <- function(tau, phi,
##                        struct.mu = 0, unstruct.mu = 0,
##                        struct.prec, n.sims, sumtozero = FALSE) {

##   if(length(struct.mu) == 1){
##     struct.mu <- rep(struct.mu, nrow(struct.prec))
##   }
##   if(length(unstruct.mu) == 1){
##     unstruct.mu <- rep(unstruct.mu, nrow(struct.prec))
##   }

##   unstruct <- rmvnorm_prec(mu = unstruct.mu, prec = Diagonal(n = nrow(struct.prec), x = 1),
##                            n.sims = n.sims, sumtozero = sumtozero)
##   struct   <- rmvnorm_prec(mu = struct.mu, prec = struct.prec,
##                            n.sims = n.sims, sumtozero = sumtozero)

##   # combine (appropriately) into total followed by vector of struct
##   total <- 1 / sqrt(tau) * (sqrt(1 - phi) * unstruct + sqrt(phi) * struct)
##   return(rbind(total, struct))

## }

## ## Simulate a BYM2 RV, sampling simultaneously the total (unstruct +
## ## struct) and the structured components
## #
## # mu: mean vec

## # prec: joint precision matrix of bym2 (first block is unstructured,
## #    second block is structured)
## # n.sims: numer of draws
## # sumtozero: logical, if T, modify draws to sum-to-zero both the total
## #    and the structured components

## # returns a matrix of draws (columns) across the dimension of the RV# (rows)
## rbym2_simul_prec <- function(mu, prec, n.sims, sumtozero = FALSE) {

##   if(length(mu) == 1){
##     mu <- rep(mu, nrow(prec))
##   }

##   x <- rmvnorm_prec(mu = mu, prec = prec, n.sims = n.sims, sumtozero = sumtozero)

##   if(!sumtozero){

##     return(x)

##   }else{

##     # we need to adjust each draw to constrain it s.t.
##     # Ax.c = e
##     # (efficiently) make some relevant objects

##     A <- matrix(c(rep(1:0, each = nrow(prec) / 2),
##                   rep(0:1, each = nrow(prec) / 2)),
##                 nrow = 2, byrow = T)
##     e <- matrix(rep(0, 2), ncol = 1)

##     Qinv.A <- Matrix::solve(prec, t(A))

##     # this is how we do it for 1 draw
##     app.constr <- function(x){
##       x - Qinv.A %*% ((A %*% Qinv.A) ^ (-1)) %*% (A %*% x - e)
##     }

##     if(n.sims == 1){
##       x.c <- app.constr(x)
##     }else{
##       x.c <- do.call("cbind", apply(x, 2, app.constr))
##     }

##     return(list(x = x, x.c = x.c))

##   } # else(!sumtozero)

## }


## dlogitbeta <- function(logit_p, a, b, give_log = 0){

##   part1 = log(gamma(a + b)) - log(gamma(a))  - log(gamma(b));
##   part2 = (a - 1) * (logit_p - log(1 + exp(logit_p)));
##   part3 = (b - 1) * log( 1 - exp(logit_p)/(1 + exp(logit_p)));
##   part4 =  logit_p - 2 * log( 1 + exp(logit_p));

##   logres = part1 + part2 + part3 + part4;

##   if(give_log){
##     return(logres);
##   }else{
##     return(exp(logres));
##   }

## }

## plot_dlogit_beta <- function(a, b){
##   x <- seq(-5, 5, by = 0.01)
##   par(mfrow = c(3, 1))
##   plot(plogis(x), dbeta(plogis(x), a, b), main = 'beta')
##   plot(x, dlogitbeta(x, a, b), main = 'logit_beta')
##   plot(plogis(x), dlogitbeta(x, a, b), main = 'scaled logit_beta')
## }

## ## plot_dlogit_beta(2, 2)
## ## plot_dlogit_beta(.5, .5)
## ## plot_dlogit_beta(2, 5)


## dlogtgaussian <- function( log_prec,  u,  s, give_log=0){

##   part1 = -0.5 * log(8 * pi) - s;
##   part2 = -0.5 * (1 / s ^ 2) *  ((exp(-log_prec / 2) - u ) ^ 2) - log_prec / 2;

##   logres = part1 + part2;

##   if(give_log){
##     return(logres);
##   }else{
##     return(exp(logres));
##   }

## }

## plot_dlogtgaussian <- function(u, s, range_x = 5){
##   x <- seq(-range_x, range_x, by = 0.01)
##   par(mfrow = c(3, 1))
##   plot(1 / sqrt(exp(x)), dnorm(1 / sqrt(exp(x)), u, s) / .5, main = 'trunc norm stdev')
##   plot(x, dlogtgaussian(x, u, s), main = 'log prec')
##   plot(1 / sqrt(exp(x)), dlogtgaussian(x, u, s), main = 'scaled log prec')
## }

## ## plot_dlogtgaussian(0, 1)
## ## plot_dlogtgaussian(3, 5, range_x = 10)
## ## plot_dlogtgaussian(0, 10, range_x = 10)


## summary function for long data
## utility function from the interwebs
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.80, .drop=TRUE) {
  library(plyr)

  ## New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     med  = median (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     l.ci = stats::quantile(xx[[col]], probs = (1-conf.interval)/2, na.rm=na.rm),
                     u.ci = stats::quantile(xx[[col]], probs = 1-(1-conf.interval)/2, na.rm=na.rm)
                     )
                 },
                 measurevar
                 )

  ## Rename the "mean" column
  old.names <- c("mean",
                 paste0('l.ci.', (1-conf.interval)/2*100, '%'),
                 paste0('u.ci.', (1-(1-conf.interval)/2)*100, '%'))
  new.names <- c(measurevar, 'l.ci', 'u.ci')
  names(datac)[match(old.names, names(datac))] <- new.names

  ## get with width of the ci
  datac$w.ci <- datac$u.ci - datac$l.ci

  return(datac)
}
