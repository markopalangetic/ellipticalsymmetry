#' @import foreach
NULL

# This function implements an exponentiation of a symmetric and positive definite (spd) matrix.
# Parameter 'A' is a spd matrix while 'pow' is an exponent that is applied.
# The function is implemented to avoid having additional packages in the dependencies.
spd_matrix_pow = function(A, pow){
  if(mean(abs(A - t(A))) > 1e-6){
    stop("The matrix is not symmetric.")
  }

  eig = eigen(A)
  for (val in eig$values){
    if (val < 0){
      stop("The matrix is not positive semi-definite.")
    }
  }
  new_eig_vals = eig$values ^ pow
  Q = eig$vectors
  D = diag(new_eig_vals)
  res = Q %*% D %*% t(Q)
  return(res)

}


#This function implements 'seq' method that will not return anything when a > b.
seq2 = function(a,b){
  if (b >= a)  {
    return (seq(a, b))
  }
  else return()
}


#This function generates data that are uniformly distributed over a unit d-dimensional sphere.
#Parameter 'n' stands for the number of samples, while 'd' is a dimension.
#The function first generates n samples from the standard multivariate normal distribution of dimension d and then divide every sample instance with its norm.
# The function is implemented to avoid having additional packages in the dependencies.
unisphere = function(n, d){
  points = matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  norms = apply(points, 1, function(x) norm(x, type='2'))
  points = points/norms
  return (points)
}


#This function calculates the tyler shape estimator from the package ICSNP
#and modifies it in accordance with the definition of F_1 from [1].


#[1] Babic, S., Gelbgras, L., Hallin, M., & Ley, C. (2019). Optimal tests for elliptical symmetry: specified and unspecified location. arXiv preprint arXiv:1911.08171.
tyler_cov = function(X, location) {
  n = dim(X)[1]
  d = dim(X)[2]
  sigma = ICSNP::tyler.shape(X, location = location)

  #b = array(0, dim = c(1, n))

  #sigma_inv = solve(sigma)

  #for (i in 1:n) {
   # b[i] = (X[i, ] - location) %*% sigma_inv %*% (X[i, ] - location)
  #}
  sigma_root = spd_matrix_pow(sigma, -1/2)
  Z = sweep(X, 2, location)%*%sigma_root
  norms_sq = apply(Z*Z, 1, function(x) sum(x))
  sigma = sigma * mean(norms_sq) / d
  return(sigma)
}

