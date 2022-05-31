
#This function implements Pseudo-Gaussian test when the location is specified.
PseudoGaussianLK = function(X, location) {
  data_size = dim(X)
  n = data_size[1] #sample size
  d = data_size[2] #dimension

  sigma = tyler_cov(X, location) #tyler estimator of the convariance matrix
  sigma_root = spd_matrix_pow(sigma,-1 / 2) #the square root of the inverse of the covariance matrix
  Z = sweep(X, 2, location) %*% sigma_root #data standardization
  norms = apply(Z, 1, norm, type = "2") #norms of the standardized data
  # The following lines calculate the test statistic as it is described in the article.
  U = Z / norms
  Z_sign_sq = array(0, dim = c(n, d))

  Z_sign_sq = Z*Z*sign(Z)
  ones = rep(1,n)
  M4 = sum(norms^4)/n

  const = d * (d + 2) / (3 * n * M4)
  statistic = const * t(ones)%*%Z_sign_sq%*%t(Z_sign_sq)%*%ones #test statistic
  p_val = 1 - stats::pchisq(statistic, df = d) #p value
  output = list(statistic = statistic[[1]], p.value = p_val[[1]])
  return(output)
}

#This function implements Pseudo-Gaussian test when the location is unspecified.
PseudoGaussianLU = function(X) {

  data_size = dim(X)
  n = data_size[1] #sample size
  d = data_size[2] #dimension

  theta = colMeans(X) #estimator of the mean
  sigma = tyler_cov(X, theta) #tyler estimator of the convariance matrix
  sigma_root = spd_matrix_pow(sigma, -1/2) #the square root of the inverse of the covariance matrix

  Z = sweep(X, 2, theta)%*%sigma_root #data standardization
  norms = apply(Z, 1, norm, type = "2") #norms of the standardized data
  U = Z/norms # normalization of the standardized data

  # The following lines calculate the test statistic as it is described in the article.
  M1 = sum(norms)/n
  M2 = sum(norms^2)/n
  M3 = sum(norms^3)/n
  M4 = sum(norms^4)/n

  Si = Z*Z*sign(Z)

  ck = 4 * gamma(d / 2) / ((d ^ 2 - 1) * sqrt(pi) * gamma((d - 1) / 2))

  ones = rep(1, n)
  centseq = ck * (d + 1) * M1 * t(ones)%*%Z - t(ones)%*%Si
  centseq = centseq / sqrt(n)

  Gam = 1/(3 * M4 / (d * (d + 2)) - 2 * ck ^ 2 * (d + 1) * M1 * M3 +
           ck ^ 2 * (d + 1) ^ 2 / d * M1 ^ 2 * M2) * diag(d)

  statistic = centseq %*% Gam %*% t(centseq) #test statistic
  p_val = 1 - stats::pchisq(statistic, df = d) #p value
  output = list(statistic = statistic[[1]], p.value = p_val[[1]])
}

#' Pseudo-Gaussian test for elliptical symmetry
#'
#' Tests for elliptical symmetry: specified and unspecified location.
#'
#' @param X A numeric matrix.
#' @param location A vector of location parameters.
#'
#' @return An object of class \code{"htest"} containing the following components:
#'  \item{\code{statistic}}{The value of the test statistic.}
#'  \item{\code{pvalue}}{The p-value of the test.}
#'  \item{\code{alternative}}{A character string describing the alternative hypothesis.}
#'  \item{\code{method}}{A character string indicating what type of test was performed.}
#'
#' @details
#' Note that \code{location} allows the user to specify the known location.
#' The default is set to \code{NA} which means that the unspecified location test will be performed unless the user specifies location.
#'
#' @section Background:
#' Pseudo-Gaussian tests for elliptical symmetry are based on Le Cam’s theory of statistical experiments.
#' They are most efficient against a multivariate form of Fechner-type asymmetry.
#' These tests require finite moments of order 4 and they have a simple asymptotic chi-squared distribution
#' under the null hypothesis of ellipticity.
#'
#'
#' @references
#' Cassart, D., Hallin, M. & Paindaveine, D., (2008). Optimal detection of Fechner-asymmetry. \emph{Journal of Statistical Planning and Inference}, \bold{138}, 2499-2525.
#'
#' Cassart, D., (2007). Optimal tests for symmetry. Ph.D. thesis, Univ. libre de Bruxelles, Brussels.
#'
#' @examples
#'
#' ## sepal width and length of the versicolor subset of the Iris data
#' X = datasets::iris[51:100,1:2]
#'
#' PseudoGaussian(X)
#'
#' @export
PseudoGaussian = function(X, location = NA) {


  dname = deparse(substitute(X)) # get the data name

  # The following condition cheks if data have the matrix form. If not, it tries to convert data into a matrix if possible.
  if(!is.matrix(X)) {
    X = as.matrix(X)
    if (!(is.matrix(X) && length(X) > 1)){
      stop("X is not in the valid matrix form.")
    }
  }

  # The following condition checks if all matrix instances are numeric values.
  else if(!is.numeric(X)){
    stop('X has to take numeric values')
  }

  # The following condition checks if the location is specified.
  if (any(is.na(location))) {
    output = PseudoGaussianLU(X)
  }
  else{
    output = PseudoGaussianLK(X, location)
  }
  #The following lines construct htest object 'res' which is the output of this function.
  names(output$statistic) = 'statistic'
  res <-
    list(
      method = 'Pseudo-Gaussian test for elliptical symmetry',
      data.name = dname,
      statistic = output$statistic,
      p.value = output$p.value,
      alternative = 'the distribution is not elliptically symmetric'
    )
  class(res) <- "htest"
  return(res)

}
