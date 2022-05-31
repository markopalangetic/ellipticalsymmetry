#This function caluclates the moments of degree 4 using Kronecker (tensor) product. See the article for more details.
getM4 <- function(x, m) {
  scal_prod = (x - m) %*% t(x - m)
  kron_prod = scal_prod %x% scal_prod
  return(kron_prod)
}

#' Schott's test for elliptical symmetry
#'
#' @description Test for elliptical symmetry.
#'
#' @param X A numeric matrix.
#' @return An object of class \code{"htest"} containing the following components:
#'  \item{\code{statistic}}{The value of the test statistic.}
#'  \item{\code{pvalue}}{The p-value of the test.}
#'  \item{\code{alternative}}{A character string describing the alternative hypothesis.}
#'  \item{\code{method}}{A character string indicating what type of test was performed.}
#'
#'
#' @section Background:
#' A Wald-type test for elliptical symmetry based on fourth moments.
#' It compares the sample fourth moments with the expected theoretical ones under ellipticity.
#' Being based on fourth-order moments, the test is very simple to use but requires moments of order 8.
#' It has an asymptotic chi-squared distribution under the null hypothesis of ellipticity.
#'
#' @references
#' Schott, James R., (2002). Testing for elliptical symmetry in covariance-matrix-based analyses. \emph{Statistics & Probability Letters}, \bold{60}(4), 395-404.
#' @examples
#'
#' ## sepal width and length of the versicolor subset of the Iris data
#' X = datasets::iris[51:100,1:2]
#'
#' Schott(X)
#'
#' @export
Schott = function(X) {

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

  data_size = dim(X)
  n = data_size[1] #sample size
  d = data_size[2] #dimension
  theta = colMeans(X) #estimator of the mean
  sigma = stats::cov(X) #the unbiased estimator of the convariance matrix
  sigma_root = spd_matrix_pow(sigma, -1/2) #the square root of the inverse of the covariance matrix

  Z = sweep(X, 2, theta)%*%sigma_root #data standardization

  # The following lines calculate the test statistic as it is described in the article.
  M4 = matrix(0, nrow = d * d, ncol = d * d)
  for (i in seq(1, n)) {
    M4 = M4 + getM4(X[i, ], theta)
  }
  M4 = 1 / n * M4

  sigma_root = spd_matrix_pow(sigma, -1/2)
  sigma_root_kron_t = t(sigma_root) %x% t(sigma_root)
  sigma_root_kron = sigma_root %x% sigma_root
  M4_stand = sigma_root_kron_t %*% M4 %*% sigma_root_kron

  norms_sq = apply(Z*Z, 1, function(x) sum(x))



  kapa = sum(norms_sq^2) / (n * d * (d + 2))
  eta = sum(norms_sq^3) / (n * d * (d + 2) * (d + 4))
  omega = sum(norms_sq^4) / (n * d * (d + 2) * (d + 4) * (d + 6))

  beta1 = omega ^ (-1) / 24
  a = omega + kapa ^ 3 - 2 * kapa * eta
  beta2 = -3 * a * (24 * omega ^ 2 + 12 * (d + 4) * a * omega) ^ (-1)

  M4_stand_sq = M4_stand %*% M4_stand
  vec_id = as.vector(diag(d))
  M4_trace = sum(diag(M4_stand_sq))

  statistic = n * (beta1 * M4_trace + beta2 * t(vec_id) %*% M4_stand_sq %*% vec_id -
                      (3 * beta1 + (d + 2) * beta2) * d * (d + 2) * kapa ^ 2) #test statistic

  statistic = statistic[[1]]
  df = d ^ 2  + d * (d - 1) * (d ^ 2 + 7 * d - 6) / 24 - 1 #This line calculates the degree of freedom of the chi-squared distribution.
  p_val = 1 - stats::pchisq(statistic, df = df)[[1]] #p value
  #The following lines construct htest object 'res' which is the output of this function.
  names(statistic) = 'statistic'
  res <- list(method = 'Schott test for elliptical symmetry',
              data.name = dname,
              statistic = statistic,
              p.value = p_val,
              alternative = 'the distribution is not elliptically symmetric')
  class(res) <- "htest"
  return(res)
}

