#The following two functions are auxilary functions to calculate phi_t.
f = function(v, k, x) {
  val = (1 + x * x / v) ^ (-(v + k) / 2)
  return(val)
}

fprime = function(v, k, x) {
  # derivative of f
  val = (-(v + k) / 2) * (1 + x * x / v) ^ (-(v + k) / 2 - 1) * 2 * x /
    v
  return(val)
}

#The following three functions are neccessary to calculate the test statistic when "f = 't'".
phi_t = function(v, k, x)
  # phi = -f'/f
  return(-fprime(v, k, x) / f(v, k, x))

phiprime_t = function(v, k, x) {
  val = (v + k) * (v - x * x) / (v + x * x) ^ 2
  return(val)
}

getK_t = function(v, d, norms, n) {

  summ = 0
  for (i in seq(1, n)) {
    summ = summ + phiprime_t(v, d, norms[i]) + (d - 1) * phi_t(v, d, norms[i]) / norms[i]
  }
  summ = summ / n
  return(summ)
}

#The following three functions are neccessary to calculate the test statistic when "f = 'logistic'".
phi_logistic = function(x){
  val = 2*x*(1-exp(-x^2))/(1+exp(-x^2))
  return(val)
}


phiprime_logistic = function(x) {
  val = (2*(1-exp(-x^2))*(1+exp(-x^2)) + 8*x^2*exp(-x^2))/(1+exp(-x^2))^2
  return(val)
}


getK_logistic = function(d, norms, n) {

  summ = 0
  for (i in seq(1, n)) {
    summ = summ + phiprime_logistic(norms[i]) + (d - 1) * phi_logistic(norms[i]) / norms[i]
  }
  summ = summ / n
  return(summ)
}

#The following three functions are neccessary to calculate the test statistic when "f = 'phi_powerExp'".
phi_powerExp = function(lambda, x){
  val = lambda*x^(2*lambda -1)
  return(val)
}


phiprime_powerExp = function(lambda, x) {
  val = lambda*(2*lambda - 1)*x^(2*lambda - 2)
  return(val)
}


getK_powerExp = function(lambda, d, norms, n) {
  summ = 0
  for (i in seq(1, n)) {
    summ = summ + phiprime_powerExp(lambda, norms[i]) + (d - 1) * phi_powerExp(lambda, norms[i]) / norms[i]
  }
  summ = summ / n
  return(summ)
}


 #This function implements Skew Optimal test when the location is specified.
SkewOptimalLK = function(X, location){
  data_size = dim(X)
  n = data_size[1] #sample size
  d = data_size[2] #dimension


  theta = colMeans(X)  #estimator of the mean
  sigma = tyler_cov(X, location) #tyler estimator of the convariance matrix
  sigma_inv = solve(sigma)  #inverse of the estimator of the convariance matrix


  # This line calculates the test statistic as it is described in the article
  statistic = n*t(theta - location)%*%sigma_inv%*%(theta- location)

  p_val = 1 - stats::pchisq(statistic, df = d) # p value
  output = list(statistic = statistic[[1]], p.value = p_val[[1]])
  return(output)
}


#This function implements Skew Optimal test when the location is unspecified.
SkewOptimalLU = function(X, f = 't', param = NA) {

  data_size = dim(X)
  n = data_size[1] #sample size
  d = data_size[2] #dimension

  theta = colMeans(X) #estimator of the mean
  sigma = tyler_cov(X, theta) #tyler estimator of the convariance matrix

  sigma_root = spd_matrix_pow(sigma, -1/2) #the square root of the inverse of the covariance matrix
  Z = sweep(X, 2, theta)%*%sigma_root #data standardization
  #d = apply(Z, 1, norm, type = "2")
  norms = apply(Z, 1, function(x) norm(x, type='2')) #norms of the standardized data
  U = Z / norms # normalization of the standardized data


  # The following lines calculate the test statistic as it is described in the article.
  # We distingusih three cases for parameter 'f'.
  idmatr = diag(d)  # identity matrix
  gama = matrix(0, nrow = d, ncol = d)
  delta = matrix(0, nrow = d, ncol = d)
  if(f == 't'){
    if (any(is.na(param))){
      param = 4
    }
    J = getK_t(param, d, norms, n)
    for (i in 1:n) {
      phi_val = phi_t(param, d, norms[i])
      gama = gama + (norms[i] - d * phi_val / J) ^ 2 / d * idmatr
      delta = delta + (norms[i] - d * phi_val / J) * U[i,]
    }
  }

  else if(f == 'logistic'){
    if (!any(is.na(param))){
      warning('param is specified but it will not be used')
    }
    J = getK_logistic(d, norms, n)
    for (i in 1:n) {
      phi_val = phi_logistic(norms[i])
      gama = gama + (norms[i] - d * phi_val / J) ^ 2 / d * idmatr
      delta = delta + (norms[i] - d * phi_val / J) * U[i, ]
    }
  }

  else if(f == 'powerExp'){
    if (any(is.na(param))){
      param = 0.5
    }
    J = getK_powerExp(param, d, norms, n)
    for (i in 1:n) {
      phi_val = phi_powerExp(param, norms[i])
      gama = gama + (norms[i] - d * phi_val / J) ^ 2 / d * idmatr
      delta = delta + (norms[i] - d * phi_val / J) * U[i, ]
    }
  }
  else{
    stop('f has to take one of the following three values: t, logistic, powerExp')
  }

  statistic = t(delta) %*% solve(gama) %*% delta # test statistic
  p_val = 1 - stats::pchisq(statistic, df = d) # p value
  output = list(statistic = statistic[[1]], p.value = p_val[[1]])
  return(output)
}


#' Tests for elliptical symmetry by Babic et al.
#'
#' @description Tests for elliptical symmetry: specified and unspecified location.
#'
#' @param X A numeric matrix.
#' @param location A vector of location parameters.
#' @param f A string that specifies the type of the radial density upon which the test is based. Currently supported options are \code{"t"}, \code{"logistic"} and \code{"powerExp"}.
#' The default is set to \code{"t"}.
#' @param param A parameter that is used when \code{f = "t"} and \code{f = "powerExp"}.
#' The default value of \code{param} represents the degrees of freedom of the multivariate t distribution and it is set to 4.
#'
#' @return An object of class \code{"htest"} containing the following components:
#'  \item{\code{statistic}}{The value of the test statistic.}
#'  \item{\code{pvalue}}{The p-value of the test.}
#'  \item{\code{alternative}}{A character string describing the alternative hypothesis.}
#'  \item{\code{method}}{A character string indicating what type of test was performed.}
#'
#' @details
#' \code{X} and \code{location} are the only input arguments for the specified location test.
#' The default value for \code{location} is set to \code{NA} which implies that the unspecified location test will be performed
#' unless the user specifies location.
#'
#' For the unspecified location test, besides the data matrix \code{X}, the input arguments are \code{f} and \code{param}.
#' The \code{f} argument is a string that specifies the type of the radial density upon which the test is based.
#' Currently supported options are: \code{"t"} for the radial density of the multivariate t distribution,
#' \code{"logistic"} for the multivariate logistic and \code{"powerExp"} for the radial density of the multivariate power-exponential distribution.
#' Note that the default is set to \code{"t"}.
#' The role of the \code{param} argument is as follows.
#' If \code{f = "t"} then \code{param} denotes the degrees of freedom of the multivariate t distribution.
#' Given that the default radial density is \code{"t"}, it follows that the default value of \code{param}
#' represents the degrees of freedom of the multivariate t distribution and it is set to 4.
#' Note also that the degrees of freedom have to be greater than 2.
#' If \code{f = "powerExp"} then \code{param} denotes the kurtosis parameter. In that case the value of \code{param}
#' has to be different from 1, because for the multivariate power exponential distribution, a kurtosis parameter equal to 1 corresponds
#' to the multivariate Gaussian distribution (the Gaussian \code{f} is excluded due to a singular Fisher information matrix).
#' The default value is set to 0.5.
#'
#' @section Background:
#' Tests for elliptical symmetry both for specified and unspecified location. These tests are based on
#' Le Cam’s theory of statistical experiments and they are optimal against generalized skew-elliptical alternatives,
#' but they remain quite powerful under a much broader class of non-elliptical distributions.
#' They have a simple asymptotic chi-squared distribution under the null hypothesis of ellipticity,
#' they are affine-invariant, computationally fast, have a simple and intuitive form, only require finite moments of order 2.
#'
#'
#'
#' @references
#' Babic, S., Gelbgras, L., Hallin, M., & Ley, C. (2021). Optimal tests for elliptical symmetry: specified and unspecified location. Bernoulli (in press).
#'
#'
#' @examples
#'
#' ## sepal width and length of the versicolor subset of the Iris data
#' X = datasets::iris[51:100,1:2]
#'
#'
#' ## location unspecified test based on the radial density of the multivariate t distribution
#' SkewOptimal(X)
#'
#' ## location unspecified test based on the radial density of the logistic distribution
#' SkewOptimal(X, f="logistic")
#'
#' ## location unspecified test based the radial density of the power exponential distribution
#' SkewOptimal(X, f="powerExp")
#'
#' @export
SkewOptimal = function(X, location = NA, f = 't', param = NA) {

  dname = deparse(substitute(X)) # getting the data name

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
  if (any(is.na(location))){
    output = SkewOptimalLU(X, f, param)
  }
  else{
    output = SkewOptimalLK(X, location)
  }
  names(output$statistic) = 'statistic'

  #The following lines construct htest object 'res' which is the output of this function.
  res <- list(method = 'SkewOptimal test for elliptical symmetry',
              data.name = dname,
              statistic = output$statistic,
              p.value = output$p.value,
              alternative = 'the distribution is not elliptically symmetric')
  class(res) <- "htest"
  return(res)
}



