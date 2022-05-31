
# This function calculates the angle on the interval [0, 2*pi] that point (x, y) constructs with the positive x axis.
getAngle = function(x, y){
  val = atan2(y, x)
  if(val < 0){
    val = 2*pi + val
  }
  return(val)
}

# This function calculates the order of the permutation that is given as an input parameter
getPermutationIndex = function(permutations){
  perm_string = paste(permutations, collapse = '')
  string_size=nchar(perm_string)
  if(string_size == 1){
    return(1)
  }
  ind=0
  for (i in seq(1, string_size)){
    if (substr(perm_string, i, i) == paste(string_size)){
      ind=i
      break
    }
  }
  new_perm = paste(substr(perm_string, 1, ind-1), substr(perm_string, ind + 1, string_size), sep = '')
  return((ind - 1)*factorial(string_size - 1) + getPermutationIndex(new_perm))
}

#This function calculates p value in the case when the exact distribution is known. The exact distribution is a linear
#combination of three chi-squared distributions. The degrees of freedom of those three coefficients are 'df1', 'df2', 'df3'
#while the linear weights are '1', '1-astar', '1-bstar'. Since the explicit distribution of such linear combination is unkown,
#we sampled 1000000 times from the three chi-squared distribution, calculated linear combination of the obtained samples and
#calculated empirical distribution from the obtained combinations.
getExactPValue = function(d, c, statistic){
  cutoffs = stats::qchisq(seq(0, 1, length.out = c + 1), df = d)
  a = replicate(c, 0)
  b = replicate(c, 0)
  for (i in seq(1, c)){
    a[i] = stats::pchisq(cutoffs[i + 1], df = d + 1) - stats::pchisq(cutoffs[i] , df = d + 1)
    b[i] = stats::pchisq(cutoffs[i + 1], df = d + 2) - stats::pchisq(cutoffs[i] , df = d + 2)
  }
  astar = 2*c*sum(a^2)/pi
  bstar = 4*c*sum(b^2)/(pi*pi)
  df1 = c*(2^d - 1) - d*(d + 1)/2
  df2 = d
  df3 = d*(d - 1)/2
  sample_size = 1000000
  sample1 = stats::rchisq(sample_size, df1)
  sample2 = stats::rchisq(sample_size, df2)
  sample3 = stats::rchisq(sample_size, df3)
  sample = sample1 + (1 - astar)*sample2 + (1 - bstar)*sample3
  emp_cdf = stats::ecdf(sample)
  p_val = 1 - emp_cdf(statistic)
  return(p_val)
}

#This function transforms an array of 0 and 1 into an integer for which the array is the bit representation.
bitsToInt = function(x) {
  packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
}

# This function is used to calculate the test statistic.
# Since this test uses bootstrap, this function will be used to calculate the empirical distribution of the test statistic.
getstatisticHP = function(X, g, c, sector){
  data_size = dim(X)
  n = data_size[1] #sample size
  d = data_size[2] #dimension

  theta = colMeans(X) # estimator of the mean
  sigma = (n - 1)*stats::cov(X)/n # estimator of the covariance matrix. Here, 1/n is used instead of 1/(n-1) as the averaging factor.
  L = solve(chol(sigma)) # Inverse of the cholesky decomposition of 'sigma'. L is an upper triangular matrix for which 'L %*% t(L) = sigma'.
  Z = sweep(X, 2, theta)%*%L #data standardization
  norms = apply(Z, 1, function(x) norm(x, type='2')) #norms of the standardized data
  emp_cdf = stats::ecdf(norms) #empirical cumulative distribution function of the norms.

  # The following lines calculate the test statistic as it is described in the article.
  # We distingusih three cases for parameter 'sector'.
  if (sector == 'orthants'){
    counts = matrix(0, nrow=2 ^ d, ncol = c)
    for (i in seq(n)){
      sample = Z[i,]
      bin = (sign(sample) + 1)/2
      ind1 = bitsToInt(bin) + 1
      sample_norm = norm(sample, type='2')
      ind2 = min(floor(emp_cdf(sample_norm)*c) + 1, c)
      counts[ind1, ind2] = counts[ind1, ind2] + 1

    }
    expected = n/(2^d*c)
  }
  if (sector == 'bivariateangles'){
    counts = matrix(0, nrow=g, ncol = c)
    for (i in seq(n)){
      sample = Z[i,]
      angle = getAngle(sample[1], sample[2])
      ind1 = min(floor(angle/(2*pi)*g) + 1, g)
      sample_norm = norm(sample, type='2')
      ind2 = min(floor(emp_cdf(sample_norm)*c) + 1, c)
      counts[ind1, ind2] = counts[ind1, ind2] + 1
    }
    expected = n/(g*c)
  }
  if (sector == 'permutations'){
    fact = factorial(d)
    counts = matrix(0, nrow=fact, ncol = c)
    for (i in seq(n)){
      sample = Z[i,]
      ind1 = getPermutationIndex(order(sample))
      sample_norm = norm(sample, type='2')
      ind2 = min(floor(emp_cdf(sample_norm)*c) + 1, c)
      counts[ind1, ind2] = counts[ind1, ind2] + 1
    }
    expected = n/(fact*c)
  }
  statistic = 0
  for (count in counts){
    statistic = statistic + (expected - count)^2/expected
  }
  return(statistic)
}

# This function implements the bootstrapping procedure for the calculation of the empirical distribution of the test statistic.
bootstrapHP = function(X, g, c, sector, R, statistic, nJobs){
  data_size = dim(X)
  n = data_size[1] #sample size
  d = data_size[2] #dimension

  theta = colMeans(X) # estimator of the mean
  sigma = (n - 1)*stats::cov(X)/n # estimator of the covariance matrix. Here, 1/n is used instead of 1/(n-1) as the averaging factor.
  L = solve(chol(sigma)) # Inverse of the cholesky decomposition of 'sigma'. L is an upper triangular matrix for which 'L %*% t(L) = sigma'.
  norms = apply(sweep(X, 2, theta)%*%L, 1, function(x) norm(x, type='2')) #norms of the standardized data

  # bootstrap_statistic = replicate(R, 0) # The following 'for' loop calculates the bootstrap replicates.
  # for (i in seq(R)){
  #   U = unisphere(n, d=d)
  #   boot_norms = sample(norms, replace=T)
  #   Xsim = U*boot_norms
  #   bootstrap_statistic[i] = getstatisticHP(Xsim, g, c, sector)
  # }

  if (nJobs == -1){
    no_cores = parallel::detectCores() - 1
  }
  else{
    no_cores <- min(parallel::detectCores() - 1, nJobs)
  }
  cl <- parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(cl)
  `%op%` = doRNG::`%dorng%`
  #`%op%` = foreach::`%dopar%`
  bootstrap_statistic <- foreach::foreach(i=1:R) %op% {
    U = unisphere(n, d=d)
    boot_norms = sample(norms, replace=T)
    Xsim = U*boot_norms
    getstatisticHP(Xsim, g, c, sector)[[1]]
  }
  parallel::stopCluster(cl)
  bootstrap_statistic = unlist(bootstrap_statistic)


  emp_cdf = stats::ecdf(bootstrap_statistic)#we calculate the empirical cumulative distribution function (cdf) based on the bootstrap replicates
  p_val = 1 - emp_cdf(statistic)#p value
  return(p_val)
}

#' Huffer and Park's test for elliptical symmetry
#'
#' @description Pearson chi-square type test for elliptical symmetry.
#'
#' @param X A numeric matrix.
#' @param c The number of spherical shells that are used to divide the space.
#' @param R The number of bootstrap replicates.
#' @param sector A string that specifies the type of sectors used to divide the space. Currently supported options are \code{"orthants"}, \code{"permutations"} and \code{"bivariateangles"}.
#' @param g A parameter that is used if \code{sector = "bivariateangles"}. It denotes the number of regions used to divide the plane.
#' @param nJobs The number of CPU cores used for the calculation. The default value -1 indicates that all cores except one are used.
#'
#'
#' @return An object of class \code{"htest"} containing the following components:
#'  \item{\code{statistic}}{The value of the test statistic.}
#'  \item{\code{pvalue}}{The p-value of the test.}
#'  \item{\code{alternative}}{A character string describing the alternative hypothesis.}
#'  \item{\code{method}}{A character string indicating what type of test was performed.}
#'
#' @details
#' Huffer and Park (2007) propose a Pearson chi-square type test with multi-dimensional cells.
#' After dividing the space into \code{c} spherical shells  and \code{g} sectors (in total \code{gc} cells),
#' and after determining the observed cell counts, the test statistic is easily computed.
#' \code{sector} is an option that allows the user to specify the type of sectors used to divide the space.
#' Currently supported options are \code{"orthants"}, \code{"permutations"} and \code{"bivariateangles"},
#' the last one being available only in dimension 2. The \code{g} argument indicates the number of sectors.
#' The user has to choose \code{g} only if \code{sector = "bivariateangles"} and it denotes the number of regions used to divide the plane.
#' In this case, regions consist of points whose angle in polar coordinates is between \eqn{2(m-1)\pi/g} and \eqn{2m\pi/g} for \eqn{m} in \eqn{(1,..., g)}.
#' If \code{sector} is set to \code{"orthants"}, then \code{g} is fixed and equal to \eqn{2^d}, while for \code{sector = "permutations"} \code{g} is \eqn{d}!.
#' No matter what type of sectors is chosen, the user has to specify the number of spherical shells that are used to divide the space, which is \code{c}.
#' The value of \code{c} should be such that the average cell counts \eqn{n/(gc)} are not too small.
#'
#' The asymptotic distribution is available only under \code{sector = "orthants"} when the underlying distribution is close to normal.
#' Otherwise, bootstrap procedures are required and the user can freely choose the number of bootstrap replicates, denoted as \code{R}.
#' Note that by default \code{sector} is set to \code{"orthants"} and \code{R = NA}, which means that the non-bootstrap version of the test
#' will be performed unless the user specifies \code{R}.
#'
#' @references
#' Huffer, Fred W., & Park, C., (2007). A test for elliptical symmetry. \emph{Journal of Multivariate Analysis}, \bold{98}(2), 256-281.
#'
#' @examples
#'
#'
#' ## sepal width and length of the versicolor subset of the Iris data
#' X = datasets::iris[51:100,1:2]
#'
#' ## the non-bootstrap test
#' HufferPark(X, c = 2)
#'
#' ## the bootstrap tests
#' HufferPark(X, c = 2, R = 10, sector="orthants", nJobs=2)
#'
#' HufferPark(X, c = 2, R = 10, sector="bivariateangles", g = 3, nJobs=2)
#'
#' HufferPark(X, c = 2, R = 10, sector="permutations", nJobs=2)
#' @export
HufferPark = function(X, c, R = NA, sector='orthants', g=NA, nJobs = -1){


  #The following condition corrects if one puts 'R = NA' (the version that do not use bootstrap) and 'sector' different than 'orthants'
  if (is.na(R) && sector != 'orthants'){
    warning("When the non-bootstrap test is used, the value of parameter sector different than 'orthants' is not supported. sector = 'orthants' is set.")
    sector = 'orthants'
  }
  #The following condition checks if 'sector' takes one of the three possible values. This case applies if the bootstrap is used.
  else if (!sector %in% c('orthants', 'bivariateangles', 'permutations')){
    stop('when the bootstrap is chosen (R is not NA), sector has to take one of the following three values: orthants, bivariateangles, permutations')
  }

  #The following condition checks if 'g' is specified when 'sector == 'bivariateangles''.
  if (sector == 'bivariateangles' && is.na(g)){
    stop('g has to be specified.')
  }

  #The following condition warns if 'g' is specified but 'sector != 'bivariateangles'' since then 'g' will not be used.
  else if (sector != 'bivariateangles' && !is.na(g)){
    warning('g is specified but it will not be used.')
  }

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

  if (nJobs%%1 != 0){
    stop('nJobs must take an integer value')
  }
  if (nJobs < -1 || nJobs == 0){
    stop('nJobs must be either -1 or a positive integer')
  }

  statistic = getstatisticHP(X, g, c, sector)[[1]] #test statistic

  d = dim(X)[2]#dimension

  #The following two conditions check if R is specified or not. That is, if the non-bootstrap version will be used or not.
  second_part_name = ''
  if (is.na(R)){
    p_val = getExactPValue(d, c, statistic)[[1]]
    warning("Non-bootstrap tests run on a single core. nJobs parameter will not be considered")
  }
  else{
    p_val = bootstrapHP(X, g, c, sector, R, statistic, nJobs)[[1]]
    second_part_name = paste(second_part_name, 'with bootstrap p-value (based on',R,'replicates)')
  }
  #The following lines construct htest object 'res' which is the output of this function.
  method = paste('Test for elliptical symmetry by Huffer and Park', second_part_name)
  names(statistic) = 'statistic'
  res <- list(method = method,
              data.name = dname,
              statistic = statistic,
              p.value = p_val,
              alternative = 'the distribution is not elliptically symmetric')
  class(res) <- "htest"
  return(res)
}


