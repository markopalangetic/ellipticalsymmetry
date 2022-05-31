
# This function is used to calculate the test statistic.
# Since this test uses bootstrap, this function will be used to calculate the empirical distribution of the test statistic.
getstatisticKS = function(X){


  data_size = dim(X)
  n = data_size[1] #sample size
  d = data_size[2] #dimension

  theta = colMeans(X) # estimator of the mean
  sigma = (n - 1)*stats::cov(X)/n # estimator of the covariance matrix. Here, 1/n is used instead of 1/(n-1) as the averaging factor.
  sigma_root = spd_matrix_pow(sigma, -1/2)  #the square root of the inverse of the covariance matrix
  Z = sweep(X, 2, theta)%*%sigma_root #data standardization
  norms = apply(Z, 1, function(x) norm(x, type='2')) #norms of the standardized data
  norms_order = order(norms) #sort the norms

  # The following lines are used to construct the matrix where the evaluation of the spherical harmonics will be saved.
  deg1 = d
  deg2 = choose(d + 1, 2) - 1
  deg3 = choose(d + 2, 3) - d
  deg4 = choose(d + 3, 4) - choose(d + 1, 2)

  deg = deg1 + deg2 + deg3 + deg4
  func_eval = matrix(0, deg, n)

  # The following lines calculate the test statistic as it is described in the article.
  indexi = 1
  for (i in seq2(1, d)){
    indexj = 1
    for (w in norms_order){
      sample = Z[w,]
      sample_norm = norm(sample, type='2')
      func_eval[indexi,indexj] = func1deg1(sample/sample_norm, i)
      indexj = indexj + 1
    }
    indexi = indexi + 1
  }
  for (i in seq2(1, d)){
    for (j  in seq2(i + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')
        func_eval[indexi,indexj] = func1deg2(sample/sample_norm, i, j)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }

  for (k in seq2(2, d)){
    indexj = 1
    for (w in norms_order){
      sample = Z[w,]
      sample_norm = norm(sample, type='2')
      func_eval[indexi,indexj] = func2deg2(sample/sample_norm, k)
      indexj = indexj + 1
    }
    indexi = indexi + 1
  }


  for(i in seq2(1, d)){
    for (j  in seq2(i + 1, d)){
      for (k  in seq2(j + 1, d)){
        indexj = 1
        for (w in norms_order){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')

          func_eval[indexi,indexj] = func1deg3(sample/sample_norm, i, j, k)
          indexj = indexj + 1
        }
        indexi = indexi + 1
      }
    }
  }
  for(k in seq2(2, d)){
    for (r  in seq2(k + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')

        func_eval[indexi,indexj] = func2deg3(sample/sample_norm, k, r)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }

  for(j in seq2(1, d)){
    for (k  in seq2(j + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')

        func_eval[indexi,indexj] = func3deg3(sample/sample_norm, j, k)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }
  for (r  in seq2(2, d)){
    indexj = 1
    for (w in norms_order){
      sample = Z[w,]
      sample_norm = norm(sample, type='2')

      func_eval[indexi,indexj] = func4deg3(sample/sample_norm, r)
      indexj = indexj + 1
    }
    indexi = indexi + 1
  }

  for (r  in seq2(2, d)){
    indexj = 1
    for (w in norms_order){
      sample = Z[w,]
      sample_norm = norm(sample, type='2')

      func_eval[indexi,indexj] = func1deg4(sample/sample_norm, r)
      indexj = indexj + 1
    }
    indexi = indexi + 1
  }

  for(i in seq2(1, d)){
    for (r  in seq2(i + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')

        func_eval[indexi,indexj] = func2deg4(sample/sample_norm, i, r)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }
  for(i in seq2(1, d)){
    for (j  in seq2(i + 1, d)){
      for (r  in seq2(j + 1, d)){
        indexj = 1
        for (w in norms_order){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')

          func_eval[indexi,indexj] = func3deg4(sample/sample_norm, i, j, r)
          indexj = indexj + 1
        }
        indexi = indexi + 1
      }
    }
  }
  for(r in seq2(2, d)){
    for (s  in seq2(r + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')

        func_eval[indexi,indexj] = func4deg4(sample/sample_norm, r, s)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }
  for(i in seq2(1, d)){
    for (j  in seq2(i + 1, d)){
      for (r  in seq2(j + 1, d)){
        for (s  in seq2(r + 1, d)){
          indexj = 1
          for (w in norms_order){
            sample = Z[w,]
            sample_norm = norm(sample, type='2')

            func_eval[indexi,indexj] = func5deg4(sample/sample_norm, i, j, r, s)
            indexj = indexj + 1
          }
          indexi = indexi + 1
        }
      }
    }
  }
  for (k  in seq2(2, d)){
    for (r  in seq2(k + 1, d)){
      for (s  in seq2(r + 1, d)){
        indexj = 1
        for (w in norms_order){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')

          func_eval[indexi,indexj] = func6deg4(sample/sample_norm, k, r, s)
          indexj = indexj + 1
        }
        indexi = indexi + 1
      }
    }
  }
  for (j  in seq2(1, d)){
    for (r  in seq2(j + 1, d)){
      for (s  in seq2(r + 1, d)){
        indexj = 1
        for (w in norms_order){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')

          func_eval[indexi,indexj] = func7deg4(sample/sample_norm, j, r, s)
          indexj = indexj + 1
        }
        indexi = indexi + 1
      }
    }
  }
  for(r in seq2(2, d)){
    for (s  in seq2(r + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')

        func_eval[indexi,indexj] = func8deg4(sample/sample_norm, r, s)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }
  statistic = 0

  func_eval_sums = matrix(0, deg, n)

  for (i in seq(deg)){
    val = 0
    for (j in seq(n)){
      val = val + func_eval[i,j]
      func_eval_sums[i,j] = val^2
    }
  }

  statistic = sqrt(max(colSums(func_eval_sums))/n) # test statistic
  return(statistic)

}


# This function implements the bootstrap procedure for the calculation of the empirical distribution of the test statistic.
bootstrapKS = function(X, R, statistic, nJobs){
  data_size = dim(X)
  n = data_size[1] # sample size
  d = data_size[2] # dimension

  theta = colMeans(X) # estimator of the mean
  sigma = (n - 1)*stats::cov(X)/n # estimator of the covariance matrix. Here, 1/n is used instead of 1/(n-1) as the averaging factor.
  sigma_root = spd_matrix_pow(sigma, -1/2) #the square root of the inverse of the covariance matrix
  Z = sweep(X, 2, theta)%*%sigma_root #data standardization
  norms = apply(Z, 1, function(x) norm(x, type='2')) #norms of the standardized data
  # The following 'for' loop calculates the bootstrap replicates.
  # bootstrap_statistic = replicate(R, 0) # in this vector we save the boostrap replicates
  # for (i in seq(R)){
  #   U = unisphere(n, d=d)
  #   boot_norms = sample(norms, replace=T)
  #   Xsim = U*boot_norms
  #   bootstrap_statistic[i] = getstatisticKS(Xsim)
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
    getstatisticKS(Xsim)[[1]]
  }
  parallel::stopCluster(cl)
  bootstrap_statistic = unlist(bootstrap_statistic)

  emp_cdf = stats::ecdf(bootstrap_statistic) #we calculate the empirical cumulative distribution function (cdf) based on the bootstrap replicates
  p_val = 1 - emp_cdf(statistic) # p value
  return(p_val)
}

#' Koltchinskii and Sakhanenko's test for elliptical symmetry
#'
#' @description Test for elliptical symmetry.
#'
#' @param X A numeric matrix.
#' @param R The number of bootstrap replicates.
#' @param nJobs The number of CPU cores used for the calculation. The default value -1 indicates that all cores except one are used.
#'
#' @return An object of class \code{"htest"} containing the following components:
#'  \item{\code{statistic}}{The value of the test statistic.}
#'  \item{\code{pvalue}}{The p-value of the test.}
#'  \item{\code{alternative}}{A character string describing the alternative hypothesis.}
#'  \item{\code{method}}{A character string indicating what type of test was performed.}
#'
#' @section Background:
#' Koltchinskii and Sakhanenko (2000) proposed a class of omnibus bootstrap tests for elliptical symmetry
#' that are affine invariant and consistent against any fixed alternative. This test is based on spherical harmonics.
#'
#' @references
#' Koltchinskii, V., & Sakhanenko, L., (2000). Testing for ellipsoidal symmetry of a multivariate distribution. \emph{High Dimensional Probability II}, 493-510, Springer.
#'
#' Sakhanenko, L., (2008). Testing for ellipsoidal symmetry: A comparison study. \emph{Computational Statistics & Data Analysis}, \bold{53}(2), 565-581.
#'
#' @examples
#'
#' ## sepal width and length of the versicolor subset of the Iris data
#' X = datasets::iris[51:100,1:2]
#'
#' KoltchinskiiSakhanenko(X, R = 10, nJobs=2)
#' @export
KoltchinskiiSakhanenko = function(X, R=1000, nJobs=-1){

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


  statistic = getstatisticKS(X)[[1]] # test statistic
  p_val = bootstrapKS(X, R, statistic, nJobs)[[1]] # p value
  #The following lines construct htest object 'res' which is the output of this function.
  names(statistic) = 'statistic'
  second_part_name = paste('with bootstrap p-value (based on',R,'replicates)')
  method = paste('Test for elliptical symmetry by Koltchinskii and Sakhanenko', second_part_name)
  res <- list(method = method,
              data.name = dname,
              statistic = statistic,
              p.value = p_val,
              alternative = 'the distribution is not elliptically symmetric')
  class(res) <- "htest"
  return(res)
}
