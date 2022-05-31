
#' Test for elliptical symmetry by A. Manzotti, F. J. Pérez, and A. J. Quiroz.
#'
#' @description Test for elliptical symmetry based on the averages of spherical harmonics.
#'
#' @param X A numeric matrix.
#' @param epsilon A value which indicates the proportion of points close to the origin which will not be used in the calculation. Default is set to 0.05.
#'
#' @return An object of class \code{"htest"} containing the following components:
#'  \item{\code{statistic}}{The value of the test statistic.}
#'  \item{\code{pvalue}}{The p-value of the test.}
#'  \item{\code{alternative}}{A character string describing the alternative hypothesis.}
#'  \item{\code{method}}{A character string indicating what type of test was performed.}
#'
#'
#' @section Background:
#' Test for elliptical symmetry based on spherical harmonics of degrees 3 and 4. The test statistic involves the averages of spherical harmonics
#' over the projections of the scale residuals on the unit sphere. This test requires moments of order 4, it has a simple asymptotic distribution
#' but does not have asymptotic power 100% against all alternatives.
#'
#' @references
#' Manzotti, A., Perez, Francisco J., & Quiroz, Adolfo J. (2002). A statistic for testing the null hypothesis of elliptical symmetry. \emph{Journal of Multivariate Analysis}, \bold{81}(2), 274-285.
#'
#' @examples
#'
#' ## sepal width and length of the versicolor subset of the Iris data
#' X = datasets::iris[51:100,1:2]
#'
#' MPQ(X)
#'
#' @export
MPQ = function(X, epsilon = 0.05){

  dname = deparse(substitute(X)) # get the data name

  # The following condition cheks if data have the matrix form and if not, it converts data into a matrix if possible
  if(!is.matrix(X)) {
    X = as.matrix(X)
    if (!(is.matrix(X) && length(X) > 1)){
      stop("X is not in the valid matrix form.")
    }
  }

  #The following condition checks if all matrix instances are numeric values
  else if(!is.numeric(X)){
    stop('X has to take numeric values')
  }

  data_size = dim(X)
  n = data_size[1] #sample size
  d = data_size[2] #dimension

  theta = colMeans(X) #estimator of the mean
  sigma = stats::cov(X) #the unbiased estimator of the convariance matrix
  sigma_root = spd_matrix_pow(sigma, -1/2)  #the square root of the inverse of the covariance matrix
  Z = sweep(X, 2, theta)%*%sigma_root #data standardization
  norms = apply(Z, 1, function(x) norm(x, type='2')) #norms of the standardized data
  rho = stats::quantile(norms, epsilon)[[1]] #This is the treshold that is used to eliminate the norms that are smaller than it.

  # The following lines calculate the test statistic as it is described in the article.
  statistic = 0
  for(i in seq2(1, d)){
    for (j  in seq2(i + 1, d)){
      for (k  in seq2(j + 1, d)){
        tmp_statistic = 0
        for (w in seq(n)){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')
          if (sample_norm >= rho){
            tmp_statistic = tmp_statistic + func1deg3(sample/sample_norm, i, j, k)
          }
        }
        statistic = statistic + (tmp_statistic/n)^2
      }
    }
  }
  for(k in seq2(2, d)){
    for (r  in seq2(k + 1, d)){
      tmp_statistic = 0
      for (w in seq(n)){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')
        if (sample_norm >= rho){
          tmp_statistic = tmp_statistic + func2deg3(sample/sample_norm, k, r)
        }
      }
      statistic = statistic + (tmp_statistic/n)^2

    }
  }

  for(j in seq2(1, d)){
    for (k  in seq2(j + 1, d)){
      tmp_statistic = 0
      for (w in seq(n)){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')
        if (sample_norm >= rho){
          tmp_statistic = tmp_statistic + func3deg3(sample/sample_norm, j, k)
        }
      }
      statistic = statistic + (tmp_statistic/n)^2

    }
  }
  for (r  in seq2(2, d)){
    tmp_statistic = 0
    for (w in seq(n)){
      sample = Z[w,]
      sample_norm = norm(sample, type='2')
      if (sample_norm >= rho){
        tmp_statistic = tmp_statistic + func4deg3(sample/sample_norm, r)
      }
    }
    statistic = statistic + (tmp_statistic/n)^2
  }

  for (r  in seq2(2, d)){
    tmp_statistic = 0
    for (w in seq(n)){
      sample = Z[w,]
      sample_norm = norm(sample, type='2')
      if (sample_norm >= rho){
        tmp_statistic = tmp_statistic + func1deg4(sample/sample_norm, r)
      }
    }
    statistic = statistic + (tmp_statistic/n)^2
  }

  for(i in seq2(1, d)){
    for (r  in seq2(i + 1, d)){
      tmp_statistic = 0
      for (w in seq(n)){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')
        if (sample_norm >= rho){
          tmp_statistic = tmp_statistic + func2deg4(sample/sample_norm, i, r)
        }
      }
      statistic = statistic + (tmp_statistic/n)^2
    }
  }
  for(i in seq2(1, d)){
    for (j  in seq2(i + 1, d)){
      for (r  in seq2(j + 1, d)){
        tmp_statistic = 0
        for (w in seq(n)){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')
          if (sample_norm >= rho){
            tmp_statistic = tmp_statistic + func3deg4(sample/sample_norm, i, j, r)
          }
        }
        statistic = statistic + (tmp_statistic/n)^2
      }
    }
  }
  for(r in seq2(2, d)){
    for (s  in seq2(r + 1, d)){
      tmp_statistic = 0
      for (w in seq(n)){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')
        if (sample_norm >= rho){
          tmp_statistic = tmp_statistic + func4deg4(sample/sample_norm, r, s)
        }
      }
      statistic = statistic + (tmp_statistic/n)^2
    }
  }
  for(i in seq2(1, d)){
    for (j  in seq2(i + 1, d)){
      for (r  in seq2(j + 1, d)){
        for (s  in seq2(r + 1, d)){
          tmp_statistic = 0
          for (w in seq(n)){
            sample = Z[w,]
            sample_norm = norm(sample, type='2')
            if (sample_norm >= rho){
              tmp_statistic = tmp_statistic + func5deg4(sample/sample_norm, i, j, r, s)
            }
          }
          statistic = statistic + (tmp_statistic/n)^2

        }
      }
    }
  }
  for (k  in seq2(2, d)){
    for (r  in seq2(k + 1, d)){
      for (s  in seq2(r + 1, d)){
        tmp_statistic = 0
        for (w in seq(n)){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')
          if (sample_norm >= rho){
            tmp_statistic = tmp_statistic + func6deg4(sample/sample_norm, k, r, s)
          }
        }
        statistic = statistic + (tmp_statistic/n)^2

      }
    }
  }
  for (j  in seq2(1, d)){
    for (r  in seq2(j + 1, d)){
      for (s  in seq2(r + 1, d)){
        tmp_statistic = 0
        for (w in seq(n)){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')
          if (sample_norm >= rho){
            tmp_statistic = tmp_statistic + func7deg4(sample/sample_norm, j, r, s)
          }
        }
        statistic = statistic + (tmp_statistic/n)^2
      }
    }
  }
  for(r in seq2(2, d)){
    for (s  in seq2(r + 1, d)){
      tmp_statistic = 0
      for (w in seq(n)){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')
        if (sample_norm > rho){
          tmp_statistic = tmp_statistic + func8deg4(sample/sample_norm, r, s)
        }
      }
      statistic = statistic + (tmp_statistic/n)^2
    }
   }
  statistic= n*statistic
  #the following three lines calculate the degree of freedom of the chi-squared distribution
  deg3 = choose(d + 2, 3) - d
  deg4 = choose(d + 3, 4) - choose(d + 1, 2)
  df = deg3 + deg4
  statistic = statistic[[1]]
  p_val = 1 - stats::pchisq(statistic/(1-epsilon), df = df)[[1]]  #p value

  names(statistic) = 'statistic'
  res <- list(method = 'Test for elliptical symmetry by Manzzoti et al.',
              data.name = dname,
              statistic = statistic,
              p.value = p_val,
              alternative = 'the distribution is not elliptically symmetric')
  class(res) <- "htest"
  return(res)
}


