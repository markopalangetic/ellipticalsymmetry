#In this file, the orthonormal basis of spherical harmonics is implemented.
#All functions from the basis up to degree 4 are implemented here.
#The functions can be found in Tables 1 and 2 in [1].
#In every method, 'x' stands for a function argument while 'i','j','k','r','s' are different parameters. The parameters have the same notation as in [1].

## [1]Manzotti, A., & Quiroz, A. J. (2001). Spherical harmonics in quadratic forms for testing multivariate normality. Test, 10(1), 87-104.



func1deg1 = function(x, i) {
  d = length(x)
  c1 = d
  const = sqrt(c1)
  res = x[i]
  return (const * res)
}

func1deg2 = function(x, i, j) {
  d = length(x)
  c1 = d * (d + 2)
  const = sqrt(c1)
  res = x[i] * x[j]
  return (const * res)
}

func2deg2 = function(x, k) {
  d = length(x)
  c1 = d * (d + 2)
  S = sum(x[1:k - 1] ^ 2)
  const = sqrt(c1 / (2 * k * (k - 1)))
  res = S - (k - 1) * x[k] ^ 2
  return (const * res)
}


func1deg3 = function(x, i, j, k) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4)
  const = sqrt(c1)
  res = x[i] * x[j] * x[k]
  return (const * res)
}

func2deg3 = function(x, k, r) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4)
  S = sum(x[1:k - 1] ^ 2)
  const = sqrt(c1 / (2 * k * (k - 1)))
  res = x[r] * (S - (k - 1) * x[k] ^ 2)
  return(const * res)
}

func3deg3 = function(x, j, k) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4)
  S = sum(x[1:k - 1] ^ 2)
  const = sqrt(c1 / (2 * (k + 1) * (k + 2)))
  res = x[j] * (S - (k + 1) * x[k] ^ 2)
  return(const * res)
}

func4deg3 = function(x, r) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4)
  S = sum(x[1:r - 1] ^ 2)
  const = sqrt(c1 / (6 * (r - 1) * (r + 2)))
  res = (r - 1) * x[r] ^ 3 - 3 * x[r] * S
  return(const * res)
}

func1deg4 = function(x, r) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(x[1:r - 1] ^ 2)
  const = sqrt(c1 / (24 * (r - 1) * (r + 1) * (r + 2) * (r + 4)))
  res = 6 * (r + 1) * x[r] ^ 2 * S - (r ^ 2 - 1) * x[r] ^ 4 - 3 * S ^
    2
  return(const * res)
}

func2deg4 = function(x, i, r) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(x[1:r - 1] ^ 2)
  const = sqrt(c1 / (6 * (r + 1) * (r + 4)))
  res = x[r] * x[i] * (3 * S - (r + 1) * x[r] ^ 2)
  return(const * res)
}

func3deg4 = function(x, i, j, r) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(x[1:r - 1] ^ 2)
  const = sqrt(c1 / (2 * (r + 3) * (r + 4)))
  res = x[i] * x[j] * (S - (r + 3) * x[r] ^ 2)
  return(const * res)
}
func4deg4 = function(x, r, s) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  Sr = sum(x[1:r - 1] ^ 2)
  Ss = sum(x[1:s - 1] ^ 2)
  const = sqrt(c1 / (4 * r * (r - 1) * (s + 3) * (s + 4)))
  res = (Ss - (s + 3) * x[s] ^ 2) * (Sr - (r - 1) * x[r] ^ 2)
  return(const * res)
}

func5deg4 = function(x, i, j, r, s) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  const = sqrt(c1)
  res = x[i] * x[j] * x[r] * x[s]
  return (const * res)
}

func6deg4 = function(x, k, r, s) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(x[1:k - 1] ^ 2)
  const = sqrt(c1 / (2 * k * (k - 1)))
  res = x[r] * x[s] * (S - (k - 1) * x[k] ^ 2)
  return(const * res)
}

func7deg4 = function(x, j, r, s) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(x[1:r - 1] ^ 2)
  const = sqrt(c1 / (2 * (r + 1) * (r + 2)))
  res = x[j] * x[s] * (S - (r + 1) * x[r] ^ 2)
  return(const * res)
}

func8deg4 = function(x, r, s) {
  d = length(x)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(x[1:r - 1] ^ 2)
  const = sqrt(c1 / (6 * (r - 1) * (r + 2)))
  res = x[r] * x[s] * (3 * S - (r - 1) * x[r] ^ 2)
  return(const * res)
}
